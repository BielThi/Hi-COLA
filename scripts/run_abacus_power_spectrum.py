#!/usr/bin/env python3
"""Compute the Abacus matter power spectrum from a CompaSOHaloCatalog subsample.

Intended to be launched as a batch job (see jobs/power_spectrum_abacus.sh) rather than run
inline in a notebook: loading the halo catalog and painting the full particle subsample
onto a mesh needs several GB, more than is comfortably available on a shared login-node
session.

CompaSOHaloCatalog(..., subsamples=dict(A=True)) only loads the halo-bound part of
Abacus's fair particle subsample; the field (unbound) part of that same subsample lives
in a sibling `field_rv_<letter>/` directory and has to be read separately with
abacusnbody.data.read_abacus.read_asdf. Both halo-bound and field particles are combined
here so the result is a genuine matter power spectrum, not a halo-tracer spectrum.

Particles are streamed through the mesh one source file at a time (rather than
concatenated into a single array first) so peak memory stays roughly independent of the
total particle count -- this matters once subsample B is enabled too (3% + 7% = 10%) or
bigger boxes are used.
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import sys
import time

import numpy as np

import hicola_paths as hp


def log(message: str) -> None:
    elapsed = time.strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{elapsed}] {message}", flush=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--redshift", type=float, default=0.5, help="Used for the default --path-to-sim")
    parser.add_argument(
        "--path-to-sim",
        default=None,
        help="CompaSOHaloCatalog redshift directory (default: the Abacus_Small raw data at --redshift)",
    )
    parser.add_argument("--load-A", type=int, choices=(0, 1), default=1, help="Load subsample A (3%%)")
    parser.add_argument("--load-B", type=int, choices=(0, 1), default=0, help="Load subsample B (7%%)")
    parser.add_argument("--nmesh", type=int, default=512)
    parser.add_argument(
        "--kbins",
        type=int,
        default=60,
        help=(
            "Number of output k-bins, independent of --nmesh. A 500 Mpc/h-scale box has "
            "few independent Fourier modes near k_min, so using nmesh as the bin count "
            "(the old default) packs hundreds of bins into a regime with only a handful "
            "of modes each, giving mostly-empty, noisy low-k bins."
        ),
    )
    parser.add_argument("--nthread", type=int, default=None, help="Default: all available cores")
    parser.add_argument("--model", default="LCDM")
    parser.add_argument(
        "--output",
        default=None,
        help=(
            "Output file (and .json metadata next to it). Default: the Abacus_Small reference "
            "file Abacus_Small_<model>_Box1_z<redshift>_PS.txt"
        ),
    )
    args = parser.parse_args()
    args.path_to_sim = args.path_to_sim or str(hp.abacus_halo_dir(args.redshift))
    return args


def main() -> int:
    args = parse_args()
    nthread = args.nthread or os.cpu_count()

    t0 = time.time()
    from abacusnbody.data.compaso_halo_catalog import CompaSOHaloCatalog
    from abacusnbody.data.read_abacus import read_asdf
    import abacusnbody.analysis.power_spectrum as ps

    log(f"Loading halo catalog + subsamples from {args.path_to_sim}")
    hcat = CompaSOHaloCatalog(
        args.path_to_sim,
        fields=["N"],
        subsamples=dict(A=bool(args.load_A), B=bool(args.load_B)),
        cleaned=True,
    )
    header = hcat.header
    Lbox = header["BoxSize"]
    loaded_letters = [l for l, on in (("A", args.load_A), ("B", args.load_B)) if on]
    log(f"Loaded catalog. Subsamples used: {loaded_letters}")

    field_files = [
        fn
        for letter in loaded_letters
        for fn in sorted(glob.glob(os.path.join(args.path_to_sim, f"field_rv_{letter}", "*.asdf")))
    ]

    nmesh = args.nmesh
    paste = "TSC"
    d = Lbox / nmesh
    field = np.zeros((nmesh, nmesh, nmesh), dtype=np.float32)
    field_shift = np.zeros((nmesh, nmesh, nmesh), dtype=np.float32)
    total_n = 0

    def paint_chunk(chunk: np.ndarray, label: str) -> None:
        nonlocal total_n
        t = time.time()
        ps.tsc_parallel(chunk, field, Lbox, nthread=nthread, offset=0.0)
        ps.tsc_parallel(chunk, field_shift, Lbox, nthread=nthread, offset=0.5 * d)
        total_n += len(chunk)
        log(f"painted {label}: n={len(chunk):,} ({time.time() - t:.1f}s)")

    halo_pos = np.array(hcat.subsamples["pos"])
    paint_chunk(halo_pos, "halo-bound subsample")
    del halo_pos, hcat

    for fn in field_files:
        chunk = np.asarray(read_asdf(fn, load=["pos"])["pos"])
        paint_chunk(chunk, os.path.basename(fn))
        del chunk

    log(f"Painted {total_n:,} particles total ({time.time() - t0:.1f}s elapsed)")

    ps.normalize_field(field, tot_weight=total_n, inplace=True, nthread=nthread)
    ps.normalize_field(field_shift, tot_weight=total_n, inplace=True, nthread=nthread)

    field_fft = ps.rfftn(field, workers=nthread)
    del field
    field_shift_fft = ps.rfftn(field_shift, workers=nthread)
    del field_shift
    ps.shift_field_fft(field_fft, field_shift_fft, nmesh, Lbox, d)
    del field_shift_fft
    log("FFT + interlacing combine done")

    W = ps.get_W_compensated(Lbox, nmesh, paste, interlaced=True)
    field_fft /= (
        W[:, np.newaxis, np.newaxis]
        * W[np.newaxis, :, np.newaxis]
        * W[np.newaxis, np.newaxis, : (nmesh // 2 + 1)]
    )

    k_max = np.pi * nmesh / Lbox
    kbin_edges, mubins = ps.get_k_mu_edges(Lbox, k_max, args.kbins, 1, True)
    binned = ps.calc_pk_from_deltak(field_fft, Lbox, kbin_edges, mubins, nthread=nthread)
    k = (kbin_edges[1:] + kbin_edges[:-1]) * 0.5
    Nmodes = np.asarray(binned["N_mode"])

    # calc_power (and by extension this hand-rolled equivalent) returns the raw power,
    # which still includes shot noise; Bart's comparison files are shot-noise-corrected,
    # so subtract it here to match.
    Pk_shot = Lbox**3 / total_n
    Pk = np.asarray(binned["power"]) - Pk_shot
    log("Binning done")

    out_path = args.output or str(hp.reference_file("Abacus_Small", args.model, 1, header["Redshift"], "PS"))
    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    h = header["H0"] / 100.0
    header_lines = [
        "# Real space DM power spectrum.",
        "# The power spectrum is corrected for shot-noise and window function.",
        f"# Ngrid={nmesh}\tkbins={args.kbins}\tInterlacing was set to True.\tSubsamples used: {loaded_letters}",
        f"# redshift={header['Redshift']:.6f}\tmodel={args.model}\tH0={header['H0']}\t"
        f"Om={header['Omega_M']}\tOl={header['Omega_DE']}\tOb={header['omega_b'] / h**2:.6f}\t"
        f"sigma8=NA\tns={header['n_s']}\tTcmb=0\tN_ur={header['N_ur']}\tN_ncdm={header['N_ncdm']}",
        f"# n_dm={total_n:.6e}[h/Mpc]^3",
        "# k[h/Mpc]\t\tPk_dm[Mpc^3/h^3]\t\tNmodes",
    ]
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("\n".join(header_lines) + "\n")
        for ki, pi, ni in zip(k, Pk, Nmodes):
            f.write(f"{ki:.18e} {pi:.18e} {ni:.18e}\n")
    log(f"Wrote {out_path}")

    metadata = {
        "path_to_sim": args.path_to_sim,
        "output": out_path,
        "model": args.model,
        "redshift": header["Redshift"],
        "boxsize_mpc_h": Lbox,
        "subsamples_used": loaded_letters,
        "n_particles": int(total_n),
        "nmesh": nmesh,
        "kbins": args.kbins,
        "paste": paste,
        "interlaced": True,
        "compensated": True,
        "shot_noise_subtracted_mpc3_h3": Pk_shot,
        "cosmology": {
            "H0": header["H0"],
            "Omega_M": header["Omega_M"],
            "Omega_DE": header["Omega_DE"],
            "omega_b": header["omega_b"],
            "n_s": header["n_s"],
            "N_ur": header["N_ur"],
            "N_ncdm": header["N_ncdm"],
        },
        "elapsed_seconds": time.time() - t0,
    }
    with open(out_path + ".json", "w", encoding="utf-8") as f:
        json.dump(metadata, f, indent=2, sort_keys=True)
    log(f"Wrote {out_path}.json")

    return 0


if __name__ == "__main__":
    sys.exit(main())
