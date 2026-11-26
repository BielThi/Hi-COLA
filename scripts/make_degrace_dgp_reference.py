"""Reference matter P(k) files of the DEGRACE-DGP boxes (MG-GLAM LCDM and nDGP N1) for the pipeline.

Input: the PowerDM.R.<step>.<box>.dat files of DEGRACE-DGP (raw GLAM estimate: CIC of all particles
on an Ngrid^3 grid, no window deconvolution, no interlacing, no shot-noise subtraction; linear bins
of width k_f = 2 pi / L, each mode in the bin of the nearest integer multiple of k_f; k = mean |k| of
the modes of the bin; N = number of modes, printed with 5 significant digits).

Output, per box, model and step (hicola_data/README.md layout):
    <ref-root>/DEGRACE_DGP/<model>/Box<n>/DEGRACE_DGP_<model>_Box<n>_z<z>_PS.txt (+ .json)
    columns  k [h/Mpc]  P [(Mpc/h)^3]  N_modes (exact count)  P_raw [(Mpc/h)^3]
    P = P_raw / <W^2>_bin - 1/nbar,  W(k) = prod_i sinc^2(pi n_i / Ngrid)  (CIC window, n_i = k_i / k_f),
    <W^2>_bin = mean of W^2 over the modes of the bin, computed exactly on the full Ngrid^3 grid
    (n_i in [-Ngrid/2, Ngrid/2 - 1]) with the same binning; aliasing is not corrected.
    model: LCDM -> LCDM, N1 -> nDGP_N1.0; z: snapnum_to_z.txt with 3 decimals; Box<n> = boxNNNN.
The script checks its mode counts against the N column of the files (to their 5 digits) and its
mean |k| against the k column, and refuses to overwrite existing files.

    python scripts/make_degrace_dgp_reference.py --raw ~/hicola_data/DEGRACE-DGP --ref-root ~/hicola_data/reference \
        --boxsize 1024 --npart-1d 2048 --ngrid 4096 --boxes 1 100 --threads 8
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numba
import numpy as np

MODELS = {"LCDM": "LCDM", "N1": "nDGP_N1.0"}
SIM = "DEGRACE_DGP"


@numba.njit(parallel=True, cache=False)
def _bin_window(ng):
    """Per bin b = 1..ng/2: number of modes, sum of W^2, sum of |n|, over n_i in [-ng/2, ng/2-1]."""
    h = ng // 2
    nb = h
    acc = np.zeros((h + 1, 3, nb + 1))
    s = np.empty(h + 1)
    for n in range(h + 1):
        x = np.pi * n / ng
        s[n] = 1.0 if n == 0 else (np.sin(x) / x) ** 4          # sinc^4 per axis = (sinc^2)^2
    for nx in numba.prange(h + 1):
        wx = 1.0 if (nx == 0 or nx == h) else 2.0               # n and -n, except 0 and -ng/2
        for ny in range(h + 1):
            wy = 1.0 if (ny == 0 or ny == h) else 2.0
            for nz in range(h + 1):
                wz = 1.0 if (nz == 0 or nz == h) else 2.0
                r = np.sqrt(nx * nx + ny * ny + nz * nz)
                b = int(np.rint(r))
                if b < 1 or b > nb:
                    continue
                w = wx * wy * wz
                acc[nx, 0, b] += w
                acc[nx, 1, b] += w * s[nx] * s[ny] * s[nz]
                acc[nx, 2, b] += w * r
    return acc.sum(axis=0)


def read_powerdm(path: Path):
    lines = path.read_text().splitlines()
    d = np.loadtxt(lines[3:])
    return lines[:3], d[:, 1], d[:, 2], d[:, 3]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--raw", required=True, help="DEGRACE-DGP folder with boxNNNN/{LCDM,N1}/PowerDM.R.*.dat and DEGRACE-DGP-common/")
    parser.add_argument("--ref-root", required=True, help="hicola_data/reference")
    parser.add_argument("--boxsize", type=float, required=True, help="box size [Mpc/h]")
    parser.add_argument("--npart-1d", type=int, required=True, help="particles per dimension (for 1/nbar)")
    parser.add_argument("--ngrid", type=int, required=True, help="CIC grid of the PowerDM estimate")
    parser.add_argument("--boxes", type=int, nargs=2, required=True, metavar=("FIRST", "LAST"))
    parser.add_argument("--threads", type=int, required=True)
    args = parser.parse_args()
    numba.set_num_threads(args.threads)

    raw = Path(args.raw).expanduser()
    out_root = Path(args.ref_root).expanduser() / SIM
    kf = 2 * np.pi / args.boxsize
    shot = args.boxsize**3 / args.npart_1d**3
    snap = np.loadtxt(raw / "DEGRACE-DGP-common" / "snapnum_to_z.txt", skiprows=1)
    step_z = {int(s): z for s, a, z, _ in snap}

    print(f"bin-averaged CIC window on {args.ngrid}^3 (numba, {args.threads} threads) ...", flush=True)
    cnt, w2sum, rsum = _bin_window(args.ngrid)[:, 1:]
    w2 = w2sum / cnt
    kmean = rsum / cnt * kf
    print(f"<W^2>: bin 1 {w2[0]:.6f}, k=1 h/Mpc {np.interp(1.0, kmean, w2):.5f}, k=5 {np.interp(5.0, kmean, w2):.4f}, last bin {w2[-1]:.4f}")

    boxes = range(args.boxes[0], args.boxes[1] + 1)
    first = raw / f"box{boxes[0]:04d}" / "LCDM"
    _, k0, n0, _ = read_powerdm(next(first.glob("PowerDM.R.*.dat")))
    reln = np.abs(n0 / cnt - 1)
    if reln.max() > 5.01e-5 or np.abs(kmean - k0).max() > 1.1e-5:
        raise SystemExit(f"mode counts or mean k differ from the files: max rel dN {reln.max():.2e}, max dk {np.abs(kmean-k0).max():.2e}")
    print(f"mode counts match the files (max rel. diff {reln.max():.1e}, file precision 5 digits); mean |k| matches (max diff {np.abs(kmean-k0).max():.1e})")

    nwritten = 0
    for b in boxes:
        for gmodel, model in MODELS.items():
            for f in sorted((raw / f"box{b:04d}" / gmodel).glob("PowerDM.R.*.dat")):
                step = int(f.name.split(".")[2])
                z = step_z[step]
                hdr, k, n, praw = read_powerdm(f)
                if not (np.array_equal(k, k0) and np.array_equal(n, n0)):
                    raise SystemExit(f"{f}: k or N columns differ from {first}")
                p = praw / w2 - shot
                odir = out_root / model / f"Box{b}"
                odir.mkdir(parents=True, exist_ok=True)
                out = odir / f"{SIM}_{model}_Box{b}_z{z:.3f}_PS.txt"
                if out.exists():
                    raise SystemExit(f"{out} exists: not overwriting reference files")
                header = (f"DEGRACE-DGP (MG-GLAM) {model} box{b:04d} step {step} z = {z:.5f}: real-space matter P(k) from {f.name}\n"
                          f"P = P_raw / <W_CIC^2>_bin - 1/nbar (CIC on {args.ngrid}^3, window averaged over the modes of each bin, "
                          f"no aliasing correction, 1/nbar = {shot:g} (Mpc/h)^3); P_raw = GLAM's raw estimate\n"
                          f"k[h/Mpc]  P[(Mpc/h)^3]  N_modes (exact count; GLAM's file gives it to 5 digits)  P_raw[(Mpc/h)^3]")
                np.savetxt(out, np.column_stack([k, p, cnt, praw]), fmt=["%.5f", "%.6e", "%d", "%.5e"], header=header)
                meta = {"sim": SIM, "model": model, "glam_model": gmodel, "box": b, "step": step, "redshift": z,
                        "source": str(f), "source_header": hdr[:2], "boxsize_mpc_h": args.boxsize, "ngrid_cic": args.ngrid,
                        "npart_1d": args.npart_1d, "shot_noise_subtracted": shot,
                        "window": "CIC, prod_i sinc^2(pi n_i/Ngrid), <W^2> averaged over the modes of each bin, divided out",
                        "interlacing": False, "aliasing_corrected": False,
                        "made_by": "scripts/make_degrace_dgp_reference.py"}
                Path(str(out) + ".json").write_text(json.dumps(meta, indent=1) + "\n")
                nwritten += 1
    print(f"wrote {nwritten} PS files under {out_root}")


if __name__ == "__main__":
    main()
