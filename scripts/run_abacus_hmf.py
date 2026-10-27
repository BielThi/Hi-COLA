#!/usr/bin/env python3
"""Halo mass function of the Abacus_Small reference box, from its CompaSO halo catalog.

Run by jobs/hmf_abacus.sh. Masses are CompaSO masses: N (the number of particles of the L1
halo, a spherical overdensity of about 200 times the mean density) x the particle mass. With
--cleaned 1 the merger-tree cleaned catalog is used (the one of the DESI HOD fits: halos merged
into others have N = 0 and are dropped). The HMF is computed as for the Hi-COLA runs
(halo_pipeline.compute_hmf: cumulative n(>M) and dn/dlog10M with Poisson errors) on the masses
of --mass-grid, so the two compare bin by bin.

Output: $HICOLA_DATA_ROOT/reference/Abacus_Small/<model>/Box1/Abacus_Small_<model>_Box1_z<z>_HMF.txt
(+ .json), columns as the Hi-COLA HMFs with the mass column M_CompaSO.

    python scripts/run_abacus_hmf.py --redshift 0.5 --model LCDM --cleaned 1 --mass-grid 1e11 1e16 60
"""
from __future__ import annotations

import argparse
import time

import numpy as np

import hicola_paths as hp
from halo_pipeline import compute_hmf, save_hmf_table


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--redshift", type=float, required=True, help="redshift of the CompaSO catalog, e.g. 0.5")
    parser.add_argument("--model", required=True, help="model name of the output file; the Abacus_Small box is LCDM")
    parser.add_argument("--cleaned", type=int, choices=(0, 1), required=True, help="use the cleaned catalog")
    parser.add_argument("--mass-grid", type=float, nargs=3, required=True, metavar=("MIN", "MAX", "N"),
                        help="N log-spaced masses from MIN to MAX [Msun/h], as HMF_MASS_GRID of the halo finder")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    from abacusnbody.data.compaso_halo_catalog import CompaSOHaloCatalog

    path = hp.abacus_halo_dir(args.redshift)
    start = time.time()
    catalog = CompaSOHaloCatalog(str(path), fields=["N"], cleaned=bool(args.cleaned))
    header = catalog.header
    n_particles = np.asarray(catalog.halos["N"], dtype="i8")
    kept = n_particles > 0  # the cleaned catalog keeps merged halos with N = 0
    particle_mass = float(header["ParticleMassHMsun"])
    masses = n_particles[kept] * particle_mass
    boxsize = float(header["BoxSize"])
    print(f"Loaded {path} in {time.time() - start:.0f} s: {len(n_particles):,} halos, {kept.sum():,} with N > 0 "
          f"(cleaned {bool(args.cleaned)}), particle mass {particle_mass:.4e} Msun/h, box {boxsize} Mpc/h, "
          f"N >= {n_particles[kept].min()} (M >= {masses.min():.3e} Msun/h)")

    m_min, m_max, n = args.mass_grid
    grid = np.logspace(np.log10(m_min), np.log10(m_max), int(n))
    hmf = compute_hmf(masses, volume=boxsize**3, mass_grid=grid)
    output = hp.reference_file("Abacus_Small", args.model, 1, args.redshift, "HMF")
    save_hmf_table(
        hmf,
        output,
        {
            "label": "AbacusSummit small_c000_ph3000 (CompaSO)",
            "source": str(path),
            "cleaned": bool(args.cleaned),
            "mass_definition": "CompaSO: N (L1 halo, SO ~200 x mean density) x particle mass"
                               + ("; merger-tree cleaned, N = 0 halos dropped" if args.cleaned else ""),
            "particle_mass_msun_h": particle_mass,
            "min_particles": int(n_particles[kept].min()),
            "min_mass_msun_h": float(masses.min()),
            "mass_units": "Msun/h",
            "volume_units": "(Mpc/h)^3",
            "volume": boxsize**3,
            "redshift": float(header["Redshift"]),
            "mass_grid_source": f"{int(n)} log-spaced masses from {m_min:g} to {m_max:g} Msun/h",
            "hmf_definition": "cumulative n(>M) and differential dn/dlog10M from the same mass grid",
            "poisson_errors": "sqrt(N) / volume, or sqrt(N_bin) / (volume dlog10M)",
            "n_halos": int(hmf["n_halos"]),
        },
        mass_column="M_CompaSO",
    )
    print(f"Saved {output} (+ .json)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
