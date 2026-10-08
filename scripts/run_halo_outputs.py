#!/usr/bin/env python3
"""Write the HMF and the HODDIES halo array from a halo catalog made by run_halo_catalog.py.

The mass grid of the HMF is the one of the target's reference HMF at the same redshift
(Box1, LCDM) when it exists, so both can be compared bin by bin.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

import hicola_paths as hp
from halo_pipeline import compute_hmf, halo_catalog_to_hoddies, save_hmf_table, validate_halo_catalog


def default_hoddies_output(halo_catalog: Path) -> Path:
    return halo_catalog.with_name("halos200c_hoddies.npy")


def default_hmf_output(halo_catalog: Path) -> Path:
    return halo_catalog.with_name("halos200c_hmf.txt")


def default_mass_grid(target: str, redshift: float) -> tuple[np.ndarray, str]:
    reference = hp.reference_file(target, "LCDM", 1, redshift, "HMF")
    if reference.exists():
        return np.loadtxt(reference, comments="#")[:, 0], str(reference)
    return np.logspace(11.0, 16.0, 60), "np.logspace(11, 16, 60)"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Write HMF and HODDIES outputs from a halo catalog.")
    parser.add_argument("--target", required=True, choices=sorted(hp.TARGETS))
    parser.add_argument("--model", required=True, help="e.g. LCDM, nDGP_N1.2")
    parser.add_argument("--realization", required=True, help="e.g. s1234, s1234-fixed")
    parser.add_argument("--redshift", type=float, required=True)
    parser.add_argument("--halo-catalog", type=Path, default=None, help="Default: from the run")
    parser.add_argument("--boxsize", type=float, default=None, help="Default: the target's box size")
    parser.add_argument("--hmf-output", type=Path, default=None)
    parser.add_argument("--hoddies-output", type=Path, default=None)
    args = parser.parse_args()

    run = (args.target, args.model, args.realization, args.redshift)
    args.halo_catalog = args.halo_catalog or hp.halo_catalog_file(*run)
    args.boxsize = args.boxsize or hp.TARGETS[args.target]["boxsize"]
    return args


def main() -> int:
    args = parse_args()
    catalog = np.load(args.halo_catalog)
    validate_halo_catalog(catalog, boxsize=args.boxsize)

    hoddies_output = args.hoddies_output or default_hoddies_output(args.halo_catalog)
    hoddies = halo_catalog_to_hoddies(catalog)
    np.save(hoddies_output, hoddies)
    with open(str(hoddies_output) + ".json", "w", encoding="utf-8") as handle:
        json.dump(
            {
                "source_catalog": str(args.halo_catalog),
                "output": str(hoddies_output),
                "schema": "HODDIES halo array",
                "mass_units": "Msun/h",
                "radius_units": "kpc/h",
                "position_units": "Mpc/h",
                "velocity_units": "km/s",
                "n_halos": int(len(hoddies)),
            },
            handle,
            indent=2,
            sort_keys=True,
        )

    mass_grid, mass_grid_source = default_mass_grid(args.target, args.redshift)
    hmf = compute_hmf(catalog["m200c"], volume=args.boxsize**3, mass_grid=mass_grid)
    hmf_output = args.hmf_output or default_hmf_output(args.halo_catalog)
    save_hmf_table(
        hmf,
        hmf_output,
        {
            "label": f"Hi-COLA {hp.simulation_name(args.target, args.model, args.realization)}",
            "target": args.target,
            "model": args.model,
            "realization": args.realization,
            "catalog": str(args.halo_catalog),
            "output": str(hmf_output),
            "mass_field": "m200c",
            "mass_units": "Msun/h",
            "volume_units": "(Mpc/h)^3",
            "volume": args.boxsize**3,
            "redshift": args.redshift,
            "mass_grid_source": mass_grid_source,
            "hmf_definition": "cumulative n(>M) and differential dn/dlog10M from the same mass grid",
            "poisson_errors": "sqrt(N) / volume, or sqrt(N_bin) / (volume dlog10M)",
            "n_halos": int(hmf["n_halos"]),
        },
    )

    print(f"Saved halo catalog: {args.halo_catalog}")
    print(f"Saved HODDIES catalog: {hoddies_output}")
    print(f"Saved HMF: {hmf_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
