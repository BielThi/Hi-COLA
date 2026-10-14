#!/usr/bin/env python3
"""Write the HMF and the HODDIES halo array from a halo catalog made by run_halo_catalog.py.

The HMF mass grid is chosen with --mass-grid:
    reference        the masses of the target's reference HMF at the same redshift (Box1, LCDM),
                     so both can be compared bin by bin; fails if there is no reference HMF
    MIN MAX N        N log-spaced masses from MIN to MAX [Msun/h], e.g. 1e11 1e16 60
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

import hicola_paths as hp
from halo_pipeline import (
    VELOCITY_DEFINITION,
    compute_hmf,
    halo_catalog_to_hoddies,
    load_catalog_metadata,
    require_peculiar_velocities,
    save_hmf_table,
    validate_halo_catalog,
)


def default_hoddies_output(halo_catalog: Path) -> Path:
    return halo_catalog.with_name("halos200c_hoddies.npy")


def default_hmf_output(halo_catalog: Path) -> Path:
    return halo_catalog.with_name("halos200c_hmf.txt")


def mass_grid(spec: list[str], target: str, redshift: float) -> tuple[np.ndarray, str]:
    """Masses of the HMF and a description of where they come from (see the module docstring)."""
    if spec == ["reference"]:
        reference = hp.reference_file(target, "LCDM", 1, redshift, "HMF")
        if not reference.exists():
            raise SystemExit(f"--mass-grid reference: there is no reference HMF {reference}")
        return np.loadtxt(reference, comments="#")[:, 0], str(reference)
    if len(spec) == 3:
        m_min, m_max, n = float(spec[0]), float(spec[1]), int(spec[2])
        return np.logspace(np.log10(m_min), np.log10(m_max), n), f"{n} log-spaced masses from {m_min:g} to {m_max:g} Msun/h"
    raise SystemExit(f"--mass-grid must be 'reference' or 'MIN MAX N', got {' '.join(spec)!r}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Write HMF and HODDIES outputs from a halo catalog.")
    parser.add_argument("--target", required=True, help="e.g. COLA_NM4096, Abacus_Small")
    parser.add_argument("--model", required=True, help="e.g. LCDM, nDGP_N1.2")
    parser.add_argument("--realization", required=True, help="e.g. s1234, s1234-fixed")
    parser.add_argument("--redshift", type=float, required=True)
    parser.add_argument("--halo-catalog", type=Path, default=None, help="Default: from the run")
    parser.add_argument("--boxsize", type=float, required=True, help="Mpc/h")
    parser.add_argument("--mass-grid", nargs="+", required=True, metavar="reference | MIN MAX N",
                        help="HMF masses, see the module docstring")
    parser.add_argument("--hmf-output", type=Path, default=None)
    parser.add_argument("--hoddies-output", type=Path, default=None)
    args = parser.parse_args()

    run = (args.target, args.model, args.realization, args.redshift)
    args.halo_catalog = args.halo_catalog or hp.halo_catalog_file(*run)
    return args


def main() -> int:
    args = parse_args()
    catalog = np.load(args.halo_catalog)
    require_peculiar_velocities(load_catalog_metadata(args.halo_catalog), args.halo_catalog)
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
                "velocity_definition": VELOCITY_DEFINITION,
                "n_halos": int(len(hoddies)),
            },
            handle,
            indent=2,
            sort_keys=True,
        )

    masses, mass_grid_source = mass_grid(args.mass_grid, args.target, args.redshift)
    print(f"HMF mass grid: {mass_grid_source}")
    hmf = compute_hmf(catalog["m200c"], volume=args.boxsize**3, mass_grid=masses)
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
