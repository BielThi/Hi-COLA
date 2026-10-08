#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

PROJECT_DIR = Path(__file__).resolve().parents[1]
if str(PROJECT_DIR) not in sys.path:
    sys.path.insert(0, str(PROJECT_DIR))

from scripts.halo_pipeline import compute_hmf, halo_catalog_to_hoddies, save_hmf_table, validate_halo_catalog


def default_hoddies_output(halo_catalog: Path) -> Path:
    return halo_catalog.with_name("halos200c_hoddies.npy")


def default_hmf_output(halo_catalog: Path) -> Path:
    return halo_catalog.with_name("halos200c_hmf.txt")


def default_mass_grid(project_dir: Path, snapshot_index: int) -> tuple[np.ndarray, str]:
    reference = project_dir / "Comparison" / "Box1" / f"Bart_LCDM_hmf_{snapshot_index}.txt"
    if reference.exists():
        return np.loadtxt(reference, comments="#")[:, 0], str(reference)
    return np.logspace(11.0, 16.0, 60), "np.logspace(11, 16, 60)"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Write HMF and HODDIES outputs from a halo catalog.")
    parser.add_argument("--halo-catalog", required=True, type=Path)
    parser.add_argument("--boxsize", type=float, default=1024.0)
    parser.add_argument("--redshift", type=float, default=1.650)
    parser.add_argument("--model", default="LCDM")
    parser.add_argument("--project-dir", type=Path, default=Path("/home/biel/Hi-COLA"))
    parser.add_argument("--snapshot-index", type=int, default=None)
    parser.add_argument("--hmf-output", type=Path, default=None)
    parser.add_argument("--hoddies-output", type=Path, default=None)
    return parser.parse_args()


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

    snapshot_index_by_redshift = {1.650: 0, 1.400: 1, 1.200: 2, 1.000: 3}
    snapshot_index = args.snapshot_index
    if snapshot_index is None:
        snapshot_index = snapshot_index_by_redshift.get(round(args.redshift, 3), 0)

    mass_grid, mass_grid_source = default_mass_grid(args.project_dir, snapshot_index)
    hmf = compute_hmf(catalog["m200c"], volume=args.boxsize**3, mass_grid=mass_grid)
    hmf_output = args.hmf_output or default_hmf_output(args.halo_catalog)
    save_hmf_table(
        hmf,
        hmf_output,
        {
            "label": f"Hi-COLA {args.model}",
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
