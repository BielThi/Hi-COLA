#!/usr/bin/env python3
"""Make a HODDIES galaxy mock from the halo array of a run and write it in the EMC FITS format.

Input:  halo_catalog/halos200c_hoddies.npy of the run (written by jobs/halo_finder.sh).
Output: halo_catalog/hod/<name>.fits, the format of the EMC AbacusHOD mocks
        (see save_galaxy_fits in halo_pipeline.py): one table with X Y Z [Mpc/h, in
        [-L/2, L/2), real space], VX VY VZ [km/s], MASS [Msun/h] and ID of the host halo,
        IS_CENT; centrals first; header N_CENT, GAL_TYPE, the HOD parameters of the tracer
        and the settings of the mock.

Velocities: the halo arrays carry Gadget velocities, v_pec / sqrt(a) (Hi-COLA writes its snapshots
in that convention and nbodykit reads them unchanged), while the EMC format needs peculiar
velocities in km/s. The halo velocities (vx, vy, vz) and dispersions (Vrms) are therefore
multiplied by sqrt(a) before the mock is made, so the satellites get consistent velocities too.
If the halo finder is changed to write peculiar velocities, remove this conversion
(to_peculiar_velocities).

The HOD model and its parameters come from the HODDIES parameter file. If it sets the tracer's
`density`, HODDIES rescales Ac (and the satellites with it) to reach that number density.
The mock depends on the seed and on the number of threads (HODDIES seeds each block of halos).

    python scripts/run_hod_mock.py --target Abacus_Small --model LCDM --realization s1234-fixed \
        --redshift 0.5 --boxsize 500 --hod-param-file Notebooks/HOD_parameter_file_z0.500.yaml \
        --tracer LRG --seed 42 --name LRG_seed42 --nthreads 32
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

import hicola_paths as hp
from halo_pipeline import load_catalog_metadata, save_galaxy_fits


def import_hoddies():
    """HOD class of the HODDIES clone next to this repository (not the copy installed in the env)."""
    sys.path.insert(0, str(hp.HODDIES_ROOT))
    import HODDIES
    from HODDIES import HOD

    print(f"HODDIES:        {HODDIES.__file__}")
    return HOD


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--target", required=True, help="e.g. COLA_NM4096, Abacus_Small")
    parser.add_argument("--model", required=True, help="e.g. LCDM, nDGP_N1.2")
    parser.add_argument("--realization", required=True, help="e.g. s1234, s1234-fixed")
    parser.add_argument("--redshift", type=float, required=True)
    parser.add_argument("--boxsize", type=float, required=True, help="Mpc/h")
    parser.add_argument("--hod-param-file", type=Path, required=True, help="HODDIES parameter file (.yaml)")
    parser.add_argument("--tracer", required=True, help="tracer of the parameter file, e.g. LRG")
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--name", required=True, help="output file name without .fits, e.g. LRG_seed42")
    parser.add_argument("--nthreads", type=int, default=32, help="HODDIES threads (the mock depends on it)")
    parser.add_argument("--halo-array", type=Path, default=None, help="Default: halos200c_hoddies.npy of the run")
    parser.add_argument("--output", type=Path, default=None, help="Default: halo_catalog/hod/<name>.fits of the run")
    args = parser.parse_args()

    run = (args.target, args.model, args.realization, args.redshift)
    args.halo_array = args.halo_array or hp.hoddies_catalog_file(*run)
    args.output = args.output or hp.hod_mock_file(*run, args.name)
    return args


def to_peculiar_velocities(halos: np.ndarray, redshift: float) -> float:
    """Convert the Gadget velocities of a halo array (v_pec / sqrt(a)) to peculiar km/s, in place."""
    factor = np.sqrt(1.0 / (1.0 + redshift))
    for field in ("vx", "vy", "vz", "Vrms"):
        halos[field] = halos[field] * factor
    return factor


def main() -> int:
    args = parse_args()
    HOD = import_hoddies()

    halos = np.load(args.halo_array)
    factor = to_peculiar_velocities(halos, args.redshift)
    print(f"Velocities:     halo vx, vy, vz, Vrms multiplied by sqrt(a) = {factor:.4f} (Gadget -> peculiar km/s)")
    metadata = load_catalog_metadata(hp.halo_catalog_file(args.target, args.model, args.realization, args.redshift))
    if "boxsize" in metadata and not np.isclose(float(metadata["boxsize"]), args.boxsize):
        raise SystemExit(f"--boxsize {args.boxsize} differs from the halo catalog's {metadata['boxsize']}")
    print(f"Halos:          {args.halo_array} ({len(halos):,} halos)")

    hod = HOD(
        hcat_file=halos,
        param_file=str(args.hod_param_file),
        hcat={"boxsize": args.boxsize, "z_simu": args.redshift},
        nthreads=args.nthreads,
        seed=args.seed,
    )
    if args.tracer not in hod.args["tracers"]:
        raise SystemExit(f"Tracer {args.tracer} is not in {args.hod_param_file} (tracers: {hod.args['tracers']})")
    hod_parameters = hod.args[args.tracer]
    print(f"HOD parameters: {json.dumps(hod_parameters, default=str)}")

    mock = hod.make_mock_cat(tracers=args.tracer, fix_seed=args.seed)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    n_cent = save_galaxy_fits(
        mock,
        args.boxsize,
        args.tracer,
        hod_parameters,
        args.output,
        extra_header={
            "SIM_NAME": hp.simulation_name(args.target, args.model, args.realization),
            "REDSHIFT": args.redshift,
            "BOXSIZE": args.boxsize,
            "SEED": args.seed,
            "NTHREADS": hod.args["nthreads"],
            "HOD_FILE": args.hod_param_file.name,
            "VEL_DEF": "peculiar km/s",
        },
    )
    n_gal = len(mock["x"])
    print(f"Galaxies:       {n_gal:,} ({n_cent:,} centrals, satellite fraction {1 - n_cent / n_gal:.3f}), "
          f"n = {n_gal / args.boxsize**3:.3e} (h/Mpc)^3")
    print(f"Saved:          {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
