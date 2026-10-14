#!/usr/bin/env python3
"""Run the expensive nbodykit FOF halo-finding step.

This script is intended to be launched with MPI, normally through jobs/halo_finder.sh,
for example:

    srun -n 72 python scripts/run_halo_catalog.py --target Abacus_Small --model LCDM \
        --realization s1234 --redshift 0.5 --boxsize 500 --linking-length 0.2 --nmin 50 \
        --mdef 200c --h 0.6736 --omega-cdm 0.26589 --omega-b 0.049302 --omega-lambda 0.684808 \
        --tcmb 2.726 --neff 3.046 --n-s 0.9649

All halo settings and the cosmology (taken from the run file by jobs/halo_finder.sh) are given
explicitly; the snapshot and output file follow from the run identity (see hicola_paths.py)
unless --snapshot/--output are given.

It writes one flat structured NumPy halo catalog. Downstream notebooks can use
that catalog directly for HMFs, or convert it to HODDIES' expected field names
when running HOD models.
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import struct
import sys
import time

import numpy as np
from mpi4py import MPI

import hicola_paths as hp
from halo_pipeline import VELOCITY_DEFINITION


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--target", required=True, help="e.g. COLA_NM4096, Abacus_Small")
    parser.add_argument("--model", required=True, help="e.g. LCDM, nDGP_N1.2")
    parser.add_argument("--realization", required=True, help="e.g. s1234, s1234-fixed")
    parser.add_argument("--redshift", type=float, required=True)
    parser.add_argument("--snapshot", default=None, help="Gadget snapshot glob (default: from the run)")
    parser.add_argument("--output", default=None, help="Output halo catalog .npy file (default: from the run)")
    parser.add_argument("--boxsize", type=float, required=True, help="Mpc/h")
    parser.add_argument("--linking-length", type=float, required=True, help="in units of the mean particle separation")
    parser.add_argument("--nmin", type=int, required=True, help="minimum number of particles per halo")
    parser.add_argument("--mdef", required=True, choices=["200c"],
                        help="mass definition; only 200c, the catalog fields and files are named for it")
    parser.add_argument("--particle-type", type=int, default=1, help="Gadget particle type of the dark matter")
    cosmo = parser.add_argument_group("cosmology (of the simulation, used for the halo masses)")
    cosmo.add_argument("--h", type=float, required=True)
    cosmo.add_argument("--omega-cdm", type=float, required=True)
    cosmo.add_argument("--omega-b", type=float, required=True)
    cosmo.add_argument("--omega-mnu", type=float, default=0.0)
    cosmo.add_argument("--omega-lambda", type=float, required=True, help="recorded in the metadata (the cosmology is flat)")
    cosmo.add_argument("--tcmb", type=float, required=True, help="K")
    cosmo.add_argument("--neff", type=float, required=True)
    cosmo.add_argument("--n-s", type=float, required=True)
    args = parser.parse_args()

    run = (args.target, args.model, args.realization, args.redshift)
    args.snapshot = args.snapshot or hp.snapshot_glob(*run)
    args.output = args.output or str(hp.halo_catalog_file(*run))
    return args


def run_cosmology(args: argparse.Namespace):
    """Astropy cosmology of the simulation and its parameters as stored in the metadata."""
    from astropy import units as u
    from astropy.cosmology import FlatLambdaCDM

    omega_m = args.omega_cdm + args.omega_b + args.omega_mnu
    astropy_cosmology = FlatLambdaCDM(
        H0=100.0 * args.h, Om0=omega_m, Ob0=args.omega_b, Tcmb0=args.tcmb * u.K, Neff=args.neff
    )
    parameters = {
        "H0_km_s_Mpc": 100.0 * args.h,
        "h": args.h,
        "Omega_m": omega_m,
        "Omega_cdm": args.omega_cdm,
        "Omega_b": args.omega_b,
        "Omega_Lambda": args.omega_lambda,
        "Tcmb0_K": args.tcmb,
        "Neff": args.neff,
        "n_s": args.n_s,
    }
    return astropy_cosmology, parameters


def log(comm: MPI.Comm, message: str) -> None:
    if comm.rank == 0:
        elapsed = time.strftime("%Y-%m-%d %H:%M:%S")
        print(f"[{elapsed}] {message}", flush=True)


def read_particle_mass(snapshot_glob: str, particle_type: int) -> float:
    files = sorted(glob.glob(snapshot_glob))
    if not files:
        raise FileNotFoundError(f"No files matched {snapshot_glob!r}")

    with open(files[0], "rb") as handle:
        handle.read(4)
        header = handle.read(256)

    npart = struct.unpack("<6I", header[:24])
    massarr = struct.unpack("<6d", header[24:72])

    if particle_type < 0 or particle_type >= 6:
        raise ValueError("--particle-type must be between 0 and 5")
    if npart[particle_type] <= 0:
        raise ValueError(f"Particle type {particle_type} has zero particles in {files[0]}")
    if massarr[particle_type] <= 0:
        raise ValueError(
            f"Particle type {particle_type} has no constant Massarr entry; "
            "extend this script to read one Mass value from the catalog."
        )

    return float(massarr[particle_type])


def read_scale_factor(snapshot_glob: str) -> float:
    """Scale factor a stored in the Gadget header (field `time`)."""
    with open(sorted(glob.glob(snapshot_glob))[0], "rb") as handle:
        handle.read(4)
        header = handle.read(256)
    return float(struct.unpack("<d", header[72:80])[0])


def compute_local_halo_catalog_array(halos, vrms: np.ndarray, length: np.ndarray, boxsize: float) -> np.ndarray:
    position = halos["Position"].compute()
    velocity = halos["Velocity"].compute()
    mass = halos["Mass"].compute()
    radius = halos["Radius"].compute()
    concentration = halos["Concentration"].compute()

    dtype = [
        ("x", "f8"),
        ("y", "f8"),
        ("z", "f8"),
        ("vx", "f8"),
        ("vy", "f8"),
        ("vz", "f8"),
        ("m200c", "f8"),
        ("log10_m200c", "f8"),
        ("r200c", "f8"),
        ("rs", "f8"),
        ("c", "f8"),
        ("vrms", "f8"),
        ("length", "i8"),
        ("halo_id", "i8"),
    ]

    arr = np.zeros(len(mass), dtype=dtype)
    # nbodykit's FOF center-of-mass wraps positions into [0, boxsize) in
    # float64 but stores them as float32, so values within float32 precision
    # of the upper edge (e.g. 1023.9999998) round up to exactly `boxsize`.
    # Re-wrap after upcasting to float64 to keep positions strictly in range.
    arr["x"], arr["y"], arr["z"] = (position.T.astype("f8")) % boxsize
    arr["vx"], arr["vy"], arr["vz"] = velocity.T
    arr["m200c"] = mass
    arr["log10_m200c"] = np.log10(mass)
    arr["r200c"] = 1000.0 * radius
    arr["c"] = concentration
    arr["rs"] = arr["r200c"] / concentration
    arr["vrms"] = vrms
    arr["length"] = length
    return arr


def compute_local_fof_stats(fof, source, local_size: int) -> tuple[np.ndarray, np.ndarray]:
    """Compute per-halo FoF member counts and 3D RMS velocity."""
    comm = fof.comm
    labels = np.asarray(fof.labels)
    velocity = np.asarray(source.compute(source["Velocity"]), dtype="f8")

    nlabels = max(comm.allgather(int(labels.max()))) + 1
    count = np.bincount(labels, minlength=nlabels).astype("f8")
    comm.Allreduce(MPI.IN_PLACE, count, op=MPI.SUM)

    sum_v = np.vstack(
        [np.bincount(labels, weights=velocity[:, i], minlength=nlabels) for i in range(3)]
    ).T
    sum_v2 = np.vstack(
        [np.bincount(labels, weights=velocity[:, i] ** 2, minlength=nlabels) for i in range(3)]
    ).T
    comm.Allreduce(MPI.IN_PLACE, sum_v, op=MPI.SUM)
    comm.Allreduce(MPI.IN_PLACE, sum_v2, op=MPI.SUM)

    valid = count > 0
    mean_v = np.zeros_like(sum_v)
    mean_v2 = np.zeros_like(sum_v2)
    mean_v[valid] = sum_v[valid] / count[valid, None]
    mean_v2[valid] = sum_v2[valid] / count[valid, None]

    variance_3d = np.maximum(np.sum(mean_v2 - mean_v**2, axis=1), 0.0)
    vrms = np.sqrt(variance_3d)
    length = count.astype("i8")

    # nbodykit's halo catalog drops label 0 and scatters the remaining labels
    # evenly across ranks, so slice the global arrays with the local offset.
    local_start = sum(comm.allgather(local_size)[: comm.rank])
    local_slice = slice(local_start, local_start + local_size)
    return vrms[1:][local_slice], length[1:][local_slice]


def main() -> int:
    args = parse_args()
    comm = MPI.COMM_WORLD

    os.environ.setdefault("OMP_NUM_THREADS", "1")
    os.environ.setdefault("MKL_NUM_THREADS", "1")
    os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

    log(comm, f"Starting FOF halo finding on {comm.size} MPI ranks")
    log(comm, f"Snapshot: {args.snapshot}")

    particle_mass_code = read_particle_mass(args.snapshot, args.particle_type)
    scale_factor = read_scale_factor(args.snapshot)
    if not np.isclose(scale_factor, 1.0 / (1.0 + args.redshift), rtol=1e-3):
        raise SystemExit(f"--redshift {args.redshift} does not match the snapshot (a = {scale_factor:.5f})")
    log(comm, f"Velocities: Gadget velocities x sqrt(a) = {np.sqrt(scale_factor):.5f} -> peculiar km/s")
    particle_mass = particle_mass_code * 1.0e10
    log(comm, f"Particle mass from Gadget header: {particle_mass_code} x 1e10 Msun/h")

    from nbodykit.algorithms.fof import FOF
    from nbodykit.cosmology.cosmology import Cosmology
    from nbodykit.lab import Gadget1Catalog

    astropy_cosmology, cosmology_parameters = run_cosmology(args)
    cosmology = Cosmology.from_astropy(astropy_cosmology, n_s=args.n_s)
    log(comm, f"Using the simulation's cosmology: {cosmology_parameters}")

    log(comm, "Opening Gadget catalog")
    gadget = Gadget1Catalog(args.snapshot, comm=comm)
    # Gadget snapshots store v_pec / sqrt(a) (nbodykit's GadgetVelocity); convert to peculiar km/s
    gadget["Velocity"] = gadget["GadgetVelocity"] * np.sqrt(scale_factor)
    gadget.attrs["BoxSize"] = np.array([args.boxsize, args.boxsize, args.boxsize])

    log(comm, "Running FOF; note that FOF(...) calls run() internally")
    fof = FOF(gadget, linking_length=args.linking_length, nmin=args.nmin)

    log(comm, f"Converting FOF groups to HaloCatalog with mdef={args.mdef}")
    halos = fof.to_halos(
        particle_mass=particle_mass,
        cosmo=cosmology,
        redshift=args.redshift,
        mdef=args.mdef,
    )

    log(comm, "Computing local halo fields")
    local_size = len(halos["Mass"].compute())
    local_vrms, local_length = compute_local_fof_stats(fof, gadget, local_size)
    local_arr = compute_local_halo_catalog_array(halos, local_vrms, local_length, args.boxsize)
    local_count = np.array(len(local_arr), dtype="i8")
    offsets = np.cumsum([0] + comm.allgather(int(local_count)))[:-1]
    local_arr["halo_id"] = np.arange(offsets[comm.rank], offsets[comm.rank] + len(local_arr), dtype="i8")

    gathered = comm.gather(local_arr, root=0)
    if comm.rank == 0:
        outdir = os.path.dirname(args.output)
        if outdir:
            os.makedirs(outdir, exist_ok=True)
        arr = np.concatenate(gathered) if gathered else local_arr
        np.save(args.output, arr)

        metadata = {
            "snapshot": args.snapshot,
            "output": args.output,
            "boxsize": args.boxsize,
            "redshift": args.redshift,
            "linking_length": args.linking_length,
            "nmin": args.nmin,
            "mdef": args.mdef,
            "particle_type": args.particle_type,
            "particle_mass_code_units": particle_mass_code,
            "particle_mass_msun_h": particle_mass,
            "target": args.target,
            "model": args.model,
            "realization": args.realization,
            "cosmology": cosmology_parameters,
            "cosmology_source": "the simulation's run file",
            "schema": "flat halo catalog",
            "mass_field": "m200c",
            "mass_units": "Msun/h",
            "radius_units": "kpc/h",
            "position_units": "Mpc/h",
            "velocity_units": "km/s",
            "velocity_definition": VELOCITY_DEFINITION,
            "vrms_definition": "sqrt(sum_i(<v_i^2> - <v_i>^2)) over FOF member particles",
            "mpi_ranks": comm.size,
            "halos": int(len(arr)),
            "fields": list(arr.dtype.names),
        }
        with open(args.output + ".json", "w", encoding="utf-8") as handle:
            json.dump(metadata, handle, indent=2, sort_keys=True)
        log(comm, f"Saved {len(arr)} halos to {args.output}")
        log(comm, f"Saved metadata to {args.output}.json")

    comm.barrier()
    log(comm, "Done")
    return 0


if __name__ == "__main__":
    sys.exit(main())
