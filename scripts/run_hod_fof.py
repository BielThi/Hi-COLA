#!/usr/bin/env python3
"""Run the expensive nbodykit FOF/HaloCatalog step for HOD.ipynb.

This script is intended to be launched with MPI, for example:

    srun -n 72 python /home/biel/Hi-COLA/scripts/run_hod_fof.py

It writes a structured NumPy array with the same fields that the notebook's
HODDIES setup expects.
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


DEFAULT_SNAPSHOT = "/scratch/biel/LCDM_Bart/snapshot_LCDM_Bart_z0.000/gadget_z0.000.*"
DEFAULT_OUTPUT = "/scratch/biel/LCDM_Bart/hod_fof/halos200c_hod_arr.npy"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--snapshot", default=DEFAULT_SNAPSHOT, help="Gadget snapshot glob")
    parser.add_argument("--output", default=DEFAULT_OUTPUT, help="Output .npy file")
    parser.add_argument("--boxsize", type=float, default=1024.0)
    parser.add_argument("--redshift", type=float, default=0.0)
    parser.add_argument("--linking-length", type=float, default=0.2)
    parser.add_argument("--nmin", type=int, default=50)
    parser.add_argument("--mdef", default="200c")
    parser.add_argument("--particle-type", type=int, default=1)
    return parser.parse_args()


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


def compute_local_hod_array(halos) -> np.ndarray:
    mass = halos["Mass"].compute()
    position = halos["Position"].compute()
    velocity = halos["Velocity"].compute()
    radius = halos["Radius"].compute()
    concentration = halos["Concentration"].compute()

    dtype = [
        ("x", "f8"),
        ("y", "f8"),
        ("z", "f8"),
        ("vx", "f8"),
        ("vy", "f8"),
        ("vz", "f8"),
        ("Mh", "f8"),
        ("log10_Mh", "f8"),
        ("Rh", "f8"),
        ("Rs", "f8"),
        ("c", "f8"),
        ("Vrms", "f8"),
        ("halo_id", "i8"),
    ]

    arr = np.zeros(len(mass), dtype=dtype)
    arr["x"], arr["y"], arr["z"] = position.T
    arr["vx"], arr["vy"], arr["vz"] = velocity.T
    arr["Mh"] = mass
    arr["log10_Mh"] = np.log10(mass)
    arr["Rh"] = 1000.0 * radius
    arr["c"] = concentration
    arr["Rs"] = arr["Rh"] / concentration
    return arr


def compute_local_fof_vrms(fof, source, local_size: int) -> np.ndarray:
    """Compute HODDIES' 3D per-halo RMS velocity from FOF particle labels."""
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

    # nbodykit's halo catalog drops label 0 and scatters the remaining labels
    # evenly across ranks, so slice the global Vrms array with the local offset.
    local_start = sum(comm.allgather(local_size)[: comm.rank])
    return vrms[1:][local_start : local_start + local_size]


def main() -> int:
    args = parse_args()
    comm = MPI.COMM_WORLD

    os.environ.setdefault("OMP_NUM_THREADS", "1")
    os.environ.setdefault("MKL_NUM_THREADS", "1")
    os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

    log(comm, f"Starting HOD FOF on {comm.size} MPI ranks")
    log(comm, f"Snapshot: {args.snapshot}")

    particle_mass_code = read_particle_mass(args.snapshot, args.particle_type)
    particle_mass = particle_mass_code * 1.0e10
    log(comm, f"Particle mass from Gadget header: {particle_mass_code} x 1e10 Msun/h")

    from astropy.cosmology import Planck18
    from nbodykit.algorithms.fof import FOF
    from nbodykit.cosmology.cosmology import Cosmology
    from nbodykit.lab import Gadget1Catalog

    planck18 = Cosmology.from_astropy(Planck18)

    log(comm, "Opening Gadget catalog")
    gadget = Gadget1Catalog(args.snapshot, comm=comm)
    gadget["Velocity"] = gadget["GadgetVelocity"]
    gadget.attrs["BoxSize"] = np.array([args.boxsize, args.boxsize, args.boxsize])

    log(comm, "Running FOF; note that FOF(...) calls run() internally")
    fof = FOF(gadget, linking_length=args.linking_length, nmin=args.nmin)

    log(comm, f"Converting FOF groups to HaloCatalog with mdef={args.mdef}")
    halos = fof.to_halos(
        particle_mass=particle_mass,
        cosmo=planck18,
        redshift=args.redshift,
        mdef=args.mdef,
    )

    log(comm, "Computing local halo fields")
    local_arr = compute_local_hod_array(halos)
    local_arr["Vrms"] = compute_local_fof_vrms(fof, gadget, len(local_arr))
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
            "mass_units": "Msun/h",
            "radius_units": "kpc/h",
            "velocity_units": "km/s",
            "Vrms_definition": "sqrt(sum_i(<v_i^2> - <v_i>^2)) over FOF member particles",
            "mpi_ranks": comm.size,
            "halos": int(len(arr)),
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
