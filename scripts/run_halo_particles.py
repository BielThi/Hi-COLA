#!/usr/bin/env python3
"""FoF halo finding that also keeps a subsample of the halo member particles (tryout for HODDIES'
use_particles: satellites placed on member particles instead of an NFW profile).

Launched with MPI by jobs/halo_finder_particles.sh. Same FoF as scripts/run_halo_catalog.py (whose
functions it uses), run again because the particle labels are needed. It writes, in the run's
halo_catalog/ folder, three files that share the same halo_id (= FoF label - 1, labels sorted by size):

    halos200c_ptcl_catalog.npy    the halo catalog (same fields as halos200c_catalog.npy)
    halos200c_ptcl_hoddies.npy    its HODDIES array, with Rh and Rs in COMOVING kpc/h (r200c x (1+z));
                                  halos200c_hoddies.npy has the physical radii of nbodykit
    halos200c_ptcl_particles.npy  a random subsample (fraction --particle-fraction) of the member
                                  particles of every halo: x y z [Mpc/h], vx vy vz [peculiar km/s], halo_id

each with a .json of metadata. The files of the main halo finder are not touched.

    srun -n 512 python scripts/run_halo_particles.py --target Abacus_Small --model LCDM --realization s1234 \\
        --redshift 0.5 --boxsize 500 --linking-length 0.2 --nmin 50 --mdef 200c --h 0.6736 ... \\
        --particle-fraction 0.1 --particle-seed 1
"""

from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
from mpi4py import MPI

import hicola_paths as hp
from halo_pipeline import VELOCITY_DEFINITION, halo_catalog_to_hoddies
from run_halo_catalog import (compute_local_fof_stats, compute_local_halo_catalog_array, log, read_particle_mass,
                              read_scale_factor, run_cosmology)

PARTICLE_DTYPE = [("x", "f4"), ("y", "f4"), ("z", "f4"), ("vx", "f4"), ("vy", "f4"), ("vz", "f4"), ("halo_id", "i8")]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--target", required=True)
    parser.add_argument("--model", required=True)
    parser.add_argument("--realization", required=True)
    parser.add_argument("--redshift", type=float, required=True)
    parser.add_argument("--boxsize", type=float, required=True, help="Mpc/h")
    parser.add_argument("--linking-length", type=float, required=True, help="in units of the mean particle separation")
    parser.add_argument("--nmin", type=int, required=True, help="minimum number of particles per halo")
    parser.add_argument("--mdef", required=True, choices=["200c"])
    parser.add_argument("--particle-type", type=int, default=1, help="Gadget particle type of the dark matter")
    cosmo = parser.add_argument_group("cosmology (of the simulation)")
    cosmo.add_argument("--h", type=float, required=True)
    cosmo.add_argument("--omega-cdm", type=float, required=True)
    cosmo.add_argument("--omega-b", type=float, required=True)
    cosmo.add_argument("--omega-mnu", type=float, default=0.0)
    cosmo.add_argument("--omega-lambda", type=float, required=True)
    cosmo.add_argument("--tcmb", type=float, required=True)
    cosmo.add_argument("--neff", type=float, required=True)
    cosmo.add_argument("--n-s", type=float, required=True)
    particles = parser.add_argument_group("member particles")
    particles.add_argument("--particle-fraction", type=float, required=True,
                           help="fraction of the member particles kept (random, the same in every halo)")
    particles.add_argument("--particle-seed", type=int, required=True, help="seed of the subsampling")
    args = parser.parse_args()
    if not 0 < args.particle_fraction <= 1:
        raise SystemExit("--particle-fraction must be in (0, 1]")
    run = (args.target, args.model, args.realization, args.redshift)
    args.snapshot = hp.snapshot_glob(*run)
    args.catalog_output = hp.ptcl_halo_catalog_file(*run)
    args.hoddies_output = hp.ptcl_hoddies_file(*run)
    args.particles_output = hp.ptcl_particles_file(*run)
    return args


def main() -> int:
    args = parse_args()
    comm = MPI.COMM_WORLD
    for name in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
        os.environ.setdefault(name, "1")

    log(comm, f"FoF with member particles on {comm.size} MPI ranks: {args.snapshot}")
    particle_mass = read_particle_mass(args.snapshot, args.particle_type) * 1.0e10
    scale_factor = read_scale_factor(args.snapshot)
    if not np.isclose(scale_factor, 1.0 / (1.0 + args.redshift), rtol=1e-3):
        raise SystemExit(f"--redshift {args.redshift} does not match the snapshot (a = {scale_factor:.5f})")

    from nbodykit.algorithms.fof import FOF
    from nbodykit.cosmology.cosmology import Cosmology
    from nbodykit.lab import Gadget1Catalog

    astropy_cosmology, cosmology_parameters = run_cosmology(args)
    cosmology = Cosmology.from_astropy(astropy_cosmology, n_s=args.n_s)
    gadget = Gadget1Catalog(args.snapshot, comm=comm)
    gadget["Velocity"] = gadget["GadgetVelocity"] * np.sqrt(scale_factor)  # v_pec / sqrt(a) -> peculiar km/s
    gadget.attrs["BoxSize"] = np.array([args.boxsize] * 3)

    log(comm, f"Running FoF (linking length {args.linking_length}, at least {args.nmin} particles)")
    fof = FOF(gadget, linking_length=args.linking_length, nmin=args.nmin)
    halos = fof.to_halos(particle_mass=particle_mass, cosmo=cosmology, redshift=args.redshift, mdef=args.mdef)

    # Halo catalog, as scripts/run_halo_catalog.py
    local_size = len(halos["Mass"].compute())
    local_vrms, local_length = compute_local_fof_stats(fof, gadget, local_size)
    nbodykit_length = np.rint(np.asarray(halos["Mass"].compute()) / particle_mass).astype("i8")  # Mass = particle mass x Length
    if not np.array_equal(local_length, nbodykit_length):
        raise SystemExit("The member counts per FoF label differ from nbodykit's halo lengths: halo_id != label - 1")
    local_arr = compute_local_halo_catalog_array(halos, local_vrms, local_length, args.boxsize)
    offsets = np.cumsum([0] + comm.allgather(len(local_arr)))[:-1]
    local_arr["halo_id"] = np.arange(offsets[comm.rank], offsets[comm.rank] + len(local_arr), dtype="i8")

    # Member particles: a random subsample of the particles with a FoF label (label 0 = no halo)
    labels = np.asarray(fof.labels)
    rng = np.random.default_rng([args.particle_seed, comm.rank])
    keep = (labels > 0) & (rng.random(labels.size) < args.particle_fraction)
    position = np.asarray(gadget.compute(gadget["Position"]), dtype="f8")[keep] % args.boxsize
    velocity = np.asarray(gadget.compute(gadget["Velocity"]), dtype="f8")[keep]
    local_particles = np.empty(int(keep.sum()), dtype=PARTICLE_DTYPE)
    local_particles["x"], local_particles["y"], local_particles["z"] = position.T
    for axis in "xyz":  # float32 rounding can turn a position just below the box size into the box size
        local_particles[axis][local_particles[axis] >= args.boxsize] = 0.0
    local_particles["vx"], local_particles["vy"], local_particles["vz"] = velocity.T
    local_particles["halo_id"] = labels[keep] - 1
    n_members = comm.allreduce(int(np.sum(labels > 0)), op=MPI.SUM)
    log(comm, f"Member particles: {n_members:,} in halos, keeping a fraction {args.particle_fraction}")

    halo_parts = comm.gather(local_arr, root=0)
    particle_parts = comm.gather(local_particles, root=0)
    if comm.rank == 0:
        catalog = np.concatenate(halo_parts)
        particles = np.concatenate(particle_parts)
        del particle_parts
        if particles["halo_id"].max() >= len(catalog):
            raise SystemExit("A particle's halo_id is beyond the halo catalog")
        per_halo = np.bincount(particles["halo_id"], minlength=len(catalog))
        expected = catalog["length"] * args.particle_fraction
        log(comm, f"Particles kept: {len(particles):,} ({len(particles) / n_members:.4f} of the members); "
                  f"per halo / (fraction x length): median {np.median(per_halo / expected):.3f}")

        hoddies = halo_catalog_to_hoddies(catalog)
        hoddies["Rh"] *= 1 + args.redshift  # nbodykit's r200c is physical; HODDIES places satellites in comoving units
        hoddies["Rs"] *= 1 + args.redshift

        common = {
            "snapshot": args.snapshot, "target": args.target, "model": args.model, "realization": args.realization,
            "redshift": args.redshift, "boxsize": args.boxsize, "linking_length": args.linking_length,
            "nmin": args.nmin, "mdef": args.mdef, "particle_mass_msun_h": particle_mass,
            "cosmology": cosmology_parameters, "mpi_ranks": comm.size, "halos": int(len(catalog)),
            "halo_id": "FoF label - 1 (labels sorted by halo size); the same in the three ptcl files",
            "velocity_definition": VELOCITY_DEFINITION, "velocity_units": "km/s", "position_units": "Mpc/h",
        }
        outputs = (
            (args.catalog_output, catalog, {"schema": "flat halo catalog", "mass_field": "m200c", "mass_units": "Msun/h",
                                            "radius_units": "kpc/h (physical r200c of nbodykit)",
                                            "vrms_definition": "sqrt(sum_i(<v_i^2> - <v_i>^2)) over FOF member particles"}),
            (args.hoddies_output, hoddies, {"schema": "HODDIES halo array", "radius_units": "kpc/h comoving (r200c x (1+z))"}),
            (args.particles_output, particles, {"schema": "member particles", "particle_fraction": args.particle_fraction,
                                                "particle_seed": args.particle_seed, "particles": int(len(particles)),
                                                "member_particles": n_members,
                                                "subsampling": "random per particle, seed (particle_seed, MPI rank)"}),
        )
        for path, array, extra in outputs:
            path.parent.mkdir(parents=True, exist_ok=True)
            np.save(path, array)
            with open(str(path) + ".json", "w", encoding="utf-8") as handle:
                json.dump({**common, **extra, "output": str(path), "fields": list(array.dtype.names)}, handle,
                          indent=2, sort_keys=True)
            log(comm, f"Saved {path} ({len(array):,} rows)")
    comm.barrier()
    log(comm, "Done")
    return 0


if __name__ == "__main__":
    sys.exit(main())
