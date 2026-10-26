#!/usr/bin/env python3
"""Make HODDIES galaxy mocks from the halo array of a run and measure their redshift-space P(k)
and marked P(k).

Input:  halo_catalog/halos200c_hoddies.npy of the run (written by jobs/halo_finder.sh), with
        peculiar velocities (catalogs made before the velocity fix of 2026-10-14 are refused).
Output: halo_catalog/hod/<name>_PS.txt (+ .h5, .txt.json) in the run folder, one per seed
        (--name contains {seed}, e.g. LRG_seed{seed} -> LRG_seed42_PS.txt): the power spectrum
        multipoles of the mock in redshift space, measured as in the EMC/acm pipeline
        (galaxy_power_spectrum.py): columns k [h/Mpc], N_modes, P_ell [(Mpc/h)^3] shot-noise
        subtracted; the .h5 is jaxpower's own format, the same as the EMC measurements.
        halo_catalog/hod/<name>_MPS.txt (+ .h5, .txt.json): the marked power spectrum of the same
        mock, same columns and estimator settings (galaxy_power_spectrum.compute_marked_spectrum,
        as acm), mark m = sum_n c_n delta_R^n with --mark-rcut R and --mark-coefficients c_0 c_1 ...;
        the .h5 has the mark's m_bar, m_min, m_max as attributes, as in acm.

Satellites (--satellites): on random member particles of their halo (halos200c_particles.npy,
written by jobs/halo_finder.sh with PARTICLE_FRACTION > 0; halos with fewer stored particles than
satellites place the rest on the NFW profile), or on HODDIES' NFW profile, as in jobs/hod_fit.sh.

HODDIES: the clone next to this repository (version: PIPELINE.md, section 4), whose commit is printed and saved in the metadata.
The HOD model and its parameters come from the HODDIES parameter file. If it sets the tracer's
`density`, HODDIES rescales Ac (and the satellites with it) to reach that number density; if
that needs Ac > 1, it keeps Ac and warns (the density is then not reached).
The mock depends on the seed and on the number of threads (HODDIES seeds each block of halos).
Each seed gets its own halo catalogue and HOD object (the halo array is loaded once), so a seed
gives the same mock whether it runs alone or in a list.
Redshift-space positions: x_los + v_los / (a H(z)), with H(z) from the run's cosmology.

    python scripts/run_hod_mock.py --target Abacus_Small --model LCDM --realization s1234 \\
        --redshift 0.5 --boxsize 500 --h 0.6736 --omega-cdm 0.26589 --omega-b 0.049302 \\
        --tcmb 2.726 --neff 3.046 --hod-param-file Notebooks/HOD_parameter_file_z0.500.yaml \\
        --tracer LRG --seeds 42 43 44 --nthreads 32 --name "LRG_seed{seed}" --satellites particles \\
        --los z --meshsize 512 --ells 0 2 4 --k-step 0.001 --mark-rcut 10 --mark-coefficients 1 -0.5
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

import hicola_paths as hp
from galaxy_power_spectrum import compute_marked_spectrum, compute_spectrum, redshift_space_positions, save_spectrum
from halo_pipeline import load_catalog_metadata, require_peculiar_velocities


def import_hoddies() -> dict:
    """HOD and Base_catalogue of the HODDIES clone next to this repository (not the copy installed
    in the env), with its git commit."""
    sys.path.insert(0, str(hp.HODDIES_ROOT))
    import HODDIES
    from HODDIES import HOD
    from HODDIES.sim_loader import Base_catalogue

    git = ["git", "-C", str(hp.HODDIES_ROOT)]
    commit = subprocess.run(git + ["describe", "--always", "--dirty"], capture_output=True, text=True).stdout.strip()
    branch = subprocess.run(git + ["rev-parse", "--abbrev-ref", "HEAD"], capture_output=True, text=True).stdout.strip()
    version = f"{branch} {commit}".strip() or "unknown"
    print(f"HODDIES:        {HODDIES.__file__} ({version})")
    return {"HOD": HOD, "Base_catalogue": Base_catalogue, "version": version}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--target", required=True, help="e.g. COLA_NM4096, Abacus_Small")
    parser.add_argument("--model", required=True, help="e.g. LCDM, nDGP_N1.2")
    parser.add_argument("--realization", required=True, help="e.g. s1234, s1234-fixed")
    parser.add_argument("--redshift", type=float, required=True)
    parser.add_argument("--boxsize", type=float, required=True, help="Mpc/h")
    cosmo = parser.add_argument_group("cosmology (of the simulation, for H(z) in the redshift-space shift)")
    cosmo.add_argument("--h", type=float, required=True)
    cosmo.add_argument("--omega-cdm", type=float, required=True)
    cosmo.add_argument("--omega-b", type=float, required=True)
    cosmo.add_argument("--omega-mnu", type=float, default=0.0)
    cosmo.add_argument("--tcmb", type=float, required=True, help="K")
    cosmo.add_argument("--neff", type=float, required=True)
    hod = parser.add_argument_group("HOD mock")
    hod.add_argument("--hod-param-file", type=Path, required=True, help="HODDIES parameter file (.yaml)")
    hod.add_argument("--tracer", required=True, help="tracer of the parameter file, e.g. LRG")
    hod.add_argument("--seeds", type=int, nargs="+", required=True, help="one mock and P(k) per seed")
    hod.add_argument("--nthreads", type=int, required=True, help="HODDIES threads (the mock depends on it)")
    hod.add_argument("--satellites", required=True, choices=["particles", "nfw"],
                     help="satellites on member particles (halos200c_particles.npy) or on HODDIES' NFW profile")
    hod.add_argument("--name", required=True,
                     help="output name with {seed}, e.g. LRG_seed{seed} -> hod/LRG_seed42_PS.txt")
    pk = parser.add_argument_group("power spectrum")
    pk.add_argument("--los", required=True, choices=["x", "y", "z"], help="line of sight (box axis)")
    pk.add_argument("--meshsize", type=int, required=True)
    pk.add_argument("--ells", type=int, nargs="+", required=True, help="multipoles, e.g. 0 2 4")
    pk.add_argument("--k-step", type=float, required=True, help="k-bin width [h/Mpc]")
    mark = parser.add_argument_group("marked power spectrum (same settings as the P(k))")
    mark.add_argument("--mark-rcut", type=float, required=True,
                      help="Gaussian smoothing radius R of the density contrast in the mark [Mpc/h]")
    mark.add_argument("--mark-coefficients", type=float, nargs="+", required=True,
                      help="c_0 c_1 ... of the mark m = sum_n c_n delta_R^n")
    parser.add_argument("--halo-array", type=Path, default=None, help="Default: halos200c_hoddies.npy of the run")
    args = parser.parse_args()

    if "{seed}" not in args.name:
        raise SystemExit(f"--name must contain {{seed}} so that every seed has its own file, got {args.name!r}")
    if not args.hod_param_file.is_file():
        raise SystemExit(f"No HOD parameter file {args.hod_param_file}")
    if len(set(args.seeds)) != len(args.seeds):
        raise SystemExit(f"--seeds has duplicates: {args.seeds}")
    run = (args.target, args.model, args.realization, args.redshift)
    args.halo_array = args.halo_array or hp.hoddies_catalog_file(*run)
    return args


def efunc(args: argparse.Namespace) -> float:
    """H(z)/H0 of the simulation's cosmology (flat, with radiation from T_CMB and N_eff)."""
    from astropy import units as u
    from astropy.cosmology import FlatLambdaCDM

    cosmology = FlatLambdaCDM(
        H0=100.0 * args.h, Om0=args.omega_cdm + args.omega_b + args.omega_mnu, Ob0=args.omega_b,
        Tcmb0=args.tcmb * u.K, Neff=args.neff,
    )
    return float(cosmology.efunc(args.redshift))


def main() -> int:
    args = parse_args()
    hoddies = import_hoddies()

    halos = np.load(args.halo_array)
    catalog_file = hp.halo_catalog_file(args.target, args.model, args.realization, args.redshift)
    metadata = load_catalog_metadata(catalog_file)
    require_peculiar_velocities(metadata, catalog_file)
    if "boxsize" in metadata and not np.isclose(float(metadata["boxsize"]), args.boxsize):
        raise SystemExit(f"--boxsize {args.boxsize} differs from the halo catalog's {metadata['boxsize']}")
    print(f"Halos:          {args.halo_array} ({len(halos):,} halos)")
    particles = None
    if args.satellites == "particles":
        particle_file = hp.particles_file(args.target, args.model, args.realization, args.redshift)
        if not particle_file.exists():
            raise SystemExit(f"No member particles {particle_file}: run jobs/halo_finder.sh with PARTICLE_FRACTION > 0")
        particles = np.load(particle_file)
        particle_info = json.loads(Path(str(particle_file) + ".json").read_text())
        print(f"Satellites:     particles, {particle_file} ({len(particles):,} particles, "
              f"fraction {particle_info['particle_fraction']})")
    else:
        print("Satellites:     nfw")

    e_z = efunc(args)
    print(f"Seeds:          {args.seeds} ({len(args.seeds)} mocks)")
    for seed in args.seeds:
        print(f"---------------- seed {seed}", flush=True)
        measure_seed(hoddies, halos, particles, seed, e_z, catalog_file, args)
    return 0


def measure_seed(hoddies: dict, halos: np.ndarray, particles, seed: int, e_z: float, catalog_file: Path,
                 args: argparse.Namespace) -> None:
    """One HOD mock with this seed, its redshift-space P(k) and marked P(k) multipoles, written to
    hod/<name>_PS.txt and hod/<name>_MPS.txt."""
    name = args.name.format(seed=seed)
    run = (args.target, args.model, args.realization, args.redshift, name)
    output = hp.hod_spectrum_file(*run, "PS")
    marked_output = hp.hod_spectrum_file(*run, "MPS")
    base_catalog = hoddies["Base_catalogue"](
        {column: halos[column].copy() for column in halos.dtype.names},  # a dict: a structured array is refused
        particle_data=None if particles is None else {column: particles[column].copy() for column in particles.dtype.names},
        boxsize=args.boxsize, z_simu=args.redshift, nthreads=args.nthreads,
    )
    hod = hoddies["HOD"](base_catalog, param_file=str(args.hod_param_file), nthreads=args.nthreads, seed=seed,
                           use_particles=particles is not None)
    if (particles is not None) != (hod.args["use_particles"] and hod.part_subsamples is not None):
        raise SystemExit("HODDIES' use of the particles differs from --satellites")
    if hod.nthreads != args.nthreads:
        raise SystemExit(f"HODDIES uses {hod.nthreads} threads instead of --nthreads {args.nthreads} (it is capped by "
                         f"numba's thread count, NUMBA_NUM_THREADS); the mock depends on it")
    if args.tracer not in hod.args["tracers"]:
        raise SystemExit(f"Tracer {args.tracer} is not in {args.hod_param_file} (tracers: {hod.args['tracers']})")
    hod_parameters = hod.args[args.tracer]
    print(f"HOD parameters: {json.dumps(hod_parameters, default=str)}")

    mock = hod.make_mock_cat(tracers=args.tracer, fix_seed=seed)
    n_gal = len(mock["x"])
    n_cent = int(np.sum(np.asarray(mock["Central"]) == 1))
    density = n_gal / args.boxsize**3
    print(f"Galaxies:       {n_gal:,} ({n_cent:,} centrals, satellite fraction {1 - n_cent / n_gal:.3f}), "
          f"n = {density:.3e} (h/Mpc)^3")

    positions = redshift_space_positions(
        mock["x"], mock["y"], mock["z"], mock["vx"], mock["vy"], mock["vz"],
        boxsize=args.boxsize, redshift=args.redshift, efunc=e_z, los=args.los,
    )
    print(f"Redshift space: along {args.los}, H(z)/H0 = {e_z:.6f}")

    metadata = {
        "run": hp.simulation_name(args.target, args.model, args.realization),
        "redshift": args.redshift,
        "boxsize": args.boxsize,
        "halo_catalog": str(catalog_file),
        "hod_param_file": str(args.hod_param_file),
        "tracer": args.tracer,
        "satellites": args.satellites,
        "hod_parameters": hod_parameters,
        "seed": seed,
        "nthreads": args.nthreads,
        "hoddies": hoddies["version"],
        "n_galaxies": n_gal,
        "n_centrals": n_cent,
        "number_density_h3_Mpc3": density,
        "redshift_space": {"los": args.los, "efunc": e_z, "shift": "v_los / (a H(z))"},
    }
    estimator = {
        "code": "jaxpower, as acm scripts/emc/measurements/abacus_small.py compute_spectrum",
        "meshsize": args.meshsize, "resampler": "tsc", "interlacing": 3, "compensate": True,
        "ells": args.ells, "k_step": args.k_step, "shot_noise": "subtracted",
    }
    output.parent.mkdir(parents=True, exist_ok=True)

    start = time.time()
    spectrum = compute_spectrum(positions, args.boxsize, args.meshsize, args.ells, args.los, args.k_step)
    print(f"Power spectrum: done in {time.time() - start:.1f} s (mesh {args.meshsize}^3, ells {args.ells}, k-step {args.k_step})")
    save_spectrum(spectrum, args.ells, output, {**metadata, "statistic": "power spectrum", "estimator": estimator})
    print(f"Saved:          {output} (+ .h5, .json)")

    start = time.time()
    marked_spectrum, mark_summary = compute_marked_spectrum(
        positions, args.boxsize, args.meshsize, args.ells, args.los, args.k_step,
        smoothing_radius=args.mark_rcut, mark_coefficients=args.mark_coefficients,
    )
    print(f"Marked P(k):    done in {time.time() - start:.1f} s (R = {args.mark_rcut} Mpc/h, "
          f"coefficients {args.mark_coefficients}, m_bar {mark_summary['m_bar']:.4f}, "
          f"m in [{mark_summary['m_min']:.4f}, {mark_summary['m_max']:.4f}])")
    mark = {
        "definition": "m = sum_n c_n delta_R^n, delta_R: galaxy density contrast smoothed with a Gaussian of radius R; "
                      "normalised by m_bar = <m n_g> / <n_g>",
        "rcut": args.mark_rcut, "coefficients": args.mark_coefficients, **mark_summary,
    }
    marked_estimator = {
        **estimator,
        "code": "jaxpower, as acm compute_marked_spectrum_acm (MarkedPowerSpectrumMultipoles, box without randoms)",
        "normalization": "box (unmarked galaxies)",
    }
    save_spectrum(
        marked_spectrum, args.ells, marked_output,
        {**metadata, "statistic": "marked power spectrum", "mark": mark, "estimator": marked_estimator},
        h5_attrs=mark_summary,
    )
    print(f"Saved:          {marked_output} (+ .h5, .json)")

if __name__ == "__main__":
    raise SystemExit(main())
