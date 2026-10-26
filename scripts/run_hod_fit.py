#!/usr/bin/env python3
"""Fit the HOD parameters on the halo catalog of a run to measured clustering, with HODDIES.

Run by jobs/hod_fit.sh (MPI: one HOD object per rank; the ranks share the evaluations of each
iteration of the minimizer).

Data:   wp(rp) and/or the multipoles xi_ell(s) of pycorr two-point estimators (e.g. DESI Y1 LRG),
        rebinned and restricted to the chosen separation ranges (bins entirely inside them), which
        are also the bins of the correlation matrix. The fit can be restricted further to
        --fit-r-min <= rp, s <= --fit-r-max (bins entirely inside; 'None' = no limit on that side).
        Covariance: the jackknife errors of the data times a correlation matrix (e.g. from mocks)
        covering the statistics --correlation-statistics (its blocks of the fitted ones are
        used), or diagonal ("none").
Model:  HODDIES mocks on the run's halos (jobs/halo_finder.sh: halos200c_hoddies.npy, Rh
        comoving), all with the same seed, measured with the bins of the data (pycorr, redshift
        space along the line of sight with H(z) of the run's cosmology); chi2 as HODDIES'
        func_stochopy. Satellites (--satellites):
          particles  a random member particle of the halo (halos200c_particles.npy, written by
                     jobs/halo_finder.sh with PARTICLE_FRACTION > 0), velocity
                     v_halo + f_sigv (v_particle - v_halo); halos with fewer stored particles than
                     satellites place the rest on the NFW profile (fraction recorded)
          nfw        HODDIES' NFW profile (Rh, concentration c), velocity dispersion from Vrms x f_sigv
Fit:    a stochopy minimizer (e.g. cmaes) over the free HOD parameters within their priors; the
        other HOD parameters are those of the HOD parameter file.

Output (hp.hod_fit_dir, kept in hicola_data next to the run's measurements), in a folder named
after what is fitted: <tracer>_<HOD model>_<n free>p_<satellites>[_<statistics>] (HODDIES' fit model
name, ptcl or nfw, the statistics only when not all of wp + xi_0 + xi_2, and _r<min>-<max> when
the fit range is restricted), e.g. LRG_SHOD_6p_ptcl, LRG_SHOD_6p_nfw_wp, LRG_SHOD_6p_ptcl_r0.1-10,
LRG_SHOD_6p_ptcl_rmin-10; a new fit of the same thing replaces it:
    settings.json              the settings of the fit, the run and the code versions
    <HOD parameter file>.yaml  copy of the HOD parameter file
    data.npz                   the fitted data: separations, values, errors, data vector, covariance
    evaluations/rank_<k>.txt   every evaluation made by rank k: parameters, chi2, galaxies (written as it goes)
    iterations.txt             per iteration of the minimizer: its best chi2 and parameters (as it goes)
    best_fit.json              best-fit parameters (best evaluation of the whole fit), chi2, minimizer status
    best_fit.yaml              HOD parameter file of the best fit, for the mocks of jobs/hod_mock.sh
    best_fit_model.npz         the model at the best fit (fit seed and extra seeds) and its HOD occupation
    reference_model.npz        the same for the HOD of a reference parameter file (if given), for comparison

    srun --mpi=pmix -n 18 python scripts/run_hod_fit.py --target Abacus_Small --model LCDM \\
        --realization s1234 --redshift 0.5 --boxsize 500 ... (see jobs/hod_fit.sh)
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import time
import traceback
from pathlib import Path

import numpy as np

import hicola_paths as hp
from halo_pipeline import load_catalog_metadata, require_peculiar_velocities

STATISTICS = ("wp", "xi_ells")  # statistics this script can fit (and their order in the data vector)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    run = parser.add_argument_group("run")
    run.add_argument("--target", required=True)
    run.add_argument("--model", required=True)
    run.add_argument("--realization", required=True)
    run.add_argument("--redshift", type=float, required=True)
    run.add_argument("--boxsize", type=float, required=True, help="Mpc/h")
    cosmo = parser.add_argument_group("cosmology of the run (H(z) of the redshift-space distortions)")
    cosmo.add_argument("--h", type=float, required=True)
    cosmo.add_argument("--omega-cdm", type=float, required=True)
    cosmo.add_argument("--omega-b", type=float, required=True)
    cosmo.add_argument("--omega-mnu", type=float, required=True)
    cosmo.add_argument("--tcmb", type=float, required=True, help="K")
    cosmo.add_argument("--neff", type=float, required=True)
    cosmo.add_argument("--n-s", type=float, required=True)
    hod = parser.add_argument_group("HOD")
    hod.add_argument("--hod-param-file", type=Path, required=True, help="HODDIES parameter file (HOD model and fixed values)")
    hod.add_argument("--tracer", required=True)
    hod.add_argument("--priors", nargs="+", required=True, help="free parameters as name:min:max")
    hod.add_argument("--nthreads", type=int, required=True, help="HODDIES and pycorr threads per rank (the mock depends on it)")
    data = parser.add_argument_group("data")
    data.add_argument("--statistics", nargs="+", required=True, choices=STATISTICS)
    data.add_argument("--correlation", required=True, help="correlation matrix file (text), or 'none' (diagonal)")
    data.add_argument("--correlation-statistics", nargs="+", required=True, choices=STATISTICS,
                      help="statistics covered by the correlation matrix, in this order (wp, then xi_ell for each ell), "
                           "with the same bins as the data; the blocks of the fitted ones are used")
    wp = parser.add_argument_group("wp (required when fitting wp or when the correlation matrix covers it)")
    wp.add_argument("--data-rppi", type=Path, help="pycorr estimator in (rp, pi)")
    wp.add_argument("--pimax", type=float)
    wp.add_argument("--rp-rebin", type=int)
    wp.add_argument("--rp-range", type=float, nargs=2, metavar=("MIN", "MAX"))
    xi = parser.add_argument_group("xi_ell (required when fitting xi_ells or when the correlation matrix covers it)")
    xi.add_argument("--data-smu", type=Path, help="pycorr estimator in (s, mu)")
    xi.add_argument("--ells", type=int, nargs="+")
    xi.add_argument("--s-rebin", type=int)
    xi.add_argument("--s-range", type=float, nargs=2, metavar=("MIN", "MAX"))
    data.add_argument("--fit-r-min", required=True, help="smallest rp and s fitted [Mpc/h], or 'None'")
    data.add_argument("--fit-r-max", required=True, help="largest rp and s fitted [Mpc/h], or 'None'")
    data.add_argument("--los", required=True, choices=["x", "y", "z"], help="line of sight of the model (box axis)")
    data.add_argument("--satellites", required=True, choices=["particles", "nfw"],
                      help="satellites on member particles or on HODDIES' NFW profile")
    fit = parser.add_argument_group("minimizer")
    fit.add_argument("--method", required=True, help="stochopy method, e.g. cmaes, de, pso, cpso, vdcma, na")
    fit.add_argument("--maxiter", type=int, required=True)
    fit.add_argument("--popsize", type=int, required=True)
    fit.add_argument("--xtol", type=float, required=True)
    fit.add_argument("--start", required=True, choices=["random", "param_file"],
                     help="starting point: drawn in the priors (with --start-seed) or the values of the HOD parameter file")
    fit.add_argument("--start-seed", type=int, required=True, help="seed of the starting point and of the minimizer")
    fit.add_argument("--mock-seed", type=int, required=True, help="HOD seed of every mock of the fit")
    fit.add_argument("--model-seeds", type=int, nargs="+", required=True, help="extra seeds for the best-fit model")
    fit.add_argument("--reference-param-file", required=True,
                     help="HODDIES parameter file whose HOD is modelled at the end for comparison, or 'none'")
    parser.add_argument("--job-log", default="", help="SLURM log of the job (recorded in settings.json)")
    args = parser.parse_args()

    args.priors = parse_priors(args.priors)
    if len(set(args.statistics)) != len(args.statistics):
        raise SystemExit(f"--statistics has duplicates: {args.statistics}")
    args.statistics = [s for s in STATISTICS if s in args.statistics]  # the order of the data vector
    for name in ("fit_r_min", "fit_r_max"):
        value = getattr(args, name)
        try:
            setattr(args, name, None if value == "None" else float(value))
        except ValueError:
            raise SystemExit(f"--{name.replace('_', '-')} must be a number or 'None', got {value!r}")
    if None not in (args.fit_r_min, args.fit_r_max) and not args.fit_r_min < args.fit_r_max:
        raise SystemExit("--fit-r-min must be below --fit-r-max")
    args.correlation_statistics = [s for s in STATISTICS if s in args.correlation_statistics]
    if args.correlation != "none" and not set(args.statistics) <= set(args.correlation_statistics):
        raise SystemExit(f"The correlation matrix covers {args.correlation_statistics}, not all the fitted {args.statistics}")
    args.loaded_statistics = [s for s in STATISTICS if s in args.statistics
                              or (args.correlation != "none" and s in args.correlation_statistics)]
    needed = {"wp": ["data_rppi", "pimax", "rp_rebin", "rp_range"], "xi_ells": ["data_smu", "ells", "s_rebin", "s_range"]}
    for statistic in args.loaded_statistics:
        missing = ["--" + name.replace("_", "-") for name in needed[statistic] if getattr(args, name) is None]
        if missing:
            raise SystemExit(f"{statistic} (fitted or covered by the correlation matrix) needs {' '.join(missing)}")
    if args.omega_mnu != 0:
        raise SystemExit("Massive neutrinos are not handled in the H(z) of the fit (omega_mnu must be 0)")
    return args


def parse_priors(items: list[str]) -> dict[str, list[float]]:
    priors = {}
    for item in items:
        try:
            name, low, high = item.split(":")
            low, high = float(low), float(high)
        except ValueError:
            raise SystemExit(f"Prior {item!r} is not name:min:max")
        if not low < high:
            raise SystemExit(f"Prior {item!r}: min must be below max")
        if name in priors:
            raise SystemExit(f"Parameter {name} has two priors")
        priors[name] = [low, high]
    return priors


# ---------------------------------------------------------------- data


def select_bins(edges: np.ndarray, separation_range) -> slice:
    """Bins of a pycorr estimator (first axis) entirely inside [min, max]."""
    low, high = separation_range
    inside = np.flatnonzero((edges[:-1] >= low * (1 - 1e-9)) & (edges[1:] <= high * (1 + 1e-9)))
    if inside.size == 0:
        raise SystemExit(f"No bin between {low} and {high} (edges {edges[0]:.4g} to {edges[-1]:.4g})")
    return slice(inside[0], inside[-1] + 1)


def fit_bins(edges: np.ndarray, args: argparse.Namespace) -> slice:
    """Bins (of the data selected by RP_RANGE or S_RANGE) entirely inside the fit range."""
    low = -np.inf if args.fit_r_min is None else args.fit_r_min
    high = np.inf if args.fit_r_max is None else args.fit_r_max
    return select_bins(edges, (low, high))


def load_data(args: argparse.Namespace) -> dict:
    """Data vector (fitted statistics in the order of args.statistics, bins in the fit range), the
    measured statistics with their errors and bin edges (fit range), and the covariance (block of
    the fitted statistics and bins)."""
    from pycorr import TwoPointCorrelationFunction

    data = {}
    values, errors = {}, {}
    full_size, kept = {}, {}   # per statistic: size in the correlation matrix's layout, indices of the fitted bins
    if "wp" in args.loaded_statistics:
        estimator = TwoPointCorrelationFunction.load(str(args.data_rppi))[::args.rp_rebin]
        estimator = estimator[select_bins(estimator.edges[0], args.rp_range)]
        n_rp = estimator.shape[0]
        in_range = fit_bins(estimator.edges[0], args)
        estimator = estimator[in_range]
        full_size["wp"], kept["wp"] = n_rp, np.arange(n_rp)[in_range]
        if args.pimax > estimator.edges[1][-1]:
            raise SystemExit(f"--pimax {args.pimax} beyond the pi edges of the data ({estimator.edges[1][-1]})")
        rp, wp, wp_err = estimator(pimax=args.pimax, return_sep=True, return_std=True)
        data.update(rp=rp, wp=wp, wp_err=wp_err, edges_rp=estimator.edges[0], edges_pi=estimator.edges[1])
        values["wp"], errors["wp"] = wp, wp_err
    if "xi_ells" in args.loaded_statistics:
        estimator = TwoPointCorrelationFunction.load(str(args.data_smu))[::args.s_rebin]
        estimator = estimator[select_bins(estimator.edges[0], args.s_range)]
        n_s = estimator.shape[0]
        in_range = fit_bins(estimator.edges[0], args)
        estimator = estimator[in_range]
        full_size["xi_ells"] = n_s * len(args.ells)
        kept["xi_ells"] = np.concatenate([i * n_s + np.arange(n_s)[in_range] for i in range(len(args.ells))])
        s, xi, cov = estimator.get_corr(ells=tuple(args.ells), return_sep=True, return_cov=True)
        xi_err = np.sqrt(np.diag(cov)).reshape(len(args.ells), -1)
        data.update(s=s, xi=np.asarray(xi), xi_err=xi_err, ells=np.array(args.ells),
                    edges_s=estimator.edges[0], edges_mu=estimator.edges[1])
        values["xi_ells"], errors["xi_ells"] = np.ravel(xi), np.ravel(xi_err)
    data_vec = np.concatenate([values[s] for s in args.statistics])
    data_err = np.concatenate([errors[s] for s in args.statistics])

    if args.correlation == "none":
        correlation = np.eye(data_vec.size)
    else:
        full = np.loadtxt(args.correlation)
        sizes = {s: full_size[s] for s in args.correlation_statistics}
        if full.shape != (sum(sizes.values()),) * 2:
            raise SystemExit(f"The correlation matrix {args.correlation} is {full.shape}, but "
                             f"{' + '.join(f'{s} ({n})' for s, n in sizes.items())} have {sum(sizes.values())} values: "
                             "the ranges, rebinning or ells differ from those of the matrix")
        start = dict(zip(sizes, np.cumsum([0] + list(sizes.values()))[:-1]))
        rows = np.concatenate([start[s] + kept[s] for s in args.statistics])
        correlation = full[np.ix_(rows, rows)]
    data.update(data_vec=data_vec, data_err=data_err, correlation=correlation,
                cov=np.outer(data_err, data_err) * correlation)
    return data


# ---------------------------------------------------------------- HODDIES


def import_hoddies() -> dict:
    """HODDIES of the clone next to this repository, with its version."""
    sys.path.insert(0, str(hp.HODDIES_ROOT))
    import HODDIES
    from HODDIES import HOD
    from HODDIES.fit_functions import compute_chi2, flatten_stats, values_in_boundaries
    from HODDIES.sim_loader import Base_catalogue

    return {
        "HOD": HOD, "Base_catalogue": Base_catalogue, "compute_chi2": compute_chi2,
        "flatten_stats": flatten_stats, "values_in_boundaries": values_in_boundaries,
        "file": HODDIES.__file__, "version": git_version(hp.HODDIES_ROOT),
    }


def git_version(folder: Path) -> str:
    git = ["git", "-C", str(folder)]
    commit = subprocess.run(git + ["describe", "--always", "--dirty"], capture_output=True, text=True).stdout.strip()
    branch = subprocess.run(git + ["rev-parse", "--abbrev-ref", "HEAD"], capture_output=True, text=True).stdout.strip()
    return f"{branch} {commit}".strip() or "unknown"


def run_cosmology(args: argparse.Namespace):
    """The run's cosmology (cosmoprimo), for H(z) in HODDIES' redshift-space distortions."""
    from cosmoprimo.fiducial import Cosmology

    return Cosmology(h=args.h, Omega_cdm=args.omega_cdm, Omega_b=args.omega_b, T_cmb=args.tcmb,
                     N_eff=args.neff, n_s=args.n_s, m_ncdm=[], engine="class")


def build_hod(hoddies: dict, halos: np.ndarray, particles, data: dict, args: argparse.Namespace, verbose: bool):
    """HOD object of the run with the fit settings: satellites (on the member particles if given),
    priors, statistics, bins of the data, data and covariance."""
    base_catalog = hoddies["Base_catalogue"](
        {column: halos[column] for column in halos.dtype.names},  # a dict: a structured array is refused
        particle_data=None if particles is None else {column: particles[column] for column in particles.dtype.names},
        boxsize=args.boxsize, z_simu=args.redshift, nthreads=args.nthreads, setup_logger=verbose,
    )
    hod = hoddies["HOD"](base_catalog, param_file=str(args.hod_param_file), nthreads=args.nthreads,
                         seed=args.mock_seed, setup_logger=verbose, use_particles=particles is not None)
    if (particles is not None) != (hod.args["use_particles"] and hod.part_subsamples is not None):
        raise SystemExit("HODDIES' use of the particles differs from --satellites")
    if hod.nthreads != args.nthreads:
        raise SystemExit(f"HODDIES uses {hod.nthreads} threads instead of {args.nthreads} (capped by numba's "
                         "thread count, NUMBA_NUM_THREADS); the mocks depend on it")
    if hod._tracers() != [args.tracer]:
        raise SystemExit(f"The HOD parameter file defines the tracers {hod._tracers()}; this fit needs only {args.tracer}")
    unknown = [name for name in args.priors if name not in hod.args[args.tracer]]
    if unknown:
        raise SystemExit(f"Free parameters {unknown} are not parameters of {args.tracer} in {args.hod_param_file}")
    hod.cosmo = run_cosmology(args)

    settings = hod.args["clustering_settings"]
    settings["rsd"] = True
    settings["los"] = args.los
    if "wp" in args.statistics:
        settings["xi_rppi"]["edges_rppi"] = (data["edges_rp"], data["edges_pi"])
        settings["xi_rppi"]["pimax"] = args.pimax
    if "xi_ells" in args.statistics:
        settings["xi_smu"]["edges_smu"] = (data["edges_s"], data["edges_mu"])
        settings["xi_smu"]["multipole_index"] = list(args.ells)

    fit_param = hod.args["fit_param"]
    fit_param["priors"] = {args.tracer: {name: list(bounds) for name, bounds in args.priors.items()}}  # replaces the file's
    fit_param["fit_statistics"] = list(args.statistics)
    fit_param["add_poisson_noise"] = False
    fit_param["hartlap_factor"] = 0
    hod.initialize_fit(data["data_vec"], data["cov"])
    return hod


FULL_STATISTICS = (["wp", "xi_ells"], [0, 2])  # fits of all of wp + xi_0 + xi_2 have no statistics in their name


def fit_name(hod, args: argparse.Namespace) -> str:
    """Name of the fit from what is fitted: HODDIES' fit model name (tracer, HOD model, extensions,
    number of free parameters), the satellites (ptcl or nfw) and, unless all of wp + xi_0 + xi_2
    are fitted, the statistics: e.g. LRG_SHOD_6p_ptcl, LRG_SHOD_6p_nfw_wp."""
    name = f"{hod.fit_model_name}_{'ptcl' if args.satellites == 'particles' else 'nfw'}"
    if (args.statistics, args.ells if "xi_ells" in args.statistics else None) != FULL_STATISTICS:
        labels = ["wp" if s == "wp" else "xi" + "".join(str(ell) for ell in args.ells) for s in args.statistics]
        name += "_" + "+".join(labels)
    if args.fit_r_min is not None or args.fit_r_max is not None:
        bound = lambda value, word: word if value is None else f"{value:g}"
        name += f"_r{bound(args.fit_r_min, 'min')}-{bound(args.fit_r_max, 'max')}"
    return name


class Objective:
    """chi2 of one point of the parameter space, as HODDIES' func_stochopy, recorded in a file."""

    def __init__(self, hod, hoddies: dict, mock_seed: int, record: Path, start_time: float):
        self.hod, self.hoddies, self.mock_seed, self.start_time = hod, hoddies, mock_seed, start_time
        self.record = open(record, "w", buffering=1)
        self.record.write("# index time_s " + " ".join(hod.name_params) + " chi2 n_gal n_cen failed eval_s\n")
        self.count, self.failed = 0, 0
        self.best_chi2, self.best_x = np.inf, None   # best evaluation of this rank

    def __call__(self, x: np.ndarray) -> float:
        start = time.time()
        hod = self.hod
        chi2, n_gal, n_cen, failed = np.inf, 0, 0, 0
        params = np.array([x], dtype="float64")
        if self.hoddies["values_in_boundaries"](params, hod.priors):
            params.dtype = [(name, "float64") for name in hod.name_params]
            try:
                hod.update_new_param(params, hod.name_params)
                cat = hod.make_mock_cat(fix_seed=self.mock_seed, verbose=False)
                stats = hod.compute_stats(cat, stat=hod._stats_to_fit, tracers=hod._tracer_to_fit, verbose=False)
                model = self.hoddies["flatten_stats"](stats, hod._stats_to_fit, hod._comb_tr_list)
                chi2 = float(self.hoddies["compute_chi2"](model, hod.data, hod.inv_cov2))
                n_gal, n_cen = int(cat.size), int(np.sum(cat["Central"] == 1))
            except Exception:
                failed = 1
                self.failed += 1
                print(f"Evaluation failed at {dict(zip(hod.name_params, np.ravel(x)))}:\n{traceback.format_exc()}",
                      file=sys.stderr, flush=True)
        self.record.write(f"{self.count} {start - self.start_time:.1f} " + " ".join(f"{v:.8g}" for v in np.ravel(x))
                          + f" {chi2:.8g} {n_gal} {n_cen} {failed} {time.time() - start:.2f}\n")
        self.count += 1
        if chi2 < self.best_chi2:
            self.best_chi2, self.best_x = chi2, np.array(x, dtype=float)
        return chi2


def hod_occupation(hod, tracer: str, log_mass: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Mean central and satellite occupation at the current parameters, with the density rescaling of
    Ac as in make_mock_cat."""
    p_cen, p_sat, _ = hod._HOD__init_hod_param(tracer)
    link = hod.args[tracer]["link_sat_to_central"]
    rescale = hod.get_ds_fac(tracer)
    if p_cen[0] * rescale <= 1:
        p_cen[0] *= rescale
        if p_sat is not None and not link:
            p_sat[0] *= rescale
    n_cen = np.array([hod._fun_cHOD[tracer](m, p_cen) for m in log_mass])
    if p_sat is None:
        return n_cen, np.zeros_like(n_cen)
    n_sat = np.nan_to_num([hod._fun_sHOD[tracer](m, p_sat) for m in log_mass])
    return n_cen, n_sat * (n_cen if link else 1)


def best_fit_model(hod, hoddies: dict, x: np.ndarray, args: argparse.Namespace) -> dict:
    """Model statistics at the best fit for the fit seed and the extra seeds, and the HOD occupation."""
    params = np.array([x], dtype="float64")
    params.dtype = [(name, "float64") for name in hod.name_params]
    hod.update_new_param(params, hod.name_params)
    return model_statistics(hod, hoddies, args)


def reference_hod(path: Path, tracer: str) -> dict:
    """HOD of the tracer in a reference parameter file (all its values, not only the free ones)."""
    import yaml

    content = yaml.safe_load(Path(path).read_text())
    if tracer not in (content or {}):
        raise SystemExit(f"The reference parameter file {path} has no {tracer} section")
    return content[tracer]


def write_best_fit_yaml(fit_dir: Path) -> Path:
    """best_fit.yaml in the fit folder, for jobs/hod_mock.sh: the fit's HOD parameter file with the
    tracer's HOD as used in the fit (settings.json) and its free parameters at the best fit
    (best_fit.json), and use_particles as the fit's satellites. Only reads files of the fit folder."""
    import yaml

    settings = json.loads((fit_dir / "settings.json").read_text())
    best = json.loads((fit_dir / "best_fit.json").read_text())
    tracer = settings["tracer"]
    content = yaml.safe_load((fit_dir / settings["hod_param_file"]).read_text())
    content[tracer] = dict(settings["hod_parameters"])
    for name in settings["priors"]:
        content[tracer][name] = best["parameters"][f"{name}_{tracer}"]
    content["use_particles"] = settings["satellites"] == "particles"
    path = fit_dir / "best_fit.yaml"
    header = (f"# Best-fit HOD of the fit {settings['fit_name']} on {settings['run']} at z = {settings['redshift']} "
              f"(jobs/hod_fit.sh, {settings['started']}):\n"
              f"# {settings['hod_param_file']} with {tracer} as in the fit and its free parameters "
              f"{', '.join(settings['priors'])} at the best fit\n"
              f"# (chi2 {best['chi2']:.2f} for {best['dof']} dof). Satellites of the fit: {settings['satellites']} "
              f"(SATELLITES in jobs/hod_mock.sh).\n")
    path.write_text(header + yaml.safe_dump(content, sort_keys=False))
    return path


def same_hod_value(key: str, a, b) -> bool:
    """Whether two values of a HOD parameter have the same effect (no smearing: None or 0; no
    assembly bias: all coefficients 0, whatever the variables listed)."""
    if key == "vsmear":
        return (a in (None, 0)) and (b in (None, 0)) or a == b
    if key == "assembly_bias" and isinstance(a, dict) and isinstance(b, dict):
        zero = lambda d: all(np.all(np.asarray(v) == 0) for v in d.values())
        return (zero(a) and zero(b)) or a == b
    return a == b


def model_statistics(hod, hoddies: dict, args: argparse.Namespace) -> dict:
    """Model statistics at the current HOD for the fit seed and the extra seeds, chi2, and the HOD occupation."""
    model = {"seeds": np.array([args.mock_seed, *args.model_seeds])}
    vectors, n_gal, n_cen, nfw = [], [], [], []
    particles = hod.args["use_particles"] and hod.part_subsamples is not None
    if particles:
        stored = np.bincount(np.asarray(hod.part_subsamples["halo_id"]), minlength=len(hod.hcat["halo_id"]))
    for seed in model["seeds"]:
        cat = hod.make_mock_cat(fix_seed=int(seed), verbose=False)
        hosts = np.asarray(cat["halo_id"])[np.asarray(cat["Central"]) == 0]
        if particles and hosts.size:  # satellites beyond the stored particles of their halo go on the NFW profile
            nfw.append(np.maximum(np.bincount(hosts, minlength=stored.size) - stored, 0).sum() / hosts.size)
        stats = hod.compute_stats(cat, stat=hod._stats_to_fit, tracers=hod._tracer_to_fit, verbose=False)
        vectors.append(hoddies["flatten_stats"](stats, hod._stats_to_fit, hod._comb_tr_list))
        n_gal.append(cat.size)
        n_cen.append(int(np.sum(cat["Central"] == 1)))
        if seed == args.mock_seed:
            pair = hod._comb_tr_list[0]
            if "wp" in stats:
                model["rp"] = np.asarray(stats["wp"][pair][0])
            if "xi_ells" in stats:
                model["s"] = np.asarray(stats["xi_ells"][pair][0])
    model.update(model_vec=np.array(vectors), n_gal=np.array(n_gal), n_cen=np.array(n_cen),
                 sat_nfw_fraction=np.array(nfw))
    model["chi2"] = np.array([hoddies["compute_chi2"](v, hod.data, hod.inv_cov2) for v in vectors])
    model["log_mass"] = np.linspace(11.0, 15.5, 200)
    model["n_cen_hod"], model["n_sat_hod"] = hod_occupation(hod, args.tracer, model["log_mass"])
    return model


# ---------------------------------------------------------------- main


def main() -> int:
    args = parse_args()
    from mpi4py import MPI
    from stochopy.optimize import minimize

    comm = MPI.COMM_WORLD
    rank, size = comm.Get_rank(), comm.Get_size()
    root = rank == 0
    start_time = time.time()

    def say(message=""):
        if root:
            print(message, flush=True)

    run = (args.target, args.model, args.realization, args.redshift)
    hoddies = import_hoddies()
    # Halos (and member particles) of jobs/halo_finder.sh
    halo_array = hp.hoddies_catalog_file(*run)
    catalog_file = hp.halo_catalog_file(*run)
    particle_file = hp.particles_file(*run) if args.satellites == "particles" else None

    # Checks, data and output folder on rank 0, then shared
    problem, data = None, None
    if root:
        try:
            missing = [f for f in (halo_array, particle_file) if f is not None and not f.exists()]
            if missing:
                raise SystemExit(f"Missing {', '.join(map(str, missing))}: run jobs/halo_finder.sh on this run and redshift"
                                 + (" with PARTICLE_FRACTION > 0" if particle_file is not None else ""))
            metadata = load_catalog_metadata(catalog_file)
            require_peculiar_velocities(metadata, catalog_file)
            if "boxsize" in metadata and not np.isclose(float(metadata["boxsize"]), args.boxsize):
                raise SystemExit(f"--boxsize {args.boxsize} differs from the halo catalog's {metadata['boxsize']}")
            data = load_data(args)
            if args.reference_param_file != "none":
                reference_hod(args.reference_param_file, args.tracer)
        except SystemExit as exc:
            problem = str(exc)
    problem = comm.bcast(problem, root=0)
    if problem:
        say(problem)
        return 1
    data = comm.bcast(data, root=0)

    say(f"HODDIES:        {hoddies['file']} ({hoddies['version']})")
    say(f"Halos:          {halo_array}")
    particle_info = json.loads(Path(str(particle_file) + ".json").read_text()) if particle_file else None
    say(f"Satellites:     {args.satellites}" + (f", {particle_file} ({particle_info['particles']:,} particles, "
                                                 f"fraction {particle_info['particle_fraction']})" if particle_file else ""))
    say(f"Data vector:    {data['data_vec'].size} values ({' + '.join(args.statistics)})"
        + (f", wp: {data['rp'].size} rp bins {data['edges_rp'][0]:.4g}-{data['edges_rp'][-1]:.4g} Mpc/h" if "wp" in data else "")
        + (f", xi_ell: {data['s'].size} s bins {data['edges_s'][0]:.4g}-{data['edges_s'][-1]:.4g} Mpc/h" if "s" in data else "")
        + f"; fit range {args.fit_r_min} to {args.fit_r_max} Mpc/h")

    halos = np.load(halo_array)
    particles = np.load(particle_file) if particle_file else None
    hod = build_hod(hoddies, halos, particles, data, args, verbose=root)
    del particles
    name = fit_name(hod, args)
    fit_dir = hp.hod_fit_dir(*run, name)
    hod.args["fit_param"]["dir_output_fit"] = str(fit_dir)
    if root:
        if fit_dir.exists():
            shutil.rmtree(fit_dir)  # results of an earlier fit of the same thing
            say(f"Replaces:       the earlier results in {fit_dir}")
        (fit_dir / "evaluations").mkdir(parents=True)
    comm.Barrier()
    say(f"Fit:            {name}")
    say(f"Output:         {fit_dir}")
    import logging
    logging.getLogger("hoddies").setLevel(logging.ERROR)  # no per-evaluation messages during the fit
    say(f"{args.tracer} HOD:        {json.dumps(hod.args[args.tracer], default=str)}")
    say(f"Free:           {', '.join(f'{n} in [{lo}, {hi}]' for n, (lo, hi) in zip(hod.name_params, hod.priors))}")
    say(f"H(z)/H0:        {hod.cosmo.efunc(args.redshift):.6f} at z = {args.redshift}")

    if args.start == "random":
        x0 = hod.get_init_params(seed=args.start_seed)  # also seeds numpy's generator, used by the minimizer
    else:
        x0 = np.array([hod.args[args.tracer][name] for name in args.priors])
        np.random.seed(args.start_seed)
        if not hoddies["values_in_boundaries"](np.array([x0]), hod.priors):
            raise SystemExit(f"The starting point of the parameter file {x0} is outside the priors")
    say(f"Start:          {dict(zip(hod.name_params, np.round(x0, 5)))}")

    if root:
        settings = {
            "fit_name": name, "run": hp.simulation_name(args.target, args.model, args.realization),
            "target": args.target, "model": args.model, "realization": args.realization,
            "redshift": args.redshift, "boxsize": args.boxsize, "halo_catalog": str(halo_array), "satellites": args.satellites,
            "particles": {"file": str(particle_file), "fraction": particle_info["particle_fraction"],
                          "seed": particle_info["particle_seed"], "n": particle_info["particles"]} if particle_file else None,
            "cosmology": {"h": args.h, "Omega_cdm": args.omega_cdm, "Omega_b": args.omega_b, "T_cmb": args.tcmb,
                          "N_eff": args.neff, "n_s": args.n_s, "efunc": float(hod.cosmo.efunc(args.redshift))},
            "hod_param_file": args.hod_param_file.name, "tracer": args.tracer,
            "hod_parameters": hod.args[args.tracer], "priors": args.priors,
            "data": {"statistics": args.statistics, "rppi": str(args.data_rppi), "smu": str(args.data_smu),
                     "correlation": args.correlation, "correlation_statistics": args.correlation_statistics,
                     "pimax": args.pimax, "rp_rebin": args.rp_rebin,
                     "rp_range": args.rp_range, "fit_r_range": [args.fit_r_min, args.fit_r_max], "ells": args.ells, "s_rebin": args.s_rebin,
                     "s_range": args.s_range, "n_data": int(data["data_vec"].size)},
            "model": {"los": args.los, "rsd": True, "mock_seed": args.mock_seed, "model_seeds": args.model_seeds,
                      "reference_param_file": str(args.reference_param_file)},
            "minimizer": {"package": "stochopy", "method": args.method, "maxiter": args.maxiter,
                          "popsize": args.popsize, "xtol": args.xtol, "start": args.start,
                          "start_seed": args.start_seed, "x0": dict(zip(hod.name_params, x0.tolist()))},
            "resources": {"mpi_ranks": size, "threads_per_rank": hod.nthreads},
            "code": {"hicola": f"{hp.HICOLA_ROOT} ({git_version(hp.HICOLA_ROOT)})",
                     "hoddies": f"{hoddies['file']} ({hoddies['version']})"},
            "job_log": args.job_log, "started": time.strftime("%Y-%m-%d %H:%M:%S"),
        }
        (fit_dir / "settings.json").write_text(json.dumps(settings, indent=2, default=str))
        shutil.copy2(args.hod_param_file, fit_dir / args.hod_param_file.name)
        if args.reference_param_file != "none":
            shutil.copy2(args.reference_param_file, fit_dir / Path(args.reference_param_file).name)
        np.savez(fit_dir / "data.npz", **data)
        iterations = open(fit_dir / "iterations.txt", "w", buffering=1)
        iterations.write("# iteration n_evaluations time_s chi2 " + " ".join(hod.name_params) + "\n")

    def callback(population, state):
        if root:
            elapsed = time.time() - start_time
            iterations.write(f"{state.nit} {state.nfev} {elapsed:.1f} {state.fun:.8g} "
                             + " ".join(f"{v:.8g}" for v in state.x) + "\n")
            print(f"Iteration {state.nit:4d}: {state.nfev:6d} evaluations, best chi2 of the iteration "
                  f"{state.fun:.4f}  ({elapsed / 60:.1f} min)", flush=True)

    objective = Objective(hod, hoddies, args.mock_seed, fit_dir / "evaluations" / f"rank_{rank:03d}.txt", start_time)
    say(f"Minimizer:      stochopy {args.method}, maxiter {args.maxiter}, popsize {args.popsize}, xtol {args.xtol}, "
        f"{size} MPI ranks x {hod.nthreads} threads")
    say("================================================================")
    options = {"maxiter": args.maxiter, "popsize": args.popsize, "xtol": args.xtol, "workers": size, "backend": "mpi"}
    result = minimize(objective, bounds=hod.priors, x0=x0, method=args.method, options=options, callback=callback)
    objective.record.close()
    failed = comm.allreduce(objective.failed, op=MPI.SUM)
    # Best evaluation of all ranks: the minimizer returns the best point of its last iteration,
    # which can be worse than one found earlier
    best_chi2, best_x = min(comm.allgather((objective.best_chi2, objective.best_x)), key=lambda item: item[0])

    if root:
        fit_time = time.time() - start_time
        say("================================================================")
        say(f"Minimizer:      {result.message} after {result.nit} iterations, {result.nfev} evaluations "
            f"({failed} failed), {fit_time / 3600:.2f} h; best of the last iteration: chi2 {result.fun:.3f}")
        if best_x is None:
            say("No evaluation gave a finite chi2: no best fit")
            return 1
        model = best_fit_model(hod, hoddies, best_x, args)
        np.savez(fit_dir / "best_fit_model.npz", **model)
        n_data, n_free = int(data["data_vec"].size), len(hod.name_params)
        from scipy import stats as scipy_stats
        best = {
            "parameters": dict(zip(hod.name_params, best_x.tolist())),
            "priors": dict(zip(hod.name_params, [list(p) for p in hod.priors])),
            "chi2": float(best_chi2), "chi2_recomputed": float(model["chi2"][0]),
            "n_data": n_data, "n_free": n_free, "dof": n_data - n_free,
            "chi2_per_dof": float(best_chi2) / (n_data - n_free),
            "p_value": float(scipy_stats.chi2.sf(best_chi2, n_data - n_free)),
            "n_galaxies": int(model["n_gal"][0]), "number_density": float(model["n_gal"][0] / args.boxsize**3),
            "satellite_fraction": float(1 - model["n_cen"][0] / model["n_gal"][0]),
            "sat_nfw_fraction": float(model["sat_nfw_fraction"][0]) if model["sat_nfw_fraction"].size else None,
            "minimizer": {"success": bool(result.success), "status": int(result.status), "message": result.message,
                          "nit": int(result.nit), "nfev": int(result.nfev), "failed_evaluations": int(failed),
                          "last_iteration_best": {"chi2": float(result.fun),
                                                  "parameters": dict(zip(hod.name_params, np.asarray(result.x).tolist()))}},
            "time_h": fit_time / 3600, "finished": time.strftime("%Y-%m-%d %H:%M:%S"),
        }
        if args.reference_param_file != "none":
            fitted_hod = dict(hod.args[args.tracer])
            reference = reference_hod(args.reference_param_file, args.tracer)
            hod.args[args.tracer].update(reference)
            ref_model = model_statistics(hod, hoddies, args)
            np.savez(fit_dir / "reference_model.npz", **ref_model)
            differs = {key: [fitted_hod.get(key), value] for key, value in reference.items()
                       if key not in args.priors and not same_hod_value(key, fitted_hod.get(key), value)}
            best["reference"] = {
                "param_file": str(args.reference_param_file),
                "parameters": {f"{name}_{args.tracer}": reference.get(name) for name in args.priors},
                "other_differences": differs,
                "chi2": float(ref_model["chi2"][0]), "chi2_per_dof": float(ref_model["chi2"][0]) / (n_data - n_free),
                "n_galaxies": int(ref_model["n_gal"][0]), "number_density": float(ref_model["n_gal"][0] / args.boxsize**3),
                "satellite_fraction": float(1 - ref_model["n_cen"][0] / ref_model["n_gal"][0]),
            }
            say(f"Reference:      {Path(args.reference_param_file).name}: chi2 {ref_model['chi2'][0]:.3f} "
                f"(extra seeds {np.round(ref_model['chi2'][1:], 1).tolist()})"
                + (f"; other HOD values differ from the fit's: {differs}" if differs else ""))
        (fit_dir / "best_fit.json").write_text(json.dumps(best, indent=2, default=str))
        best_yaml = write_best_fit_yaml(fit_dir)
        say(f"Best fit:       chi2 {best_chi2:.3f} for {n_data} - {n_free} = {n_data - n_free} dof "
            f"(recomputed {model['chi2'][0]:.3f}; extra seeds {np.round(model['chi2'][1:], 1).tolist()})")
        say(f"Parameters:     {json.dumps(best['parameters'])}")
        say(f"Galaxies:       {best['n_galaxies']:,}, n = {best['number_density']:.3e} (h/Mpc)^3, "
            f"satellite fraction {best['satellite_fraction']:.3f}"
            + (f", of which {best['sat_nfw_fraction']:.1%} on the NFW profile (too few stored particles)"
               if best["sat_nfw_fraction"] is not None else ""))
        say(f"Saved:          {fit_dir}")
        say(f"For the mocks:  HOD_PARAM_FILE=\"{best_yaml}\" and SATELLITES=\"{args.satellites}\" in jobs/hod_mock.sh")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
