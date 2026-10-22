"""Loading and combining measurements for the comparison notebooks.

All functions return plain numpy arrays: x (k in h/Mpc, or M in Msun/h) and y (P(k) in
(Mpc/h)^3, or n(>M) in h^3/Mpc^3). Only the first two columns of a file are read, which is
the common format of every PS/HMF file in the pipeline (see hicola_data/README.md).

    import analysis as an
    k, P_mean, P_std = an.load_reference('COLA_NM4096', 'LCDM', 1.65, 'PS')
    k, P = an.load_hicola('COLA_NM4096', 'LCDM', ['s1234-fixed', 's1234-fixed-inv'], 1.65, 'PS')
"""

from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np

import hicola_paths as hp

PAIRED_FIXED = ("s1234-fixed", "s1234-fixed-inv")


def load_table(path: str | Path) -> tuple[np.ndarray, np.ndarray]:
    """First two columns of a PS or HMF text file."""
    table = np.loadtxt(path, comments="#", ndmin=2)
    return table[:, 0], table[:, 1]


# ---------------------------------------------------------------- reference simulations


def load_reference_boxes(sim: str, model: str, redshift: float, kind: str, boxes=None):
    """x and the stack (n_boxes, n_x) of a reference measurement over the given boxes.

    boxes: list of box numbers; None means every Box<n> folder of that sim and model.
    """
    if boxes is None:
        folders = (hp.DATA_ROOT / "reference" / sim / model).glob("Box*")
        boxes = sorted(int(f.name[3:]) for f in folders if f.name[3:].isdigit())
        if not boxes:
            raise FileNotFoundError(f"No Box<n> folders in {hp.DATA_ROOT / 'reference' / sim / model}")
    tables = [load_table(hp.reference_file(sim, model, box, redshift, kind)) for box in boxes]
    x = tables[0][0]
    for box, (xb, _) in zip(boxes, tables):
        if not np.allclose(xb, x):
            raise ValueError(f"{sim} {model} Box{box}: x grid differs from Box{list(boxes)[0]}")
    return x, np.vstack([y for _, y in tables])


def load_reference(sim: str, model: str, redshift: float, kind: str, boxes=None):
    """x, mean and standard deviation of a reference measurement over its boxes."""
    x, stack = load_reference_boxes(sim, model, redshift, kind, boxes)
    return x, stack.mean(axis=0), stack.std(axis=0)


def reference_ratio(sim: str, model_num: str, model_den: str, redshift: float, kind: str, boxes=None):
    """x, mean and standard deviation over boxes of the per-box ratio model_num / model_den.

    Bins where the denominator is not positive are ignored (NaN when no box is left).
    """
    x, den = load_reference_boxes(sim, model_den, redshift, kind, boxes)
    x_num, num = load_reference_boxes(sim, model_num, redshift, kind, boxes)
    num = np.vstack([np.interp(x, x_num, row) for row in num])
    with np.errstate(divide="ignore", invalid="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)  # all-NaN bins
        per_box = np.where(den > 0, num / den, np.nan)
        return x, np.nanmean(per_box, axis=0), np.nanstd(per_box, axis=0)


# ---------------------------------------------------------------- Hi-COLA runs


def hicola_measurement_file(target: str, model: str, realization: str, redshift: float, kind: str) -> Path:
    """The run's own file on scratch if it is still there, otherwise the kept copy."""
    run_file = (hp.pofk_file if kind == "PS" else hp.hmf_file)(target, model, realization, redshift)
    kept_file = hp.hicola_file(target, model, realization, redshift, kind)
    for path in (run_file, kept_file):
        if path.exists():
            return path
    raise FileNotFoundError(
        f"No {kind} for {target}/{model}/{realization} at z={redshift:.3f}: neither {run_file} nor "
        f"{kept_file} exists. Run jobs/backend.sh"
        + (" and jobs/halo_finder.sh" if kind == "HMF" else "")
        + f" on {hp.run_file(target, model, realization)}."
    )


def load_hicola(target: str, model: str, realizations, redshift: float, kind: str):
    """x and the mean over realizations of a Hi-COLA measurement (e.g. a paired-fixed pair)."""
    if isinstance(realizations, str):
        realizations = [realizations]
    tables = [load_table(hicola_measurement_file(target, model, r, redshift, kind)) for r in realizations]
    x = tables[0][0]
    for realization, (xr, _) in zip(realizations, tables):
        if not np.allclose(xr, x):
            raise ValueError(f"{target}/{model}/{realization}: x grid differs from {realizations[0]}")
    return x, np.mean([y for _, y in tables], axis=0)


def load_hod_spectrum(target: str, model: str, realization: str, redshift: float, name: str, dk: float,
                      statistic: str = "PS"):
    """Redshift-space multipoles of a HOD mock (jobs/hod_mock.sh), rebinned to bins of width dk.

    statistic: 'PS' (power spectrum) or 'MPS' (marked power spectrum).

    Reads the run folder on scratch, or the kept copy in hicola_data. Returns k (mode-weighted
    mean of each bin) and {ell: P_ell}, averaged over the fine bins with weights N_modes; the
    k = 0 mode and empty bins are left out.
    """
    run_file = hp.hod_spectrum_file(target, model, realization, redshift, name, statistic)
    kept_file = hp.hicola_hod_spectrum_file(target, model, realization, redshift, name, statistic)
    path = run_file if run_file.exists() else kept_file
    if not path.exists():
        raise FileNotFoundError(
            f"No HOD {statistic} '{name}' for {target}/{model}/{realization} at z={redshift:.3f}: neither "
            f"{run_file} nor {kept_file} exists. Run jobs/hod_mock.sh on {hp.run_file(target, model, realization)}."
        )
    with open(path) as handle:
        header = handle.readline()
    ells = [int(word[1:].split("[")[0]) for word in header.split() if word.startswith("P") and "[" in word]
    table = np.loadtxt(path, comments="#", ndmin=2)
    k, nmodes = table[:, 0], table[:, 1]
    use = (k > 0) & (nmodes > 0)
    index = np.floor(k[use] / dk).astype(int)
    weights = np.bincount(index, weights=nmodes[use])
    filled = weights > 0
    k_binned = np.bincount(index, weights=k[use] * nmodes[use])[filled] / weights[filled]
    poles = {
        ell: np.bincount(index, weights=table[use, 2 + i] * nmodes[use])[filled] / weights[filled]
        for i, ell in enumerate(ells)
    }
    return k_binned, poles


def load_hod_spectra(target: str, model: str, realization: str, redshift: float, name: str, seeds, dk: float,
                     statistic: str = "PS"):
    """HOD P(k) ('PS') or marked P(k) ('MPS') multipoles of several seeds: k and {ell: array (n_seeds, n_k)}.

    name contains {seed} (as NAME in jobs/hod_mock.sh, e.g. 'LRG_seed{seed}'); every seed must exist.
    """
    tables = [load_hod_spectrum(target, model, realization, redshift, name.format(seed=seed), dk, statistic)
              for seed in seeds]
    k = tables[0][0]
    for seed, (k_seed, _) in zip(seeds, tables):
        if not (len(k_seed) == len(k) and np.allclose(k_seed, k)):
            raise ValueError(f"{target}/{model}/{realization} seed {seed}: k-bins differ from seed {seeds[0]}")
    return k, {ell: np.vstack([poles[ell] for _, poles in tables]) for ell in tables[0][1]}


def load_hod_fit(target: str, model: str, realization: str, redshift: float, fit_name: str) -> dict:
    """Results of a HOD fit (jobs/hod_fit.sh) from its folder, hp.hod_fit_dir.

    Returns a dict: 'dir', 'settings' (dict), 'data' (dict of arrays), 'evaluations' (structured
    array of every evaluation, all ranks, in order of time), 'iterations' (structured array, one row
    per iteration of the minimizer), and, once the fit has finished, 'best_fit' (dict) and 'model'
    (dict of arrays); None while it runs or if it stopped before the end.
    """
    folder = hp.hod_fit_dir(target, model, realization, redshift, fit_name)
    if not (folder / "settings.json").exists():
        raise FileNotFoundError(f"No HOD fit '{fit_name}' for {target}/{model}/{realization} at z={redshift:.3f} "
                                f"({folder}/settings.json missing). Run jobs/hod_fit.sh with FIT_NAME={fit_name}.")
    import json

    def read_table(path):
        table = np.genfromtxt(path, names=True, ndmin=1)
        return table if table.size else None

    rank_tables = [read_table(path) for path in sorted((folder / "evaluations").glob("rank_*.txt"))]
    rank_tables = [table for table in rank_tables if table is not None]
    evaluations = np.concatenate(rank_tables) if rank_tables else None
    if evaluations is not None:
        evaluations = evaluations[np.argsort(evaluations["time_s"], kind="stable")]
    iterations = read_table(folder / "iterations.txt") if (folder / "iterations.txt").exists() else None
    finished = (folder / "best_fit.json").exists()
    return {
        "dir": folder,
        "settings": json.loads((folder / "settings.json").read_text()),
        "data": dict(np.load(folder / "data.npz")),
        "evaluations": evaluations,
        "iterations": iterations,
        "best_fit": json.loads((folder / "best_fit.json").read_text()) if finished else None,
        "model": dict(np.load(folder / "best_fit_model.npz")) if finished else None,
    }


# ---------------------------------------------------------------- combining


def interp_to(x_new, x, y):
    """y interpolated onto x_new, NaN outside the range of x."""
    return np.interp(x_new, x, y, left=np.nan, right=np.nan)


def ratio(x_num, y_num, x_den, y_den):
    """y_num / y_den on the x_den grid, NaN where y_den is not positive or out of range."""
    num = interp_to(x_den, x_num, y_num)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(y_den > 0, num / y_den, np.nan)


def within(x, x_min=None, x_max=None):
    """Boolean mask of x_min <= x <= x_max (either bound optional)."""
    mask = np.ones_like(x, dtype=bool)
    if x_min is not None:
        mask &= x >= x_min
    if x_max is not None:
        mask &= x <= x_max
    return mask

