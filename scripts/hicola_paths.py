"""Paths and naming shared by the Hi-COLA pipeline scripts, jobs and notebooks.

Three roots, each overridable with an environment variable (see jobs/hicola_env.sh):

    HICOLA_ROOT        this repository
    HICOLA_RUN_ROOT    simulation outputs (large, regenerable)  default: /scratch/$USER/Hi-COLA_Output
    HICOLA_DATA_ROOT   reference data + kept Hi-COLA products   default: <HICOLA_ROOT>/../hicola_data

A run is identified by (target, model, realization), e.g. ("Abacus_Small", "nDGP_N1.2", "s1234-fixed"):

    target       reference simulation the run is matched to: COLA_NM4096 | Abacus_Small
    model        LCDM | nDGP_N<rcH0>
    realization  s<seed>[-fixed[-inv]]

Its Lua file is HiCOLA/Parameter_files/Backend_parameter_files/<target>/<model>_<realization>.lua,
a complete Hi-COLA parameter file with simulation_name = "<target>_<model>_<realization>" and
output_folder = "<run root>/<target>/<model>/<realization>/" (checked by jobs/run_info.lua).
Hi-COLA writes its own file names there.

Reference measurements live in $HICOLA_DATA_ROOT/reference/<sim>/<model>/Box<n>/ and the kept
Hi-COLA measurements in $HICOLA_DATA_ROOT/hicola/<target>/<model>/<realization>/, both named
<prefix>_z<z>_<kind>.txt with kind PS (matter power spectrum) or HMF (halo mass function).
"""

from __future__ import annotations

import getpass
import importlib
import os
from pathlib import Path

HICOLA_ROOT = Path(os.environ.get("HICOLA_ROOT") or Path(__file__).resolve().parents[1])
RUN_ROOT = Path(os.environ.get("HICOLA_RUN_ROOT") or f"/scratch/{getpass.getuser()}/Hi-COLA_Output")
DATA_ROOT = Path(os.environ.get("HICOLA_DATA_ROOT") or HICOLA_ROOT.parent / "hicola_data")
HODDIES_ROOT = Path(os.environ.get("HODDIES_ROOT") or HICOLA_ROOT.parent / "HODDIES")

# Hi-COLA's own output folder: frontend tables, and the logs of each run
OUTPUT_DIR = HICOLA_ROOT / "HiCOLA" / "Output"
BACKEND_PARAM_DIR = HICOLA_ROOT / "HiCOLA" / "Parameter_files" / "Backend_parameter_files"
FRONTEND_PARAM_DIR = HICOLA_ROOT / "HiCOLA" / "Parameter_files" / "Frontend_parameter_files"
NOTEBOOK_DIR = HICOLA_ROOT / "Notebooks"

KINDS = ("PS", "HMF")

# Analysis settings per reference simulation. The simulation settings themselves
# (box size, cosmology, output redshifts) are in the run files, Backend_parameter_files/<target>/*.lua.
TARGETS = {
    "COLA_NM4096": {
        "boxsize": 1024.0,
        "n_boxes": 10,
        "reference_redshifts": (1.65, 1.4, 1.2, 1.0),
        "cosmology": "euclid",
        "halo_nmin": 30,
    },
    "Abacus_Small": {
        "boxsize": 500.0,
        "n_boxes": 1,
        "reference_redshifts": (0.5,),
        "cosmology": "abacus",
        "halo_nmin": 50,
    },
}


def zlabel(redshift: float) -> str:
    """Redshift label used everywhere (Hi-COLA's own convention): 1.65 -> 'z1.650'."""
    return f"z{float(redshift):.3f}"


# ---------------------------------------------------------------- Hi-COLA runs


def simulation_name(target: str, model: str, realization: str) -> str:
    return f"{target}_{model}_{realization}"


def run_file(target: str, model: str, realization: str) -> Path:
    return BACKEND_PARAM_DIR / target / f"{model}_{realization}.lua"


def run_dir(target: str, model: str, realization: str) -> Path:
    return RUN_ROOT / target / model / realization


def snapshot_dir(target: str, model: str, realization: str, redshift: float) -> Path:
    name = simulation_name(target, model, realization)
    return run_dir(target, model, realization) / f"snapshot_{name}_{zlabel(redshift)}"


def snapshot_glob(target: str, model: str, realization: str, redshift: float) -> str:
    return str(snapshot_dir(target, model, realization, redshift) / f"gadget_{zlabel(redshift)}.*")


def pofk_file(target: str, model: str, realization: str, redshift: float) -> Path:
    """P(k) Hi-COLA measures on the snapshot at an output redshift."""
    return snapshot_dir(target, model, realization, redshift) / f"pofk_{zlabel(redshift)}.txt"


def pofk_cb_file(target: str, model: str, realization: str, redshift: float) -> Path:
    """P_cb(k) Hi-COLA writes at every time step, with the linear P(k) alongside."""
    name = simulation_name(target, model, realization)
    return run_dir(target, model, realization) / f"pofk_{name}_cb_{zlabel(redshift)}.txt"


def halo_catalog_dir(target: str, model: str, realization: str, redshift: float) -> Path:
    return snapshot_dir(target, model, realization, redshift) / "halo_catalog"


def halo_catalog_file(target: str, model: str, realization: str, redshift: float) -> Path:
    return halo_catalog_dir(target, model, realization, redshift) / "halos200c_catalog.npy"


def hoddies_catalog_file(target: str, model: str, realization: str, redshift: float) -> Path:
    return halo_catalog_dir(target, model, realization, redshift) / "halos200c_hoddies.npy"


def hmf_file(target: str, model: str, realization: str, redshift: float) -> Path:
    return halo_catalog_dir(target, model, realization, redshift) / "halos200c_hmf.txt"


def log_dir(target: str, model: str, realization: str) -> Path:
    return OUTPUT_DIR / "COLA" / target / model / realization / "logs"


def standalone_run_dir(name: str) -> Path:
    """Output folder of Backend_parameter_files/standalone/<name>.lua."""
    return RUN_ROOT / "standalone" / name


# ---------------------------------------------------------------- kept Hi-COLA products


def hicola_dir(target: str, model: str, realization: str) -> Path:
    return DATA_ROOT / "hicola" / target / model / realization


def hicola_file(target: str, model: str, realization: str, redshift: float, kind: str) -> Path:
    _check_kind(kind)
    name = f"HiCOLA_{simulation_name(target, model, realization)}_{zlabel(redshift)}_{kind}.txt"
    return hicola_dir(target, model, realization) / name


# ---------------------------------------------------------------- reference simulations


def reference_dir(sim: str, model: str, box: int) -> Path:
    return DATA_ROOT / "reference" / sim / model / f"Box{box}"


def reference_file(sim: str, model: str, box: int, redshift: float, kind: str) -> Path:
    _check_kind(kind)
    return reference_dir(sim, model, box) / f"{sim}_{model}_Box{box}_{zlabel(redshift)}_{kind}.txt"


def reference_raw_dir(sim: str, model: str | None = None, box: int | None = None) -> Path:
    """Original data as downloaded: per box for Abacus_Small, one tree for COLA_NM4096."""
    if model is None:
        return DATA_ROOT / "reference" / sim / "raw"
    return reference_dir(sim, model, box) / "raw"


def abacus_halo_dir(redshift: float = 0.5, sim_name: str = "AbacusSummit_small_c000_ph3000") -> Path:
    """CompaSOHaloCatalog redshift folder of the Abacus_Small box."""
    raw = reference_raw_dir("Abacus_Small", "LCDM", 1)
    return raw / "small" / sim_name / "halos" / zlabel(redshift)


# ---------------------------------------------------------------- misc


def reference_cosmology(target: str):
    """Reference cosmology module (euclid_/abacus_reference_cosmology) used for a target."""
    return importlib.import_module(f"{TARGETS[target]['cosmology']}_reference_cosmology")


def _check_kind(kind: str) -> None:
    if kind not in KINDS:
        raise ValueError(f"kind must be one of {KINDS}, got {kind!r}")
