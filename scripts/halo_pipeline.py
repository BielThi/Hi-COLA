from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Mapping

import numpy as np


HALO_CATALOG_FIELDS = (
    "x",
    "y",
    "z",
    "vx",
    "vy",
    "vz",
    "m_fof",
    "log10_m_fof",
    "r_halo",
    "rs",
    "c",
    "vrms",
    "length",
    "halo_id",
)

HODDIES_DTYPE = [
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


# Halo velocities (vx, vy, vz, vrms) are peculiar velocities in km/s. Hi-COLA's Gadget snapshots
# store v_pec / sqrt(a), so run_halo_catalog.py multiplies by sqrt(a) and records this in the
# catalog metadata; catalogs made before that fix (2026-10-14) lack the key and are refused.
VELOCITY_DEFINITION = "peculiar"


def require_peculiar_velocities(metadata: Mapping[str, Any], catalog_path: str | Path) -> None:
    """Stop if a halo catalog's velocities are not marked as peculiar (made before the fix)."""
    if metadata.get("velocity_definition") != VELOCITY_DEFINITION:
        raise SystemExit(
            f"{catalog_path}: halo velocities are not marked as peculiar velocities; this catalog was "
            "made before the velocity fix of 2026-10-14 (velocities too large by sqrt(1+z)). "
            "Rerun jobs/halo_finder.sh on this run and redshift."
        )


def require_fields(catalog: np.ndarray, fields: tuple[str, ...] | list[str]) -> None:
    names = catalog.dtype.names or ()
    missing = [field for field in fields if field not in names]
    if missing:
        raise ValueError(f"Missing required catalog fields: {missing}; available fields are {names}")


def validate_halo_catalog(catalog: np.ndarray, *, boxsize: float | None = None) -> None:
    require_fields(catalog, HALO_CATALOG_FIELDS)
    if np.any(catalog["m_fof"] <= 0.0):
        raise ValueError("Halo catalog contains non-positive m_fof values")
    if not np.allclose(catalog["log10_m_fof"], np.log10(catalog["m_fof"]), rtol=1.0e-12, atol=1.0e-12):
        raise ValueError("log10_m_fof is inconsistent with m_fof")
    if np.any(catalog["r_halo"] <= 0.0):
        raise ValueError("Halo catalog contains non-positive r_halo values")
    if np.any(catalog["rs"] <= 0.0):
        raise ValueError("Halo catalog contains non-positive rs values")
    if not np.allclose(catalog["rs"], catalog["r_halo"] / catalog["c"], rtol=1.0e-12, atol=1.0e-12):
        raise ValueError("rs is inconsistent with r_halo / c")
    if np.any(catalog["vrms"] < 0.0):
        raise ValueError("Halo catalog contains negative vrms values")
    if np.any(catalog["length"] <= 0):
        raise ValueError("Halo catalog contains non-positive FoF lengths")
    if len(np.unique(catalog["halo_id"])) != len(catalog):
        raise ValueError("halo_id values are not unique")
    if boxsize is not None:
        for field in ("x", "y", "z"):
            if np.nanmin(catalog[field]) < 0.0 or np.nanmax(catalog[field]) >= boxsize:
                raise ValueError(f"{field} positions are outside [0, boxsize)")


def halo_catalog_to_hoddies(catalog: np.ndarray, redshift: float) -> np.ndarray:
    """HODDIES halo array of a halo catalog: Mh = the FoF mass, Rh = r_halo and c derived from it by
    nbodykit (as a virial mass, the catalog's radius_mdef). HODDIES places satellites in comoving
    coordinates (as with Abacus' comoving r98), while nbodykit's radius is physical: Rh and Rs are
    converted to comoving kpc/h."""
    require_fields(
        catalog,
        ["x", "y", "z", "vx", "vy", "vz", "m_fof", "log10_m_fof", "r_halo", "rs", "c", "vrms", "halo_id"],
    )
    arr = np.empty(len(catalog), dtype=HODDIES_DTYPE)
    for field in ["x", "y", "z", "vx", "vy", "vz", "c", "halo_id"]:
        arr[field] = catalog[field]
    arr["Mh"] = catalog["m_fof"]
    arr["log10_Mh"] = catalog["log10_m_fof"]
    arr["Rh"] = catalog["r_halo"] * (1 + redshift)
    arr["Rs"] = catalog["rs"] * (1 + redshift)
    arr["Vrms"] = catalog["vrms"]
    return arr


def load_catalog_metadata(catalog_path: str | Path) -> dict[str, Any]:
    metadata_path = Path(str(catalog_path) + ".json")
    if not metadata_path.exists():
        return {}
    with metadata_path.open(encoding="utf-8") as handle:
        return json.load(handle)


def geometric_edges_from_centers(centers: np.ndarray) -> np.ndarray:
    centers = np.asarray(centers, dtype="f8")
    log_centers = np.log10(centers)
    log_edges = np.empty(len(centers) + 1, dtype="f8")
    log_edges[1:-1] = 0.5 * (log_centers[:-1] + log_centers[1:])
    log_edges[0] = log_centers[0] - 0.5 * (log_centers[1] - log_centers[0])
    log_edges[-1] = log_centers[-1] + 0.5 * (log_centers[-1] - log_centers[-2])
    return 10.0**log_edges


def compute_hmf(masses: np.ndarray, volume: float, mass_grid: np.ndarray) -> dict[str, np.ndarray | int]:
    masses = np.asarray(masses, dtype="f8")
    masses = masses[np.isfinite(masses) & (masses > 0.0)]

    thresholds = np.asarray(mass_grid, dtype="f8")
    sorted_masses = np.sort(masses)
    counts_gt = len(sorted_masses) - np.searchsorted(sorted_masses, thresholds, side="right")
    n_gt = counts_gt / volume
    sigma_n_gt = np.sqrt(counts_gt) / volume

    mass_edges = geometric_edges_from_centers(thresholds)
    counts_bin, _ = np.histogram(masses, bins=mass_edges)
    dlog10m = np.diff(np.log10(mass_edges))
    dn_dlog10m = counts_bin / (volume * dlog10m)
    sigma_dn_dlog10m = np.sqrt(counts_bin) / (volume * dlog10m)

    return {
        "mass": thresholds,
        "n_gt": n_gt,
        "sigma_n_gt_poisson": sigma_n_gt,
        "dn_dlog10m": dn_dlog10m,
        "sigma_dn_dlog10m_poisson": sigma_dn_dlog10m,
        "counts_gt": counts_gt,
        "counts_bin": counts_bin,
        "n_halos": len(masses),
    }


def save_hmf_table(hmf: Mapping[str, Any], output_path: str | Path, metadata: Mapping[str, Any], mass_column: str) -> None:
    """HMF table and its .json; mass_column names the mass definition in the header, e.g. M_fof."""
    output_path = Path(output_path)
    table = np.column_stack(
        [
            hmf["mass"],
            hmf["n_gt"],
            hmf["sigma_n_gt_poisson"],
            hmf["dn_dlog10m"],
            hmf["sigma_dn_dlog10m_poisson"],
            hmf["counts_gt"],
            hmf["counts_bin"],
        ]
    )
    header = (
        f"{mass_column}[M_sun/h] n(>M)[h^3/Mpc^3] sigma_n(>M)_poisson[h^3/Mpc^3] "
        "dn/dlog10M[h^3/Mpc^3] sigma_dn/dlog10M_poisson[h^3/Mpc^3] counts(>M) counts_bin"
    )
    np.savetxt(output_path, table, header=header)
    with open(output_path.with_suffix(output_path.suffix + ".json"), "w", encoding="utf-8") as handle:
        json.dump(dict(metadata), handle, indent=2, sort_keys=True)
