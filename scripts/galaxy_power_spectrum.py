"""Redshift-space power spectrum multipoles of a galaxy mock, measured as in the EMC/acm pipeline.

compute_spectrum and the redshift-space step follow acm (github.com/epaillas/acm,
scripts/emc/measurements/abacus_small.py: compute_spectrum and get_hod_positions), so the
results compare directly with the EMC AbacusHOD measurements: jaxpower, TSC mass assignment,
interlacing 3, window compensation, k-bins of width k_step, shot noise from compute_fkp2_shotnoise.
If acm changes its estimator, this copy does not follow.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Mapping

import numpy as np


def redshift_space_positions(x, y, z, vx, vy, vz, boxsize: float, redshift: float, efunc: float, los: str) -> np.ndarray:
    """Positions shifted along the line of sight by v / (a H), wrapped into [-L/2, L/2).

    x, y, z in Mpc/h in [0, L) (or any periodic image), v in peculiar km/s,
    efunc = H(z)/H0 of the simulation's cosmology.
    """
    positions = np.c_[x, y, z].astype("f8")
    velocities = np.c_[vx, vy, vz].astype("f8")
    axis = "xyz".index(los)
    hubble = 100.0 * efunc  # km/s / (Mpc/h)
    scale_factor = 1.0 / (1.0 + redshift)
    positions[:, axis] += velocities[:, axis] / (hubble * scale_factor)
    return np.mod(positions + boxsize / 2.0, boxsize) - boxsize / 2.0


def compute_spectrum(positions: np.ndarray, boxsize: float, meshsize: int, ells, los: str, k_step: float):
    """Power spectrum multipoles of positions in [-L/2, L/2) with jaxpower (as acm's compute_spectrum)."""
    import jax
    from jax import config

    config.update("jax_enable_x64", True)
    from jaxpower import BinMesh2SpectrumPoles, MeshAttrs, ParticleField, compute_fkp2_shotnoise, compute_mesh2_spectrum

    mattrs = MeshAttrs(boxsize=boxsize, boxcenter=0.0, meshsize=meshsize)
    data = ParticleField(positions, attrs=mattrs, exchange=True, backend="jax")
    mesh = data.paint(resampler="tsc", interlacing=3, compensate=True, out="real")
    mean = mesh.mean()
    mesh = mesh - mean
    bin = BinMesh2SpectrumPoles(mesh.attrs, edges={"step": k_step}, ells=tuple(ells))
    jitted_compute_mesh2_spectrum = jax.jit(compute_mesh2_spectrum, static_argnames=["los"], donate_argnums=[0])
    spectrum = jitted_compute_mesh2_spectrum(mesh, bin=bin, los=los)
    num_shotnoise = compute_fkp2_shotnoise(data, bin=bin)
    spectrum = spectrum.clone(norm=[pole.values("norm") * mean**2 for pole in spectrum], num_shotnoise=num_shotnoise)
    jax.block_until_ready(spectrum)
    return spectrum


def save_spectrum(spectrum, ells, output_path: str | Path, metadata: Mapping[str, Any]) -> None:
    """Write <output>.txt (k, N_modes, P_ell shot-noise subtracted), <output>.h5 (jaxpower/lsstypes,
    the format of the EMC measurements) and <output>.txt.json (metadata)."""
    output_path = Path(output_path)
    poles = [spectrum.get(ell) for ell in ells]
    k = np.asarray(poles[0].coords("k"))
    columns = [k, np.asarray(poles[0].values("nmodes"))] + [np.asarray(pole.value()) for pole in poles]
    shotnoise = np.asarray(poles[0].values("num_shotnoise")) / np.asarray(poles[0].values("norm"))
    header = (
        "k[h/Mpc] N_modes " + " ".join(f"P{ell}[(Mpc/h)^3]" for ell in ells)
        + f"   (shot noise subtracted from P0: {np.nanmean(shotnoise):.6e} (Mpc/h)^3)"
    )
    np.savetxt(output_path, np.column_stack(columns), header=header)
    spectrum.write(str(output_path.with_suffix(".h5")))
    with open(str(output_path) + ".json", "w", encoding="utf-8") as handle:
        json.dump(dict(metadata), handle, indent=2, sort_keys=True, default=str)

