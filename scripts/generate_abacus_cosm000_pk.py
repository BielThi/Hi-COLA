#!/usr/bin/env python3
"""Generate the linear P_cb(k, z=0) input file for AbacusSummit's abacus_cosm000 cosmology.

Why this exists: LCDM_Abacus.lua and nDGP1dot0_Abacus.lua set the initial conditions for
the COLA runs from an external linear power spectrum file (ic_type_of_input =
"powerspectrum", ic_sigma8_normalization = false -- i.e. whatever amplitude is in the file
is used as-is, with no rescaling). That file used to be
LCDM_1306dot3219_matterpower_z0.000.dat, a generic example bundled with Hi-COLA with an
implied sigma8(z=0) ~= 0.917 -- nowhere close to AbacusSummit c000's real sigma8_cb =
0.8114 -- which was the dominant cause of the amplitude offset seen between our COLA runs
and the Abacus box's measured matter P(k).

This script computes the real thing instead: P_cb(k, z=0) for abacus_cosm000 via CLASS
(through cosmoprimo's AbacusSummit fiducial cosmology, which already carries the exact
cosm000 parameters -- h, omega_b, omega_cdm, omega_ncdm, A_s, n_s -- cross-checked against
both cosmoprimo's abacus_cosmologies.csv and the Abacus box's own CompaSOHaloCatalog
header). "cb" (CDM + baryons, excluding the massive neutrino) is used rather than total
matter because the simulated particles -- in both the Abacus box and this COLA setup --
are cb-only; neutrinos are not particle-sampled.

Output: k[h/Mpc], P(k)[(Mpc/h)^3], two columns, same format as the file it replaces, wide
enough in k (1e-5 to 50 h/Mpc) to comfortably cover the 500 Mpc/h box's fundamental mode
down to well past the ic_nmesh=1024 mesh's Nyquist frequency (~6.4 h/Mpc) with no
extrapolation needed at the edges of the range HiCOLA will actually use.

Run this in the `hicola` conda env (it has cosmoprimo + classy installed):
    conda activate hicola
    python scripts/generate_abacus_cosm000_pk.py
"""

from __future__ import annotations

import numpy as np

import hicola_paths as hp

OUTPUT_PATH = hp.BACKEND_PARAM_DIR / "abacus_cosm000_pk_cb_z0.000.dat"

# Expected sigma8_cb(z=0) for abacus_cosm000, from cosmoprimo's abacus_cosmologies.csv.
# Used only as a sanity check on the file we just wrote (see verify_sigma8 below).
EXPECTED_SIGMA8_CB = 0.811355
SIGMA8_TOLERANCE = 0.01  # absolute; CLASS precision settings can shift this by ~0.001-0.003


def compute_pk_cb(k_hmpc: np.ndarray) -> np.ndarray:
    """Linear P_cb(k, z=0) in (Mpc/h)^3 for abacus_cosm000, from CLASS via cosmoprimo."""
    import cosmoprimo.fiducial as fiducial

    cosmo = fiducial.AbacusSummit(name="000", engine="class")
    print(f"h={cosmo.h}  A_s={cosmo.A_s:.6e}  n_s={cosmo.n_s}")
    print(f"sigma8_cb(z=0) = {cosmo.sigma8_cb:.6f}  (expect ~{EXPECTED_SIGMA8_CB})")

    fo = cosmo.get_fourier()
    pk_cb_interp = fo.pk_interpolator(of="delta_cb", extrap_kmin=1e-7, extrap_kmax=1e3)
    return pk_cb_interp(k_hmpc, z=0.0)


def sigma8_from_pk(k: np.ndarray, pk: np.ndarray, R: float = 8.0) -> float:
    """sigma8 = sqrt( (1/2pi^2) int dk k^2 P(k) W_TH(kR)^2 ), integrated in log-k."""
    from scipy.integrate import simpson

    x = k * R
    W = 3.0 * (np.sin(x) - x * np.cos(x)) / x**3
    integrand = k**3 * pk * W**2
    integral = simpson(integrand, x=np.log(k))
    return np.sqrt(integral / (2.0 * np.pi**2))


def main() -> None:
    k = np.geomspace(1e-5, 50.0, 1200)
    pk = compute_pk_cb(k)

    # HiCOLA's ic_type_of_input="powerspectrum" reader (FML/FileUtils/FileUtils.cpp)
    # hardcodes nheaderlines=1 and just skips that many raw lines with no comment-char
    # awareness -- so this MUST be exactly one line, or the 2nd header line gets fed to
    # std::stod and crashes every rank at startup (ask me how I know).
    header = "k[h/Mpc]  P_cb(k)[(Mpc/h)^3]  -- AbacusSummit abacus_cosm000, z=0, CLASS via cosmoprimo (scripts/generate_abacus_cosm000_pk.py)"
    np.savetxt(OUTPUT_PATH, np.column_stack([k, pk]), fmt="%.8e", header=header)
    print(f"Wrote {OUTPUT_PATH}  ({len(k)} rows, k in [{k[0]:.2e}, {k[-1]:.2e}] h/Mpc)")

    # Sanity check: re-derive sigma8 from the file we just wrote and compare to the
    # known target. This is exactly the check that caught the old file's sigma8~0.917 bug.
    s8 = sigma8_from_pk(k, pk)
    print(f"sigma8_cb(z=0) recovered from written file: {s8:.6f}")
    if abs(s8 - EXPECTED_SIGMA8_CB) > SIGMA8_TOLERANCE:
        raise RuntimeError(
            f"sigma8 mismatch: got {s8:.6f}, expected ~{EXPECTED_SIGMA8_CB} "
            f"(tolerance {SIGMA8_TOLERANCE}). Something is off -- do not point the .lua "
            f"files at this output yet."
        )
    print("OK: matches the expected AbacusSummit c000 sigma8_cb within tolerance.")


if __name__ == "__main__":
    main()
