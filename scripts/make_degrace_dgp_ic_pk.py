"""Write GLAM's input P(k) of the DEGRACE-DGP boxes (PkTable.dat) in Hi-COLA's ic_type_of_input = "powerspectrum" format.

PkTable.dat (DEGRACE-DGP-common) is the linear matter P(k) at z = 0 that GLAM's IC generator read
(k [h/Mpc], P [(Mpc/h)^3], 544 log-spaced points, k = 1e-4 .. 5.205 h/Mpc, five '#' header lines).
Hi-COLA's reader skips exactly one header line, and its spline returns the last tabulated value
beyond the table (a flat P(k) above k_max, without warning). This script writes

- the PkTable.dat points unchanged (12 significant digits; the file has 5),
- then a power law through the last table point, P = P(k_max) (k / k_max)^n, with n the log-log
  slope fitted to the last --n-fit table points (continuous value and slope at k_max; a fit over
  the whole last decade gives n = -2.39 against a local end slope of -2.55 and starts 7% above the
  last point), evaluated on the table's own log spacing from k_max up to --k-extend,
- with one header line saying where the extrapolation starts.

Interpolation: Hi-COLA uses a natural cubic spline in (k, P) on these points (2% spacing in k);
GLAM's own interpolation of PkTable.dat is not documented. The normalisation is the table's own
(sigma8 = 0.8161 from the table vs 0.8159 in its header and Setup.dat; GLAM renormalises to sigma8).

    python scripts/make_degrace_dgp_ic_pk.py \
        --pktable ~/hicola_data/DEGRACE-DGP/DEGRACE-DGP-common/PkTable.dat --n-fit 10 --k-extend 20 \
        --output HiCOLA/Parameter_files/Backend_parameter_files/Initial_power_spectra/degrace_dgp_PkTable_z0.000.dat
"""
from __future__ import annotations

import argparse

import numpy as np


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--pktable", required=True, help="DEGRACE-DGP PkTable.dat")
    parser.add_argument("--n-fit", type=int, required=True, help="number of last table points for the end slope")
    parser.add_argument("--k-extend", type=float, required=True, help="extrapolate the power law up to this k [h/Mpc]")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    k, p = np.loadtxt(args.pktable, comments="#", usecols=(0, 1), unpack=True)
    if np.any(np.diff(k) <= 0) or np.any(p <= 0):
        raise SystemExit(f"{args.pktable}: k not strictly increasing or P <= 0")
    kmax = k[-1]
    fit = np.arange(len(k)) >= len(k) - args.n_fit
    slope, lnA_fit = np.polyfit(np.log(k[fit]), np.log(p[fit]), 1)
    resid = np.log(p[fit]) - (lnA_fit + slope * np.log(k[fit]))
    lnA = np.log(p[-1]) - slope * np.log(kmax)          # anchored at the last table point
    decade = k >= kmax / 10
    slope_decade = np.polyfit(np.log(k[decade]), np.log(p[decade]), 1)[0]
    local = np.diff(np.log(p[-6:])) / np.diff(np.log(k[-6:]))
    dlnk = np.log(k[-1] / k[-2])
    nnew = int(np.floor(np.log(args.k_extend / kmax) / dlnk))
    k_ext = kmax * np.exp(dlnk * np.arange(1, nnew + 1))
    p_ext = np.exp(lnA + slope * np.log(k_ext))

    header = (f"k[h/Mpc]  P(k)[(Mpc/h)^3] at z=0 -- DEGRACE-DGP PkTable.dat (GLAM input, {len(k)} points, k <= {kmax:g}) "
              f"then power law through the last point, P = {np.exp(lnA):.6e} k^{slope:.6f} (slope fitted on the last {fit.sum()} points, k >= {k[fit][0]:.4g}), "
              f"for {kmax:g} < k <= {k_ext[-1]:.4g} ({nnew} points) (scripts/make_degrace_dgp_ic_pk.py)")
    np.savetxt(args.output, np.column_stack([np.r_[k, k_ext], np.r_[p, p_ext]]), fmt="%.12e", header=header, comments="# ")
    print(f"wrote {args.output}: {len(k)} table points + {nnew} extrapolated, k = {k[0]:.4e} .. {k_ext[-1]:.4e} h/Mpc")
    print(f"end slope fitted on the last {fit.sum()} points (k in [{k[fit][0]:.4g}, {kmax:.4g}]): {slope:.5f}, "
          f"max |log residual| {np.abs(resid).max():.2e}; local slopes of the last 5 intervals {np.round(local, 4)}; "
          f"(slope of a fit over the last decade: {slope_decade:.5f})")


if __name__ == "__main__":
    main()
