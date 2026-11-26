"""Write the initial P(k) of an AbacusSummit simulation in Hi-COLA's ic_type_of_input = "powerspectrum" format.

zeldovich-PLT (the AbacusSummit IC code) reads CLASS_power, the linear P(k) at ZD_Pk_file_redshift
(z = 1 for AbacusSummit), and multiplies it by ZD_Pk_sigma_ratio^2, a scale-independent factor
(src/power_spectrum.cpp, PowerSpectrum::Normalize, ZD_Pk_norm > 0). The result is the P(k) of its
initial conditions at InitialRedshift (z = 99). This script writes that P(k) (k grid: see Interpolation), so that Hi-COLA run with

    ic_input_redshift = ic_initial_redshift = <InitialRedshift>

uses it unchanged (its growth-factor ratio D(z_ini) / D(z_input) is then 1).

Interpolation: zeldovich-PLT evaluates P(k) with a natural cubic spline in (ln k, ln P); Hi-COLA uses a
natural cubic spline (GSL cspline) in (k, P). On the sparse CLASS_power k points the two differ by up
to 0.17% (k ~ 1.8 h/Mpc). The output is therefore zeldovich-PLT's spline evaluated on the CLASS_power
points plus --n-dense log-spaced points over the same k range; with --n-dense 20000, Hi-COLA's spline
of the output matches zeldovich-PLT's P(k) to < 2e-9 (12-digit output) for 2 pi / 1000 < k < 11 h/Mpc (checked
November 2026).

- CLASS_power: the file of the AbacusSummit Cosmologies repository (abacus_cosm<NNN>/CLASS_power),
  as named by ZD_Pk_filename in the IC headers. Columns k [h/Mpc], P [(Mpc/h)^3]; '#' lines skipped.
- --sigma-ratio: ZD_Pk_sigma_ratio of the IC header (ic_dens_N*.asdf, header).
- Output: k range of CLASS_power (no extrapolation), one header line, then "k P" (Hi-COLA skips exactly one line).

    python scripts/make_abacus_ic_pk.py --class-power CLASS_power --sigma-ratio 0.02113950155394198 \
        --file-redshift 1.0 --initial-redshift 99.0 --n-dense 20000 --output abacus_cosm000_ic_pk_z99.000.dat
"""
from __future__ import annotations

import argparse

import numpy as np
from scipy.interpolate import CubicSpline


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--class-power", required=True, help="AbacusSummit CLASS_power file")
    parser.add_argument("--sigma-ratio", type=float, required=True, help="ZD_Pk_sigma_ratio of the IC header")
    parser.add_argument("--file-redshift", type=float, required=True, help="ZD_Pk_file_redshift (header)")
    parser.add_argument("--initial-redshift", type=float, required=True, help="InitialRedshift (header)")
    parser.add_argument("--n-dense", type=int, required=True, help="log-spaced points added to the k grid")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    k, p = np.loadtxt(args.class_power, comments="#", usecols=(0, 1), unpack=True)
    if np.any(np.diff(k) <= 0) or np.any(p <= 0):
        raise SystemExit(f"{args.class_power}: k not strictly increasing or P <= 0")
    loglog = CubicSpline(np.log(k), np.log(p), bc_type="natural")    # zeldovich-PLT's interpolation
    k = np.union1d(k, np.geomspace(k[0], k[-1], args.n_dense)[1:-1])
    p_ini = np.exp(loglog(np.log(k))) * args.sigma_ratio**2
    header = (f"k[h/Mpc]  P(k)[(Mpc/h)^3] at z={args.initial_redshift:g} -- AbacusSummit IC P(k): {args.class_power} "
              f"(z={args.file_redshift:g}) x ZD_Pk_sigma_ratio^2, log-log natural cubic spline on {args.n_dense} extra points, ratio = {args.sigma_ratio!r} "
              f"(scripts/make_abacus_ic_pk.py)")
    np.savetxt(args.output, np.column_stack([k, p_ini]), fmt="%.12e", header=header, comments="# ")
    print(f"wrote {args.output}: {k.size} points, k = {k[0]:.4e} .. {k[-1]:.4e} h/Mpc, "
          f"factor {args.sigma_ratio**2:.10e}")


if __name__ == "__main__":
    main()
