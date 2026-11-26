"""Checks of zeldovich-PLT fields regenerated for an AbacusSummit simulation (jobs/white_noise_abacus.sh).

All inputs are zeldovich-PLT density files (raw float32, C order, axis 0 = x, ZD_qdensity = 2):

check A:  --abacus (regenerated at N_A with Abacus's own parameters and P(k)) against Abacus's stored IC
          density --stored (ic_dens_N<N_A>.asdf): max and rms difference, as is and transposed (2,1,0).
          "as is" should agree to ~1e-4 of the rms (float32 / platform differences), "transposed" not at all.
check AB: --abacus against --white-low (same seed and N_A, flat P(k)): same phases, i.e. A(k) / B(k) real
          and positive for every mode non-zero in both.
check BC: --white-low against --white-high (same seed, N_C > N_A): identical complex values on the
          modes common to both grids (|n_i| < N_A / 2), i.e. zeldovich-PLT's resolution invariance.

Fourier convention: forward FFT divided by N^3 (Hi-COLA's).  k = 0 is excluded by the non-zero masks.
Reference values for highbase_c000_ph100 (N_A = 1152, N_C = 2048, Perlmutter, 2026-11-16):
  A: rms(diff)/rms = 5.2e-5;  AB: median |arg(A/B)| = 1.9e-7, Re(A/B) > 0 for 0.9983 of the modes;
  BC: rms|C - B| / rms|B| = 2.9e-7.

    python scripts/validate_abacus_white_noise.py --stored ic_dens_N576.asdf --abacus A/density576 \
        --white-low B/density576 --white-high C/density1024 --n-low 576 --n-high 1024 --threads 72
"""
from __future__ import annotations

import argparse

import numpy as np
import scipy.fft as sfft


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--stored", required=True, help="Abacus IC density ic_dens_N<n-low>.asdf")
    parser.add_argument("--abacus", required=True, help="zeldovich-PLT density, Abacus parameters, N = n-low")
    parser.add_argument("--white-low", required=True, help="zeldovich-PLT density, flat P(k), N = n-low")
    parser.add_argument("--white-high", required=True, help="zeldovich-PLT density, flat P(k), N = n-high")
    parser.add_argument("--n-low", type=int, required=True)
    parser.add_argument("--n-high", type=int, required=True)
    parser.add_argument("--threads", type=int, required=True)
    args = parser.parse_args()
    n1, n2 = args.n_low, args.n_high

    def raw(path: str, n: int) -> np.ndarray:
        a = np.fromfile(path, dtype=np.float32)
        if a.size != n**3:
            raise SystemExit(f"{path}: {a.size} values, expected {n}^3")
        return a.reshape(n, n, n)

    def fourier(field: np.ndarray) -> np.ndarray:
        return sfft.rfftn(field, workers=args.threads) / np.float32(field.shape[0]**3)

    # check A
    import asdf
    import abacusnbody.data  # noqa: F401  (registers the blosc compressor used by the Abacus files)
    with asdf.open(args.stored) as af:
        stored = np.array(af["data"]["density"], dtype=np.float32)
    if stored.shape != (n1,) * 3:
        raise SystemExit(f"{args.stored}: shape {stored.shape}, expected {n1}^3")
    a = raw(args.abacus, n1)
    for name, arr in (("raw as is", a), ("raw transposed (2,1,0)", a.transpose(2, 1, 0))):
        d = arr - stored
        print(f"A: {name:24s} rms(stored) = {stored.std():.6e}  max|diff| = {np.abs(d).max():.3e}  "
              f"rms(diff)/rms = {d.std() / stored.std():.3e}", flush=True)
    del stored, d

    # check AB
    ak = fourier(a)
    bk = fourier(raw(args.white_low, n1))
    ok = (np.abs(ak) > 0) & (np.abs(bk) > 0)
    ratio = ak[ok] / bk[ok]
    phase = np.angle(ratio)
    print(f"AB: {ok.sum()} modes non-zero in both, {(~ok).sum()} not;  |arg(A/B)|: median "
          f"{np.median(np.abs(phase)):.2e}, max {np.abs(phase).max():.2e}, fraction > 1e-3: "
          f"{(np.abs(phase) > 1e-3).mean():.2e};  Re(A/B) > 0 for {(ratio.real > 0).mean():.6f} of them", flush=True)
    del a, ak, ratio, phase, ok

    # check BC
    ck = fourier(raw(args.white_high, n2))
    h = n1 // 2                                   # common modes |n_i| < n1/2 (Nyquist of n1 left out)
    lo1, hi1 = slice(0, h), slice(n1 - h + 1, None)
    lo2, hi2 = slice(0, h), slice(n2 - h + 1, None)
    num = den = maxrel = 0.0
    for s0 in ((lo1, lo2), (hi1, hi2)):
        for s1 in ((lo1, lo2), (hi1, hi2)):
            b = bk[s0[0], s1[0], :h]
            c = ck[s0[1], s1[1], :h]
            num += float(np.sum(np.abs(c - b)**2))
            den += float(np.sum(np.abs(b)**2))
            nz = np.abs(b) > 0
            maxrel = max(maxrel, float(np.max(np.abs(c[nz] - b[nz]) / np.abs(b[nz]))))
    print(f"BC: common modes |n_i| < {h}: rms|C - B| / rms|B| = {np.sqrt(num / den):.3e}, "
          f"max relative difference {maxrel:.3e}", flush=True)


if __name__ == "__main__":
    main()
