"""Turn a zeldovich-PLT density field into the white-noise file read by Hi-COLA's
ic_random_field_type = "read_white_noise" (see PIPELINE.md).

zeldovich-PLT (the AbacusSummit IC code) run with a flat P(k) (white_noise.pow), ZD_Pk_norm = 0,
no PLT corrections and ZD_qdensity = 2 writes delta(x) = c * W(x): the white noise of its seed,
times a constant. Its phases are resolution-invariant, so the same seed at another NP gives the
same modes as the reference simulation. This script

- reads that density (raw float32, C order). zeldovich-PLT's file has the axis order of the
  AbacusSummit ic_dens_N*.asdf grids, axis 0 = x (--input-order xyz): the 576^3 test run with it
  reproduces Abacus's IC density with r(k) = 1 on large scales, while a transposed field gives
  r ~ 0 and the same P(k) (checked November 2026). zyx transposes the input to axis 0 = x,
- measures c^2 = <|delta(k)|^2> over the non-zero modes, with Hi-COLA's FFT normalization
  (forward transform divided by N^3; k = 0 and the Nyquist planes left out, as Hi-COLA does),
- writes W = delta / c (raw float32, C order, axis 0 = x), so that <|W(k)|^2> = 1,
- prints <|W(k)|^2> in k-shells (flat for white noise) and the fraction of zero modes.

    python scripts/make_white_noise.py --input density2048 --nmesh 2048 --input-order xyz \
        --output white_noise_N2048.f32 --threads 64
"""
from __future__ import annotations

import argparse
import time

import numpy as np
import scipy.fft as sfft


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", required=True, help="zeldovich-PLT density file (raw float32)")
    parser.add_argument("--nmesh", type=int, required=True, help="grid size N of the file (ppd)")
    parser.add_argument("--input-order", choices=("zyx", "xyz"), required=True,
                        help="axis order of the input: xyz for zeldovich-PLT density files (axis 0 = x)")
    parser.add_argument("--output", required=True, help="white-noise file to write (raw float32, axis 0 = x)")
    parser.add_argument("--threads", type=int, required=True, help="threads for the FFT")
    args = parser.parse_args()

    n = args.nmesh
    t0 = time.time()
    delta = np.fromfile(args.input, dtype=np.float32)
    if delta.size != n**3:
        raise SystemExit(f"{args.input}: {delta.size} values, expected {n}^3 = {n**3}")
    delta = delta.reshape(n, n, n)
    if args.input_order == "zyx":
        delta = np.ascontiguousarray(delta.transpose(2, 1, 0))  # -> [x][y][z]
    print(f"read {args.input} ({n}^3, order {args.input_order}) in {time.time() - t0:.0f} s", flush=True)

    dk = sfft.rfftn(delta, workers=args.threads)
    dk /= np.float32(n**3)                                        # Hi-COLA's forward normalization
    power = (dk.real.astype(np.float64)**2 + dk.imag.astype(np.float64)**2)

    used = np.ones(dk.shape, dtype=bool)                          # modes Hi-COLA keeps
    used[0, 0, 0] = False
    used[n // 2, :, :] = False
    used[:, n // 2, :] = False
    used[:, :, n // 2] = False
    nonzero = used & (power > 0)
    c2 = power[nonzero].mean()
    print(f"modes kept: {used.sum()}, non-zero: {nonzero.sum()} ({nonzero.sum() / used.sum():.4f}), "
          f"c^2 = <|delta(k)|^2> = {c2:.8e}", flush=True)

    # Flatness check: <|W(k)|^2> in shells of |k| (in units of the fundamental mode)
    kx = np.fft.fftfreq(n, d=1.0 / n)
    kz = np.fft.rfftfreq(n, d=1.0 / n)
    kmag = np.sqrt(kx[:, None, None]**2 + kx[None, :, None]**2 + kz[None, None, :]**2)
    edges = np.array([1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096], dtype=float)
    edges = edges[:np.searchsorted(edges, np.sqrt(3) * n / 2) + 1]   # up to the cube corner
    idx = np.digitize(kmag[nonzero], edges)
    w2 = power[nonzero] / c2
    print("<|W(k)|^2> per shell [n_min, n_max) of |k| / (2 pi / L): n_modes, mean")
    for i in range(1, len(edges)):
        sel = idx == i
        if sel.any():
            print(f"  [{edges[i - 1]:6.0f}, {edges[i]:6.0f}): {sel.sum():12d}  {w2[sel].mean():.5f}")
    del dk, power, used, nonzero, kmag, idx, w2

    white = delta / np.float32(np.sqrt(c2))
    white.astype(np.float32).tofile(args.output)
    print(f"wrote {args.output}: W(x) rms = {white.std():.5f}, total {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
