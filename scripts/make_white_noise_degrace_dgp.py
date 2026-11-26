"""White noise W(x) of a DEGRACE-DGP box (MG-GLAM, regenerated ICs) for Hi-COLA's
ic_random_field_type = "read_white_noise" (see PIPELINE.md and scripts/make_white_noise.py).

GLAM's ICs (boxNNNN/ICs/, read with DEGRACE-DGP-loader/load_degrace_dgp.py) are Zel'dovich at z_init:
2048^3 particles on a lattice, particle (i, j, k) = index i + N (j + N k), at
q = (i, j, k) L/N + q0 + psi(q). The tests of November 2026 found q0 = 0.125 Mpc/h on each axis
(a quarter of the lattice spacing; Hi-COLA's lattice is at (i, j, k) L/N) and v = s psi with
correlation 1.000000 (s = 541.722 km/s per Mpc/h for the loader's velocities). This script

1. reads all particles and puts psi on the N^3 lattice, taken from the velocities, psi = v / s
   (s fitted on the first record; velocities are stored with float32 relative precision, while
   positions near L = 1024 carry float32 rounding of ~6e-5 Mpc/h, about 1e-3 of rms(psi));
   it checks q0 = mean(pos - q) against --lattice-offset and reports rms(pos - q - q0 - v/s);
2. computes the linear density delta_G(k) = -i k . psi(k) (Zel'dovich: delta = -div psi) on GLAM's
   lattice, shifts it to Hi-COLA's lattice, delta(k) = delta_G(k) exp(-i k . q0) (q0 = --lattice-offset
   on each axis), and whitens it, W(k) = delta(k) / sqrt(P(|k|)), with P the --pk table interpolated
   as Hi-COLA does (natural cubic spline in (k, P)), so that Hi-COLA's W(k) sqrt(P(k)) with the same
   file gives back GLAM's delta(k) exactly, whatever the table's tail;
3. sets to zero the k = 0 mode, the Nyquist planes and the modes outside GLAM's IC cutoff sphere
   |k| >= --cutoff-n k_f (GLAM's ICs have no power there: 1024 for these boxes, checked: the script
   stops if the power outside is not negligible), so that Hi-COLA sees them as zero modes (also with
   ic_fix_amplitude, which only rescales non-zero modes), and normalises as make_white_noise.py does,
   c^2 = <|W(k) / N^3|^2> (Hi-COLA's forward FFT is divided by N^3), over the modes with
   k < --norm-kmax only: below the end of GLAM's own table (5.2052 h/Mpc) the P(k) shape is GLAM's,
   above it the table's tail is an extrapolation and W is not exactly white (about 8% low near the
   cutoff for the power-law tail of degrace_dgp_PkTable_z0.000.dat), which would bias c^2;
4. writes W(x) as raw float32, C order, axis 0 = x (what Hi-COLA reads), plus a .json;
5. prints <|W(k)|^2> per k-shell (flat below --norm-kmax);
6. with --validate-nmesh M > 0: CIC density of all particles (absolute positions) on an M^3 grid,
   compared mode by mode with W(k) sqrt(P) for k < 1 h/Mpc: cross-correlation r(k) (axis order and
   sign) and the shift per axis fitted to arg(delta_W delta_CIC*) = a . k (0 if the field is on
   Hi-COLA's lattice, +q0 without the shift). M must be incommensurate with the particle lattice
   (e.g. odd): with M = 512 (4 particles per cell per axis) the discrete CIC deposit itself biases
   this fit (+0.04 Mpc/h for q0 = 0.125 and rms psi = 0.1 Mpc/h; checked on a synthetic lattice,
   where an incommensurate grid gives the true offset to 1e-3).

Self-tests of the phase-shift and divergence conventions run first on a 16^3 grid.

    python scripts/make_white_noise_degrace_dgp.py --ic-dir /pscratch/sd/b/bielthi/DEGRACE-DGP/box0002/ICs \
        --loader ~/hicola_data/DEGRACE-DGP/DEGRACE-DGP-loader \
        --pk HiCOLA/Parameter_files/Backend_parameter_files/Initial_power_spectra/degrace_dgp_PkTable_z0.000.dat \
        --boxsize 1024 --npart-1d 2048 --lattice-offset 0.125 --cutoff-n 1024 --norm-kmax 5.2052 --threads 128 --read-threads 16 \
        --validate-nmesh 495 --validate-procs 28 --output white_noise_N2048.f32
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from multiprocessing import Pool

import numpy as np
import scipy.fft as sfft
from scipy.interpolate import CubicSpline

L = None   # loader module, imported in main()


def selftest() -> None:
    n, box, q0 = 16, 8.0, 0.37
    x = np.arange(n) * box / n
    m = np.array([1, 2, 3])
    kv = 2 * np.pi * m / box
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    # shift: a field sampled at x + q0, moved to x
    f_glam = np.cos(kv[0] * (X + q0) + kv[1] * (Y + q0) + kv[2] * (Z + q0))
    F = np.fft.rfftn(f_glam)
    k1 = 2 * np.pi * np.fft.fftfreq(n, box / n); k3 = 2 * np.pi * np.fft.rfftfreq(n, box / n)
    F *= np.exp(-1j * (k1[:, None, None] + k1[None, :, None] + k3[None, None, :]) * q0)
    f_back = np.fft.irfftn(F, s=(n, n, n))
    target = np.cos(kv[0] * X + kv[1] * Y + kv[2] * Z)
    assert np.abs(f_back - target).max() < 1e-10, "phase-shift convention"
    # divergence: psi = A sin(k.x) k_hat -> -div psi = -A |k| cos(k.x)
    A = 0.3; kn = np.linalg.norm(kv)
    ph = kv[0] * X + kv[1] * Y + kv[2] * Z
    psi = [A * np.sin(ph) * kv[a] / kn for a in range(3)]
    kk = [k1[:, None, None], k1[None, :, None], k3[None, None, :]]
    D = sum(-1j * kk[a] * np.fft.rfftn(psi[a]) for a in range(3))
    assert np.abs(np.fft.irfftn(D, s=(n, n, n)) - (-A * kn * np.cos(ph))).max() < 1e-10, "divergence convention"
    print("self-tests (phase shift, divergence) passed", flush=True)


def read_record(icdir, r, n):
    pos, vel = next(L.iter_ic_particles(icdir, records=[r]))
    idx0 = r * L.NPAGE
    k, rest = divmod(idx0, n * n)
    j0 = rest // n
    nj = L.NPAGE // n
    i = np.arange(n, dtype=np.float64)
    j = np.arange(j0, j0 + nj, dtype=np.float64)
    return k, j0, nj, pos, vel, i, j


def cic_worker(args):
    icdir, recs, m, box = args
    grid = np.zeros(m**3)
    for r in recs:
        pos, _ = next(L.iter_ic_particles(icdir, records=[r]))
        x = pos.astype(np.float64) / (box / m)
        i0 = np.floor(x).astype(np.int64); d = x - i0; i0 %= m; i1 = (i0 + 1) % m
        idx, w = [], []
        for cx in (0, 1):
            for cy in (0, 1):
                for cz in (0, 1):
                    ix = i1[:, 0] if cx else i0[:, 0]; iy = i1[:, 1] if cy else i0[:, 1]; iz = i1[:, 2] if cz else i0[:, 2]
                    wx = d[:, 0] if cx else 1 - d[:, 0]; wy = d[:, 1] if cy else 1 - d[:, 1]; wz = d[:, 2] if cz else 1 - d[:, 2]
                    idx.append((ix * m + iy) * m + iz); w.append(wx * wy * wz)
        grid += np.bincount(np.concatenate(idx), weights=np.concatenate(w), minlength=m**3)
    return grid


def main() -> None:
    global L
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--ic-dir", required=True, help="boxNNNN/ICs folder (PMcrd.DAT, PMcrs*.DAT)")
    parser.add_argument("--loader", required=True, help="DEGRACE-DGP-loader folder")
    parser.add_argument("--pk", required=True, help="P(k) table [k h/Mpc, P (Mpc/h)^3], one header line, for the shape")
    parser.add_argument("--boxsize", type=float, required=True)
    parser.add_argument("--npart-1d", type=int, required=True)
    parser.add_argument("--lattice-offset", type=float, required=True, help="q0 [Mpc/h] of GLAM's lattice on each axis")
    parser.add_argument("--threads", type=int, required=True, help="FFT threads")
    parser.add_argument("--read-threads", type=int, required=True)
    parser.add_argument("--cutoff-n", type=int, required=True, help="GLAM's IC cutoff: modes with |k| >= cutoff-n * k_f are zero")
    parser.add_argument("--norm-kmax", type=float, required=True, help="normalise <|W|^2> = 1 over the modes with k below this [h/Mpc]")
    parser.add_argument("--validate-nmesh", type=int, required=True, help="CIC grid for the check against the particles (0 = no check; odd, incommensurate with the lattice)")
    parser.add_argument("--validate-procs", type=int, required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    sys.path.insert(0, args.loader)
    import load_degrace_dgp as loader
    L = loader

    selftest()
    t0 = time.time()
    n, box, q0 = args.npart_1d, args.boxsize, args.lattice_offset
    h = L.read_ic_header(f"{args.ic_dir}/PMcrd.DAT")
    if int(h["NROW"]) != n or float(h["extras"][99]) != box:
        raise SystemExit(f"header NROW {h['NROW']} / box {h['extras'][99]} differ from --npart-1d / --boxsize")
    nrec = n**3 // L.NPAGE
    dq = box / n

    # slope s from the first record
    k_, j0, nj, pos, vel, i, j = read_record(args.ic_dir, 0, n)
    q = np.stack(np.broadcast_arrays(i[None, :] * dq, j[:, None] * dq, np.full((nj, n), k_ * dq)), -1).reshape(-1, 3)
    dpos = (pos - q + box / 2) % box - box / 2
    dpos_c = dpos - dpos.mean(0)
    s = np.sum(vel * dpos_c) / np.sum(dpos_c**2)
    print(f"v = s psi on record 0: s = {s:.4f} km/s per Mpc/h", flush=True)

    psi = [np.empty((n, n, n), dtype=np.float32) for _ in range(3)]   # [k][j][i] = [z][y][x]
    stats = {"sum_dpos": np.zeros(3), "sum_res2": 0.0, "sum_psi2": 0.0, "n": 0}

    def fill(r):
        k_, j0, nj, pos, vel, i, j = read_record(args.ic_dir, r, n)
        p = vel / np.float32(s)
        for a in range(3):
            psi[a][k_, j0:j0 + nj, :] = p[:, a].reshape(nj, n)
        q = np.stack(np.broadcast_arrays(i[None, :] * dq, j[:, None] * dq, np.full((nj, n), k_ * dq)), -1).reshape(-1, 3)
        dpos = (pos.astype(np.float64) - q + box / 2) % box - box / 2
        res = dpos - q0 - p
        return dpos.sum(0), float((res**2).sum()), float((p.astype(np.float64)**2).sum())

    with ThreadPoolExecutor(args.read_threads) as ex:
        for sd, r2, p2 in ex.map(fill, range(nrec)):
            stats["sum_dpos"] += sd; stats["sum_res2"] += r2; stats["sum_psi2"] += p2; stats["n"] += L.NPAGE
    q0_meas = stats["sum_dpos"] / stats["n"]
    rel_res = np.sqrt(stats["sum_res2"] / stats["sum_psi2"])
    print(f"read {nrec} records in {time.time()-t0:.0f} s; mean(pos - q) = {np.round(q0_meas, 6)} Mpc/h (--lattice-offset {q0}); "
          f"rms(pos - q - q0 - v/s) / rms(v/s) = {rel_res:.2e}", flush=True)
    if np.abs(q0_meas - q0).max() > 1e-3:
        raise SystemExit("measured lattice offset differs from --lattice-offset")

    # delta(k) = -i k.psi(k);  axes [z][y][x], rfft along x
    kf = 2 * np.pi / box
    kz = (kf * np.fft.fftfreq(n, 1.0 / n)).astype(np.float32)
    ky = kz.copy()
    kx = (kf * np.fft.rfftfreq(n, 1.0 / n)).astype(np.float32)
    kvecs = [kx[None, None, :], ky[None, :, None], kz[:, None, None]]   # psi[0] = x, psi[1] = y, psi[2] = z
    delta = None
    for a in range(3):
        F = sfft.rfftn(psi[a], workers=args.threads)
        psi[a] = None
        F *= (-1j * kvecs[a]).astype(np.complex64)
        if delta is None:
            delta = F
        else:
            delta += F
        del F
    print(f"delta(k) done at {time.time()-t0:.0f} s", flush=True)

    kt, pt = np.loadtxt(args.pk, skiprows=1, usecols=(0, 1), unpack=True)
    pspline = CubicSpline(kt, pt, bc_type="natural")     # Hi-COLA's interpolation (GSL natural cspline in (k, P))
    with open(args.pk) as fp:
        pk_header = fp.readline().strip()
    nh = n // 2
    kcut = args.cutoff_n * kf
    shell_edges_n = np.array([1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 768, 896, 960, 1024, 2048, 4096], dtype=float)
    shell_edges_n = shell_edges_n[:np.searchsorted(shell_edges_n, np.sqrt(3) * nh) + 1]
    k_edges = np.array([0, 0.1, 1.0, 3.0, args.norm_kmax, kcut, 20.0])
    sh_n = np.zeros(len(shell_edges_n) + 1); sh_w = np.zeros(len(shell_edges_n) + 1)
    kb_n = np.zeros(len(k_edges) + 1); kb_w = np.zeros(len(k_edges) + 1)
    sum_w2, nnorm, nkept, nzero_in = 0.0, 0, 0, 0
    sum_in, n_in, max_out, max_nyq, sum_nyq, n_nyq = 0.0, 0, 0.0, 0.0, 0.0, 0
    for iz in range(n):
        kmag = np.sqrt(kz[iz]**2 + ky[:, None]**2 + kx[None, :]**2)
        sl = delta[iz]
        p_before = sl.real.astype(np.float64)**2 + sl.imag.astype(np.float64)**2
        inside = (kmag < kcut) & (kmag > 0)
        nyq = np.zeros(kmag.shape, dtype=bool)
        nyq[nh, :] = True; nyq[:, nh] = True
        if iz == nh:
            nyq[...] = True
        sum_in += p_before[inside].sum(); n_in += int(inside.sum())
        outside = ~inside & ~nyq & (kmag > 0)            # Nyquist planes: zeroed below in any case (Hi-COLA too)
        if outside.any():
            max_out = max(max_out, float(p_before[outside].max()))
        max_nyq = max(max_nyq, float(p_before[nyq].max())); sum_nyq += p_before[nyq].sum(); n_nyq += int(nyq.sum())
        kk = np.where(inside, kmag, kf)
        fac = np.exp(-1j * (kz[iz] + ky[:, None] + kx[None, :]) * q0) / np.sqrt(pspline(kk))
        fac[~inside] = 0                                   # GLAM's cutoff sphere and k = 0: exactly zero
        sl *= fac.astype(np.complex64)
        if iz == nh:
            sl[...] = 0
        sl[nh, :] = 0
        sl[:, nh] = 0
        p2 = (sl.real.astype(np.float64)**2 + sl.imag.astype(np.float64)**2) / float(n)**6
        used = p2 > 0
        nkept += int(used.sum()); nzero_in += int((inside & ~used).sum())
        norm = used & (kmag < args.norm_kmax)
        sum_w2 += p2[norm].sum(); nnorm += int(norm.sum())
        b = np.digitize(kmag[used] / kf, shell_edges_n); sh_n += np.bincount(b, minlength=len(sh_n)); sh_w += np.bincount(b, weights=p2[used], minlength=len(sh_w))
        b = np.digitize(kmag[used], k_edges); kb_n += np.bincount(b, minlength=len(kb_n)); kb_w += np.bincount(b, weights=p2[used], minlength=len(kb_w))
    c2 = sum_w2 / nnorm
    out_ratio = max_out / (sum_in / n_in)
    print(f"whitened at {time.time()-t0:.0f} s: {nkept} non-zero modes, {nzero_in} zero modes inside |n| < {args.cutoff_n} "
          f"(Nyquist planes); c^2 = {c2:.8e} from {nnorm} modes with k < {args.norm_kmax}; "
          f"max |delta(k)|^2 outside the cutoff (Nyquist planes excluded) / mean inside = {out_ratio:.2e} (set to 0); "
          f"Nyquist planes (set to 0): mean {sum_nyq / n_nyq / (sum_in / n_in):.2e}, max {max_nyq / (sum_in / n_in):.2e} x mean inside", flush=True)
    if out_ratio > 1e-6:
        raise SystemExit(f"power outside |n| = {args.cutoff_n} is not negligible: --cutoff-n does not match the ICs")
    print("<|W(k)|^2> per shell [n_min, n_max) of |k| / k_f: n_modes, mean")
    for b in range(1, len(shell_edges_n)):
        if sh_n[b] > 0:
            print(f"  [{shell_edges_n[b-1]:6.0f}, {shell_edges_n[b]:6.0f}): {int(sh_n[b]):12d}  {sh_w[b]/sh_n[b]/c2:.5f}")
    print(f"<|W(k)|^2> per k range [h/Mpc] (normalised over k < {args.norm_kmax}; cutoff {kcut:.4f}):")
    for b in range(1, len(k_edges)):
        if kb_n[b] > 0:
            print(f"  [{k_edges[b-1]:7.4f}, {k_edges[b]:7.4f}): {int(kb_n[b]):12d}  {kb_w[b]/kb_n[b]/c2:.5f}")

    # Hi-COLA: W(k) = rfftn(W(x)) / N^3 must equal delta / (N^3 sqrt(c2))  ->  W(x) = irfftn(delta) / sqrt(c2)
    val = None
    if args.validate_nmesh > 0:
        # modes with |n_i| <= nlim (k < 1 h/Mpc) in the CIC grid's layout [x][y][z] (rfft along z),
        # taken from delta in [z][y][x] (rfft along x) with Hermitian symmetry
        nlim = int(1.0 / kf) + 1
        nv = np.arange(-nlim, nlim + 1); nzv = np.arange(0, nlim + 1)
        NX, NY, NZ = np.meshgrid(nv, nv, nzv, indexing="ij")
        pos_x = NX >= 0
        val = np.empty(NX.shape, dtype=np.complex128)
        val[pos_x] = delta[NZ[pos_x] % n, NY[pos_x] % n, NX[pos_x]]
        val[~pos_x] = np.conj(delta[(-NZ[~pos_x]) % n, (-NY[~pos_x]) % n, -NX[~pos_x]])
        val /= float(n)**3 * np.sqrt(c2)
        val_modes = (NX, NY, NZ)
    W = sfft.irfftn(delta, s=(n, n, n), workers=args.threads)
    del delta
    W *= np.float32(1.0 / np.sqrt(c2))
    out = np.ascontiguousarray(W.transpose(2, 1, 0))   # [z][y][x] -> [x][y][z], axis 0 = x
    del W
    out.astype(np.float32, copy=False).tofile(args.output)
    rms = float(out.std())
    del out
    print(f"wrote {args.output}: W(x) rms = {rms:.5f}, at {time.time()-t0:.0f} s", flush=True)
    meta = {"source_ic_dir": args.ic_dir, "nmesh": n, "boxsize_mpc_h": box, "lattice_offset_mpc_h": q0,
            "lattice_offset_measured": q0_meas.tolist(), "psi_from": "velocities, psi = v / s",
            "s_kms_per_mpch": s, "rel_rms_pos_minus_v_over_s": rel_res, "pk_shape_file": args.pk, "pk_header": pk_header,
            "c2_hicola_convention": c2, "norm_kmax_h_mpc": args.norm_kmax, "n_modes_norm": nnorm,
            "cutoff_n": args.cutoff_n, "n_modes_nonzero": nkept, "W_rms": rms, "order": "raw float32, C order, axis 0 = x",
            "made_by": "scripts/make_white_noise_degrace_dgp.py"}

    if val is not None:
        m = args.validate_nmesh
        icdir = args.ic_dir
        chunks = np.array_split(np.arange(nrec), args.validate_procs)
        t1 = time.time()
        with Pool(args.validate_procs) as p:
            grids = p.map(cic_worker, [(icdir, c_, m, box) for c_ in chunks], chunksize=1)
        g = sum(grids); del grids
        dcf = np.fft.rfftn(g.reshape(m, m, m) / g.mean() - 1) / m**3
        NX, NY, NZ = val_modes
        dc = dcf[NX % m, NY % m, NZ]
        KX, KY, KZ = kf * NX, kf * NY, kf * NZ
        km = np.sqrt(KX**2 + KY**2 + KZ**2)
        dw = val * np.sqrt(pspline(np.where(km > 0, km, kf)))
        ok = (km > 0) & (km < 1.0) & (np.abs(val) > 0)
        cross = (dw * np.conj(dc)).real
        edges = np.array([0, 0.05, 0.1, 0.2, 0.4, 0.7, 1.0])
        b = np.digitize(km[ok], edges)
        sx = np.bincount(b, weights=cross[ok]); s1 = np.bincount(b, weights=np.abs(dw[ok])**2); s2 = np.bincount(b, weights=np.abs(dc[ok])**2)
        with np.errstate(invalid="ignore"):                # bin 0 (below the first edge) is empty
            r = sx / np.sqrt(s1 * s2)
        phi = np.angle(dw[ok] * np.conj(dc[ok])); w = np.abs(dw[ok]) * np.abs(dc[ok])
        A = np.c_[KX[ok], KY[ok], KZ[ok]] * np.sqrt(w)[:, None]
        a_fit = np.linalg.lstsq(A, phi * np.sqrt(w), rcond=None)[0]
        print(f"validation (CIC {m}^3 of all particles, {time.time()-t1:.0f} s): r(k) in k bins {edges.tolist()}: {np.round(r[1:], 5).tolist()}; "
              f"fitted shift per axis (x, y, z) = {np.round(a_fit, 5).tolist()} Mpc/h (0 = on Hi-COLA's lattice, +q0 = {q0} without the shift; "
              f"use an M incommensurate with the lattice, e.g. odd: a grid commensurate with it biases this fit)", flush=True)
        meta["validation"] = {"nmesh_cic": m, "k_edges": edges.tolist(), "r": r[1:].tolist(), "shift_fit_mpc_h": a_fit.tolist()}
    with open(args.output + ".json", "w") as fp:
        json.dump(meta, fp, indent=1)
    print(f"done in {time.time()-t0:.0f} s")


if __name__ == "__main__":
    main()
