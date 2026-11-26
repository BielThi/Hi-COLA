"""Checks on the Abacus IC files (ic_dens / ic_disp, N576):
(i) displacement convention: delta(k) = -i k . psi(k), with psi[..., c] along axis c;
(ii) alignment of the IC density with the z = 0.5 CompaSO halos (cross-correlation r(k) at low k)."""
import asdf, numpy as np, scipy.fft as sfft
import abacusnbody.data            # registers the blosc decompressor for asdf
from abacusnbody.data.compaso_halo_catalog import CompaSOHaloCatalog

IC = "/global/cfs/cdirs/desi/public/cosmosim/AbacusSummit/ic/AbacusSummit_highbase_c000_ph100"
CFS = "/global/cfs/cdirs/desi/public/cosmosim/AbacusSummit/AbacusSummit_highbase_c000_ph100"
L, N = 1000.0, 576

with asdf.open(f"{IC}/ic_disp_N576.asdf") as af:
    psi = np.array(af["data"]["displacements"], dtype=np.float32)      # (N,N,N,3), 2.3 GB
with asdf.open(f"{IC}/ic_dens_N576.asdf") as af:
    d = np.array(af["data"]["density"], dtype=np.float32)
print("rms psi per component:", psi.reshape(-1, 3).std(axis=0))

k1 = (2*np.pi*np.fft.fftfreq(N, d=L/N)).astype(np.float32)
kr = (2*np.pi*np.fft.rfftfreq(N, d=L/N)).astype(np.float32)
kax = [k1[:, None, None], k1[None, :, None], kr[None, None, :]]
dk = sfft.rfftn(d, workers=16)
m = np.ones(dk.shape, bool); m[0, 0, 0] = False
m[N//2] = False; m[:, N//2] = False; m[:, :, -1] = False          # drop k=0 and Nyquist planes

# (i) displacement vs density
for name, comps in {"psi[...,c] -> axis c": (0, 1, 2), "reversed components": (2, 1, 0)}.items():
    div = sum(1j * kax[ax] * sfft.rfftn(np.ascontiguousarray(psi[..., c]), workers=16)
              for ax, c in enumerate(comps))
    r = div[m] / dk[m]
    print(f"{name:24s} median Re = {np.median(r.real):+.4f}  median Im = {np.median(r.imag):+.4f}  scatter = {r.real.std():.2e}")

# (ii) IC density vs z = 0.5 halos
cat = CompaSOHaloCatalog(f"{CFS}/halos/z0.500", fields=["N", "x_L2com"])
x = np.asarray(cat.halos["x_L2com"])[np.asarray(cat.halos["N"]) >= 100]
H, _ = np.histogramdd((x + L/2) % L, bins=N, range=[(0, L)] * 3)       # Abacus positions are in [-L/2, L/2)
hk = sfft.rfftn((H / H.mean() - 1).astype(np.float32), workers=16)
kk = np.sqrt(kax[0]**2 + kax[1]**2 + kax[2]**2)
bins = np.arange(0.5, 15) * 2*np.pi/L                                    # k < ~0.09 h/Mpc
for name, field in {"as is": dk,
                    "transposed": sfft.rfftn(np.ascontiguousarray(d.transpose(2, 1, 0)), workers=16)}.items():
    idx = np.digitize(kk.ravel(), bins)
    cross = np.bincount(idx, (field * hk.conj()).real.ravel(), len(bins) + 1)[1:-1]
    pd_ = np.bincount(idx, (np.abs(field)**2).ravel(), len(bins) + 1)[1:-1]
    ph = np.bincount(idx, (np.abs(hk)**2).ravel(), len(bins) + 1)[1:-1]
    print(name, "r(k) =", np.round(cross / np.sqrt(pd_ * ph), 2))
