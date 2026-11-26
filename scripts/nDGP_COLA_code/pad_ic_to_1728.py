"""Zero-pad the Abacus IC density (N1152, z = 99) to a 1728^3 grid for the phase-matched COLA run.
Modes above the 1152 Nyquist frequency (k ~ 3.6-5.4 h/Mpc) are set to zero.
Memory: ~40 GB peak; run on a full CPU node."""
import gc, asdf, numpy as np, scipy.fft as sfft
import abacusnbody.data

IC = "/global/cfs/cdirs/desi/public/cosmosim/AbacusSummit/ic/AbacusSummit_highbase_c000_ph100"
OUT = "/pscratch/sd/b/bielthi/ndgp_cola_validation/ics/delta_c000_ph100_z99_N1728.f32"
N0, N1 = 1152, 1728

with asdf.open(f"{IC}/ic_dens_N1152.asdf") as af:
    d = np.array(af["data"]["density"], dtype=np.float32)
f0 = sfft.rfftn(d, workers=64); del d; gc.collect()              # complex64, ~6 GB

h = N0 // 2                                   # 576 = Nyquist index of the 1152 grid
lo, hi = slice(0, h), slice(-(h - 1), None)   # modes n = 0..575 and -575..-1 (Nyquist dropped)
f1 = np.zeros((N1, N1, N1 // 2 + 1), dtype=np.complex64)         # ~21 GB
for s0 in (lo, hi):
    for s1 in (lo, hi):
        f1[s0, s1, :h] = f0[s0, s1, :h]
f1 *= np.float32((N1 / N0) ** 3)              # keeps the real-space amplitude of delta
del f0; gc.collect()

d1 = sfft.irfftn(f1, s=(N1,) * 3, workers=64).astype(np.float32); del f1; gc.collect()
d1.tofile(OUT)                                # raw float32, C order, axis 0 = x
print("wrote", OUT)
