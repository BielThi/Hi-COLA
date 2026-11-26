"""Read the precomputed Abacus matter P(k) at z = 0.5 (full pack9 and 10% A+B subsample)
and compare them after subtracting shot noise (none is subtracted in the files)."""
from astropy.table import Table

BASE = "/global/cfs/cdirs/desi/public/cosmosim/AbacusSummit/power/AbacusSummit_highbase_c000_ph100/z0.500"

def read(fn):
    try:
        return Table.read(fn, format="ascii.ecsv")
    except Exception:
        return Table.read(fn, format="ascii.csv", comment="#")

full = read(f"{BASE}/pack9/power_nfft2048.csv")
sub = read(f"{BASE}/AB/power_nfft2048.csv")
print(full.meta); print(full.colnames); print(full[:3])

kcol, pcol = "k", "power"            # replace with the column names printed above
L, N = 1000.0, 3456
sn_full = L**3 / N**3                # 0.024 (Mpc/h)^3
sn_ab = L**3 / (0.10 * N**3)         # 0.24  (Mpc/h)^3
ratio = (sub[pcol] - sn_ab) / (full[pcol] - sn_full)
for k, r in zip(full[kcol][::20], ratio[::20]):
    print(f"k = {k:.3f}  P_AB / P_full = {r:.4f}")
# Do not use the faulty_..._zspace multipoles: RSD along x, multipoles along z.
