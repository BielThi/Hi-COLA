"""Read the AbacusSummit highbase_c000_ph100 z = 0.5 data with abacusutils.
Needs: pip install abacusutils (already in the cosmodesi environment on NERSC)."""
import glob
from abacusnbody.data.read_abacus import read_asdf
from abacusnbody.data.compaso_halo_catalog import CompaSOHaloCatalog

SLICE = "/pscratch/sd/b/bielthi/abacus_highbase"   # where retrieve_abacus_z05.sh puts the full slice
CFS = "/global/cfs/cdirs/desi/public/cosmosim/AbacusSummit/AbacusSummit_highbase_c000_ph100"
L = 1000.0

# 1) One slab of the full particle slice
files = sorted(glob.glob(f"{SLICE}/**/*_pack9/*.asdf", recursive=True))
cat = read_asdf(files[0], load=("pos", "vel"))     # add "pid" if the matching _pid files exist
print(cat.meta["BoxSize"], cat.meta["Redshift"], len(cat))
# pos in Mpc/h, box centred on 0 (-L/2 .. +L/2); vel in km/s

# 2) Loop over all slabs, e.g. to paint a mesh or write sub-cubes for FoF
for fn in files:
    cat = read_asdf(fn, load=("pos",))
    x = cat["pos"] + L / 2          # shift to [0, L) if your code expects that
    # ... paint x onto your mesh / write out a sub-cube here ...
    del cat

# 3) CompaSO halos with the 10% A subsample (public CFS copy)
hcat = CompaSOHaloCatalog(f"{CFS}/halos/z0.500", fields=["N", "x_L2com", "v_L2com"],
                          subsamples=dict(A=True, rv=True))
halos = hcat.halos          # astropy Table
