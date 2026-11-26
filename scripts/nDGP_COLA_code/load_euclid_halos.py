"""Load the CosmoHub download (Parquet or FITS) and split it by simulation and snapshot."""
import pandas as pd
from astropy.table import Table

FN = "elephant_halos.parquet"        # or elephant_halos.fits
df = pd.read_parquet(FN) if FN.endswith(".parquet") else Table.read(FN).to_pandas()

for (sim, snap), g in df.groupby(["simulation", "snapshot_id"]):
    print(sim, snap, len(g))
# m200c: bound particles only; m200c_all: all particles inside r200c (closer to a plain SO mass, use it next to FoF)
