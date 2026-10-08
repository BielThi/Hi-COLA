# Comparison pipeline

How this fork runs Hi-COLA against reference simulations: Hi-COLA backend → FoF halos
(nbodykit) → HMF / HODDIES → comparison notebooks. Hi-COLA itself (`HiCOLA/`) keeps its
original layout and file names; this file only covers what was added around it.

## Where things live

| what | where | |
|---|---|---|
| code, parameter files, notebooks | this repository (`$HICOLA_ROOT`) | git |
| frontend tables, job logs | `HiCOLA/Output/` (Hi-COLA's own output folder) | not in git |
| simulation outputs | `$HICOLA_RUN_ROOT` = `/scratch/$USER/Hi-COLA_Output` | cleaned up automatically |
| reference data, kept measurements | `$HICOLA_DATA_ROOT` = `../hicola_data` | see its README.md |

The three variables are set in `jobs/hicola_env.sh` and read by the Lua files
(`Backend_parameter_files/common.lua`) and the Python code (`scripts/hicola_paths.py`).
Export one before running to change it, e.g. on another cluster.

## Runs

A run is matched to a reference simulation (`target`) and named by a model and a realization:

```
HiCOLA/Parameter_files/Backend_parameter_files/
├── common.lua                         naming scheme and paths
├── COLA_NM4096/base.lua               1024 Mpc/h, Euclid cosmology, Bart's initial P(k), z = 1.65 ... 1.0
├── COLA_NM4096/<model>_<realization>.lua
├── Abacus_Small/base.lua              500 Mpc/h, AbacusSummit c000, z = 1.1, 0.8, 0.5
├── Abacus_Small/<model>_<realization>.lua
├── standalone/*.lua                   other runs (tests, cubic Galileon, ...), not matched to a reference
└── Hi-COLA_params.lua                 original Hi-COLA template
```

- `model`: `LCDM` or `nDGP_N<rcH0>` (e.g. `nDGP_N1.2` → DGP with rc·H0/c = 1.2)
- `realization`: `s<seed>`, `s<seed>-fixed` (fixed amplitude), `s<seed>-fixed-inv`
  (fixed amplitude, reversed phases)

A run file only sets these two and loads `base.lua`; all runs of a target share every other
setting. To add a run, copy a run file and change the names; overrides go after the `dofile`.

`simulation_name` is `<target>_<model>_<realization>` and the output folder is
`$HICOLA_RUN_ROOT/<target>/<model>/<realization>/`, where Hi-COLA writes its usual files:

```
pofk_<simulation_name>_cb_z<z>.txt          every time step
gravitymodel_<simulation_name>_k*.txt, cosmology_<simulation_name>.txt
snapshot_<simulation_name>_z<z>/
├── gadget_z<z>.*, pofk_z<z>.txt            output redshifts
└── halo_catalog/halos200c_{catalog.npy,hmf.txt,hoddies.npy}   from the halo finder
```

## Running

From the repository, with the `hicola` conda env available (the scripts activate it):

```bash
# backend; extra arguments go to sbatch
jobs/submit_backend.sh HiCOLA/Parameter_files/Backend_parameter_files/Abacus_Small/LCDM_s1234.lua
jobs/submit_backend.sh <run.lua> -N 2 -n 128 --time=1:00:00

# halo finder + HMF + HODDIES array, one job per output redshift (or give one redshift)
jobs/submit_halo_finder.sh <run.lua>
jobs/submit_halo_finder.sh <run.lua> 0.5 --dependency=afterok:<backend job id>

# Abacus_Small reference P(k) (edit the parameters in jobs/power_spectrum_abacus.sh)
jobs/submit_power_spectrum_abacus.sh

# frontend (writes to HiCOLA/Output/)
jobs/frontend.sh <model.ini> <cosmology.ini>
```

Logs go to `HiCOLA/Output/COLA/<target>/<model>/<realization>/logs/`. At the end of the
backend and halo-finder jobs, `scripts/collect_products.py` copies the P(k) and HMF at the
output redshifts to `$HICOLA_DATA_ROOT/hicola/`, which is what to keep since scratch is
cleaned up.

## Notebooks

Notebooks find their files through `scripts/hicola_paths.py`:

```python
import hicola_paths as hp   # after the sys.path setup in the first cell

hp.reference_file('COLA_NM4096', 'LCDM', box, 1.65, 'PS')       # reference measurement
hp.hicola_file('Abacus_Small', 'nDGP_N1.2', 's1234', 0.5, 'PS')  # kept Hi-COLA measurement
hp.pofk_file('Abacus_Small', 'nDGP_N1.2', 's1234', 0.5)          # in the run folder on scratch
hp.halo_catalog_dir(...), hp.hmf_file(...), hp.standalone_run_dir('LCDM_test')
```

The HOD notebooks import HODDIES from `../HODDIES` (`$HODDIES_ROOT`).

## Setting up elsewhere

1. Clone this repository and `HODDIES` next to each other, and copy `hicola_data/` next to them
   (or export `HICOLA_DATA_ROOT`).
2. Create the `hicola` conda env and build `HiCOLA/FML/COLASolver/nbody` (its Makefile points
   to the env's FFTW, GSL and Lua).
3. Export `HICOLA_RUN_ROOT` if the scratch path differs, and adapt the `#SBATCH` lines in
   `jobs/*.sh` (partition, cores per node) to the cluster.
