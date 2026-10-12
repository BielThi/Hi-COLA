# Comparison pipeline

This fork runs Hi-COLA simulations matched to reference simulations and compares their
matter power spectra, halo mass functions and HOD galaxies:

```
frontend (optional)     backend (Hi-COLA)        halo finder (nbodykit)      analysis
.ini → expansion/force  run.lua → snapshots,     FoF → M200c catalog →       notebooks: P(k), HMF,
tables                  P(k) at each output z    HMF + HODDIES halo array    HOD / LRG mocks
```

Hi-COLA itself (`HiCOLA/`) keeps its original layout and file names. Everything described
here was added around it: `jobs/`, `scripts/`, `Notebooks/`, and the run files in
`HiCOLA/Parameter_files/Backend_parameter_files/`.

## 1. Where things live

```
~/                                  (any parent folder works, see section 9)
├── Hi-COLA/                        this repository                          $HICOLA_ROOT
│   ├── HiCOLA/                     Hi-COLA code (frontend, FML backend), original layout
│   │   ├── Output/                 frontend tables, job logs (not in git)
│   │   └── Parameter_files/        frontend .ini and backend .lua files
│   ├── jobs/                       SLURM jobs and their submit scripts
│   ├── scripts/                    Python pipeline code (paths, halo finder, measurements)
│   ├── Notebooks/                  analysis notebooks and HOD parameter files
│   └── PIPELINE.md                 this file
├── HODDIES/                        HOD code, imported by the HOD notebooks     $HODDIES_ROOT
└── hicola_data/                    reference data + kept measurements          $HICOLA_DATA_ROOT
    (see hicola_data/README.md)

/scratch/$USER/Hi-COLA_Output/      simulation outputs (snapshots, halo catalogs) $HICOLA_RUN_ROOT
```

| what | where | backed up? |
|---|---|---|
| code, run files, notebooks | `Hi-COLA/` | git (GitHub `BielThi/Hi-COLA`) |
| frontend tables, job logs | `Hi-COLA/HiCOLA/Output/` | no (regenerable) |
| snapshots, P(k) per time step, halo catalogs | `$HICOLA_RUN_ROOT` on scratch, ~85 GB per run | no: scratch is cleaned up automatically |
| reference simulations, kept P(k)/HMF of each run | `hicola_data/` | home |

The four variables are set in `jobs/hicola_env.sh` and read by the Python code
(`scripts/hicola_paths.py`) with the same defaults; export one before running to change it.
In the jobs, `HICOLA_RUN_ROOT` is taken from the `output_folder` of the run file, so the
`.lua` file decides where a run goes.

## 2. Software setup

**Conda env `hicola`** (Python 3.11). Main packages: `nbodykit` (FoF, P(k)), `mpi4py`,
`abacusutils` (Abacus catalogs), `cosmoprimo` + `pyclass` (CLASS), `numba`, `fitsio`,
`pynbody`, `astropy`, `sympy` + `configobj` (frontend), `jupyter`, `nbstripout`, plus FFTW
(with MPI and threads), GSL, Lua and an MPI compiler, used to build the backend.
The job scripts activate it themselves (`$HICOLA_CONDA_ENV`, default `hicola`).

**Backend binary** `HiCOLA/FML/COLASolver/nbody`, built with `make` in that folder with the
`hicola` env activated: the Makefile takes FFTW, GSL and Lua from the active env
(`$(CONDA_PREFIX)`) and FML from the clone it is in, so every clone builds against its own code.

**HODDIES**: the fork `github.com/BielThi/HODDIES`, cloned next to this repository. The HOD
notebooks import it from there (`$HODDIES_ROOT`), not from the `HODDIES` package installed
in the env, which is a different copy. The fork's `main` (commit `2572a69`) has two
workarounds on top of the original code: numba parallelism turned off (crashes on jed; same
results, slower), and no custom cosmology built from the `cosmo` parameters (that line raised
a `NameError`), so the notebooks set `HOD_obj.cosmo` themselves.

**hi_class**: the env's `classy` was built from `~/hi_class_public`
(`miguelzuma/hi_class_public`, commit `50f447c4`). Nothing in the pipeline imports it any
more (CLASS is used through cosmoprimo/`pyclass`), so it is not needed on a new machine.

**nbstripout** removes notebook outputs from what git stores (the files on disk keep
them). It is registered in `.gitattributes`; each new clone has to run it once:

```bash
conda activate hicola && nbstripout --install
```

## 3. Runs and their names

A run is matched to a reference simulation (the **target**: same box, cosmology and initial
P(k)) and named by a **model** and a **realization**:

| | values |
|---|---|
| target | `COLA_NM4096` (1024 Mpc/h, Euclid cosmology, Bart's initial P(k), z = 1.65, 1.4, 1.2, 1.0) · `Abacus_Small` (500 Mpc/h, AbacusSummit c000, z = 1.1, 0.8, 0.5) |
| model | `LCDM` · `nDGP_N<rcH0>`, e.g. `nDGP_N1.2` = DGP with rc·H0/c = 1.2 (Hi-COLA's built-in DGP) |
| realization | `s<seed>` · `s<seed>-fixed` (fixed amplitude) · `s<seed>-fixed-inv` (fixed amplitude, reversed phases) |

```
HiCOLA/Parameter_files/Backend_parameter_files/
├── Hi-COLA_params.lua            original Hi-COLA template
├── LCDM_1306dot3219_matterpower_z0.000.dat   original Hi-COLA example P(k) (σ8 ≈ 0.92, testing only)
├── Initial_power_spectra/        Bart_InitPk.txt (COLA_NM4096), abacus_cosm000_pk_cb_z0.000.dat (Abacus_Small)
├── COLA_NM4096/<model>_<realization>.lua
├── Abacus_Small/<model>_<realization>.lua
└── standalone/                   runs not matched to a reference, kept as examples:
    ├── cuGal.lua                 cubic Galileon through the frontend (gravity_model = "HiCOLA")
    └── nDGP_action.lua           nDGP from its Horndeski action through the frontend
```

Every run file is a complete Hi-COLA parameter file, like the original template: every
setting is written in it. Within a target, the files differ only in `simulation_name`,
`output_folder`, `gravity_model` (and `gravity_model_dgp_rcH0overc`), `ic_random_seed`,
`ic_fix_amplitude` and `ic_reverse_phases`. Input paths (initial P(k), frontend tables) are
relative to the Hi-COLA folder, where the jobs run `nbody`.

The pipeline finds a run's outputs by its name, so a run file must follow this naming (the
jobs check it and stop with a message otherwise):

| in the run file | value |
|---|---|
| location | `<target>/<model>_<realization>.lua` |
| `simulation_name` | `"<target>_<model>_<realization>"` |
| `output_folder` | `"<run root>/<target>/<model>/<realization>/"`, run root `/scratch/biel/Hi-COLA_Output` on jed |

To add a run, copy a run file, rename it and change those values (and whatever else differs).
Hi-COLA then writes its usual files to the output folder:

```
pofk_<simulation_name>_cb_z<z>.txt         P_cb(k) and linear P(k) at every time step
gravitymodel_<simulation_name>_k*.txt, cosmology_<simulation_name>.txt
snapshot_<simulation_name>_z<z>/           one per output redshift
├── gadget_z<z>.*                          particles (Gadget format)
├── pofk_z<z>.txt                          P(k) of the snapshot: this is the "PS" that is kept
└── halo_catalog/                          written by the halo finder
    ├── halos200c_catalog.npy (+ .json)    FoF halos, M200c, positions, velocities, concentration
    ├── halos200c_hmf.txt (+ .json)        HMF: this is the "HMF" that is kept
    └── halos200c_hoddies.npy (+ .json)    halo array in HODDIES' format
```

Standalone runs write to `<run root>/standalone/<file name>/` and keep their own
`simulation_name`. They need frontend tables in `HiCOLA/Output/` (section 4), which are not
in git.

## 4. Running

Each job script in `jobs/` has a **user parameters** block at the top: edit it, then submit
the script with `sbatch`, from any folder. Paths in the parameter blocks are relative to the
Hi-COLA folder; a run file given on the command line is relative to where you submit. The
`#SBATCH` lines are the production settings on jed (academic partition); edit them in the
script or override them on the command line.

| script | parameters | run with |
|---|---|---|
| `jobs/backend.sh` | `RUN_FILE` (the run's `.lua`) | `sbatch jobs/backend.sh` |
| `jobs/halo_finder.sh` | `RUN_FILE`, `REDSHIFT`, FoF settings (`LINKING_LENGTH`, `NMIN`: 30 for COLA_NM4096, 50 for Abacus_Small so far, `MDEF`) | `sbatch jobs/halo_finder.sh` |
| `jobs/power_spectrum_abacus.sh` | redshift, subsamples, mesh, k-bins of the Abacus reference P(k) | `sbatch jobs/power_spectrum_abacus.sh` |
| `jobs/frontend.sh` | `MODEL_INI`, `COSMOLOGY_INI` | `sbatch jobs/frontend.sh` (or `bash jobs/frontend.sh`, it is quick) |

```bash
sbatch ~/Hi-COLA/jobs/backend.sh                              # from anywhere
cd ~/Hi-COLA
sbatch jobs/backend.sh                                        # RUN_FILE from the script
sbatch -N 2 -n 128 --time=1:00:00 jobs/backend.sh             # other resources
sbatch jobs/backend.sh <run.lua>                              # a run file on the command line replaces RUN_FILE
sbatch --dependency=afterok:<backend job id> jobs/halo_finder.sh   # start once the backend is done
```

- The halo finder processes one redshift per job; edit `REDSHIFT` and submit again for the others.
- Logs: `HiCOLA/Output/logs/<job name>_<job id>.out` and `.err`, written to the full path in
  the `#SBATCH -o/-e` lines (`/home/biel/Hi-COLA/HiCOLA/Output/logs/`), whatever folder you
  submit from. Each log starts with what
  was run: code commit, conda env, run and data roots, and for a backend or halo-finder job
  the run file with its box, particles, steps, gravity, cosmology, initial P(k), seed and
  output redshifts (plus the FoF settings, halo cosmology and HMF mass grid for the halo finder).
- Everything else is set in the run file; the only other inputs are the data root and env
  name in `jobs/hicola_env.sh`. The halo finder takes the box size from the
  run file and the cosmology from the run's target.
- At the end of the backend and halo-finder jobs, `scripts/collect_products.py` copies the
  P(k) and HMF of each output redshift to `hicola_data/hicola/` (section 6). It can be rerun
  by hand: `python scripts/collect_products.py <target> <model> <realization>`.
- `lua jobs/run_info.lua <run.lua>` prints a run's name, output folder, box size and
  redshifts (with the `hicola` env active).

## 5. Halo finder and measurements

- **FoF** (`scripts/run_halo_catalog.py`, nbodykit, MPI): linking length, minimum members per
  halo and mass definition set in `jobs/halo_finder.sh` (so far 0.2; 30 for COLA_NM4096 and
  50 for Abacus_Small; 200c, the only definition the catalog format supports). The box size
  and the cosmology used for the masses are those of the run file.
- **HMF** (`scripts/run_halo_outputs.py`): cumulative n(>M) and dn/dlog10M with Poisson
  errors, on the mass grid set in `jobs/halo_finder.sh`: `reference` (the masses of the
  target's reference HMF at the same redshift, COLA_NM4096 only) or `MIN MAX N` log bins.
- The Python scripts have no hidden defaults: every setting is a required argument, passed by
  the job scripts.
- **P(k)**: Hi-COLA's own measurement on each snapshot (`pofk_z<z>.txt`, shot noise
  subtracted).
- **Abacus reference P(k)** (`scripts/run_abacus_power_spectrum.py`): subsamples, Nmesh and
  k-bins set in `jobs/power_spectrum_abacus.sh` (so far subsample A (3%), 1024, 60); the
  estimator is fixed in the script: halo + field particles, TSC, interlaced, window compensated,
  shot noise subtracted.
- `scripts/euclid_reference_cosmology.py` and `scripts/abacus_reference_cosmology.py` hold the
  reference cosmologies used by the HOD notebooks (they match the run files).
- `scripts/generate_abacus_cosm000_pk.py` regenerates the Abacus initial P(k) with CLASS and
  checks its σ8.

## 6. Reference data and kept measurements

`hicola_data/README.md` describes the data in detail. File names:

```
hicola_data/reference/<sim>/<model>/Box<n>/<sim>_<model>_Box<n>_z<z>_{PS,HMF}.txt
hicola_data/hicola/<target>/<model>/<realization>/HiCOLA_<target>_<model>_<realization>_z<z>_{PS,HMF}.txt
```

- `COLA_NM4096`: Bart's COLA boxes, `LCDM` and `nDGP_N1.2`, Box1–10, z = 1.65, 1.4, 1.2, 1.0,
  PS and HMF (spherical overdensity; check its definition in `papers/Bart.pdf` before
  comparing with M200c).
- `Abacus_Small`: AbacusSummit `small_c000_ph3000`, `LCDM`, Box1, z = 0.5, PS only.
- `raw/` folders hold the original data untouched.
- Redshifts are always written with 3 decimals (`z1.650`), as in Hi-COLA's own file names.

## 7. Notebooks

`Notebooks/`. Each notebook's first cell puts `scripts/` on the path and imports two modules:

- `hicola_paths as hp` gives every file location:

  ```python
  hp.reference_file('COLA_NM4096', 'LCDM', box, 1.65, 'PS')        # reference measurement
  hp.hicola_file('Abacus_Small', 'nDGP_N1.2', 's1234', 0.5, 'PS')   # kept Hi-COLA measurement
  hp.pofk_file('Abacus_Small', 'nDGP_N1.2', 's1234', 0.5)           # in the run folder on scratch
  hp.hmf_file(...), hp.halo_catalog_dir(...), hp.standalone_run_dir('cuGal')
  ```

- `analysis as an` loads and combines measurements:

  ```python
  k, P_mean, P_std = an.load_reference('COLA_NM4096', 'LCDM', 1.65, 'PS', boxes=range(1, 11))  # mean, std over boxes
  k, r_mean, r_std = an.reference_ratio('COLA_NM4096', 'nDGP_N1.2', 'LCDM', 1.65, 'PS')
  k, P = an.load_hicola('COLA_NM4096', 'LCDM', an.PAIRED_FIXED, 1.65, 'PS')       # mean over realizations
  an.ratio(x_num, y_num, x_den, y_den), an.interp_to(...), an.within(...)
  ```

  `load_hicola` reads the run folder on scratch, or the kept copy in `hicola_data/hicola/`
  once scratch has been cleaned up, so the comparisons keep working after the cleanup.

| notebook | what it does |
|---|---|
| `compare_COLA_NM4096.ipynb` | P(k) and HMF against Bart's boxes (mean ± 1σ), nDGP/LCDM ratios, (P − mean)/σ. Parameters: redshift, runs (e.g. `s1234` and/or the paired-fixed mean), reference boxes averaged |
| `compare_Abacus_Small.ipynb` | P(k) against the Abacus box (LCDM), nDGP/LCDM, HMFs. Parameters: redshift, nDGP strength, runs |
| `LRG_Abacus_Small.ipynb` | LRG mocks (HODDIES) on the paired-fixed Abacus_Small runs: P(k), bias, centrals only, RSD multipoles |
| `HOD.ipynb`, `HOD_Abacus.ipynb` | HODDIES mocks on one run; parameters in `HOD_parameter_file_z<z>.yaml` |
| `frontend.ipynb` | plots of the frontend tables in `HiCOLA/Output/` |
| `scan.ipynb` | scan of the nDGP frontend parameters (template `nDGP_scan/nDGP_scan.ini`) |
| `halos.ipynb`, `overdensity.ipynb` | older analyses kept as starting points; they still point to runs that were removed |

In the comparison notebooks, all figures follow the parameter cell at the top. A run that
doesn't exist yet is skipped with a message rather than stopping the notebook.

**Jupyter on a compute node** (heavy notebooks need more memory than the login node):

```bash
Sinteract -n 72 -c 1 -m 256G -t 04:00:00 -p academic
conda activate hicola
ipnport=$(shuf -i8000-9999 -n1)
jupyter notebook --no-browser --port=${ipnport} --ip=$(hostname -i)
```

## 8. Adding a reference simulation

1. Put its data in `hicola_data/reference/<sim>/<model>/Box<n>/` with the file names above
   (originals in `raw/`), and describe it in `hicola_data/README.md`.
2. Add the folder `Backend_parameter_files/<sim>/` with run files: copy existing ones and
   change box, cosmology, initial P(k) (in `Initial_power_spectra/`), output redshifts and the
   names (section 3).
3. For the HOD notebooks, add a `scripts/<name>_reference_cosmology.py` if the cosmology is new.

## 9. Setting up on another machine

1. Clone `BielThi/Hi-COLA` and `BielThi/HODDIES` into the same folder and copy
   `hicola_data/` next to them (or export `HICOLA_DATA_ROOT` / `HODDIES_ROOT`).
2. Create the `hicola` env (section 2), build `nbody` with it activated, run `nbstripout --install`.
3. Two machine-specific paths are written in full in the files:
   - the scratch path in the `output_folder` of the run files, and
   - the log folder in the `#SBATCH -o/-e` lines of the job scripts.

   If they differ, change them with
   ```bash
   sed -i 's#/scratch/biel/Hi-COLA_Output#<new scratch path>#' HiCOLA/Parameter_files/Backend_parameter_files/*/*.lua
   sed -i 's#/home/biel/Hi-COLA/#<new Hi-COLA path>/#' jobs/*.sh
   ```
   and export the scratch path as `HICOLA_RUN_ROOT` for the notebooks. Adapt the other
   `#SBATCH` lines in `jobs/*.sh` (partition, cores per node, `--exclusive`, memory) to the
   cluster. A second clone on jed needs the same `sed` for its log folder.
4. Check with a small run, e.g. a copy of a run file with `particle_Npart_1D = 192`,
   `timestep_nsteps = {10}` and the derived meshes overridden, submitted with
   `sbatch -N 1 -n 16 --time=0:30:00 jobs/backend.sh <that file>`.

## 10. Known issues

- Hi-COLA creates only the last level of `output_folder`; `jobs/backend.sh` creates the
  nested run folder before starting. Run `nbody` by hand only after a `mkdir -p` of it.
- `Frontend_parameter_files/nDGP/nDGP.ini` does not run (`Hubble0 = '1.0/Hubble_dS'` is
  not a number). The nDGP tables used by `standalone/nDGP_action.lua` exist only in
  `HiCOLA/Output/` on jed.
- Reference HMF of COLA_NM4096 is spherical-overdensity with an unstated Δ; Hi-COLA HMFs are
  FoF + M200c.
- Abacus_Small has no reference HMF yet.
