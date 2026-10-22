# Comparison pipeline

This fork of Hi-COLA runs simulations matched to reference simulations and compares their
matter power spectra, halo mass functions and HOD galaxies:

```
frontend (optional)      backend (Hi-COLA)          halo finder (nbodykit)       analysis (notebooks)
.ini → expansion/force   run .lua → snapshots,      FoF → M200c halo catalog →   P(k), HMF vs reference,
tables (HiCOLA/Output)   P(k) at each output z      HMF + HODDIES halo array     HOD / LRG mocks
```

Hi-COLA itself (`HiCOLA/`) keeps its original layout and file names. What was added around it:
`jobs/` (SLURM jobs), `scripts/` (Python), `Notebooks/`, and the run files in
`HiCOLA/Parameter_files/Backend_parameter_files/`. The data lives outside the repository, in
`hicola_data/` (reference data, kept results) and on scratch (simulation outputs).

**This branch, `hicola-dev-v2`**, is the pipeline adapted to HODDIES dev-v2 (fork branch
`hicola-dev-v2`). It lives next to the `main` pipeline (old HODDIES), which stays as it was:

```
~/Hi-COLA        main             ~/dev-v2/Hi-COLA     hicola-dev-v2   (git worktrees of the
~/HODDIES        main             ~/dev-v2/HODDIES     hicola-dev-v2    same two repositories)
~/hicola_data                     ~/dev-v2/hicola_data -> ~/hicola_data
```

Both share the run files, the simulations on scratch, the halo catalogs and `hicola_data`. Only
the HOD outputs are kept apart, because the mocks differ between HODDIES versions:
`halo_catalog/hod_dev-v2/` and `..._HOD-dev-v2_<name>_...` here, `halo_catalog/hod/` and
`..._HOD_<name>_...` on `main` (`HOD_LABEL` in `scripts/hicola_paths.py`). Each clone has its own
`nbody` build, frontend tables (`HiCOLA/Output/`) and job logs. Submit the jobs of this pipeline
from `~/dev-v2/Hi-COLA/jobs/`.

Contents: 1 quick start · 2 design rules · 3 where things live · 4 software setup · 5 runs and
their names · 6 running jobs · 7 halo finder and measurements · 8 data and file formats ·
9 notebooks · 10 adding runs or reference simulations · 11 another machine · 12 the jed cluster ·
13 troubleshooting · 14 known issues · 15 history

## 1. Quick start: one run, end to end

```bash
# 1. pick (or create, section 10) a run file, e.g.
#    HiCOLA/Parameter_files/Backend_parameter_files/Abacus_Small/nDGP_N1.2_s1234-fixed.lua
# 2. set it as RUN_FILE in jobs/backend.sh (or give it on the command line) and submit
sbatch ~/Hi-COLA/jobs/backend.sh                       # → job id, e.g. 41000000
# 3. set RUN_FILE, REDSHIFT, NMIN, HMF_MASS_GRID in jobs/halo_finder.sh and submit, one job per redshift
sbatch --dependency=afterok:41000000 ~/Hi-COLA/jobs/halo_finder.sh
# 4. follow the jobs and read the logs
squeue -u $USER
less ~/Hi-COLA/HiCOLA/Output/logs/backend_41000000.out
# 5. compare: open Notebooks/compare_Abacus_Small.ipynb (or compare_COLA_NM4096.ipynb), set its
#    parameter cell (redshift, runs, ...) and run it
```

The P(k) and HMF of every output redshift are copied to `hicola_data/hicola/` at the end of
the jobs, so the comparisons keep working after scratch has been cleaned up.

## 2. Design rules

These were chosen on purpose; keep them when changing the pipeline.

- **Hi-COLA stays as in the original code**: original folder layout, original output file
  names, run files that are complete Hi-COLA `.lua` parameter files (no `dofile`, no
  generated values). The only change inside `HiCOLA/` is the `nbody` Makefile (section 4).
- **Parameters are explicit; dependencies are automatic.** Everything that is a choice
  (simulation settings, FoF settings, HMF bins, estimator settings, which runs and boxes a
  plot uses) is written where it is used: in the run file, in the parameter block at the top
  of the job script, or in the parameter cell of the notebook. There are no hidden defaults:
  the Python scripts take every setting as a required argument. What follows from those
  choices (paths, the box size and cosmology of a run for the halo finder, which boxes exist)
  is worked out automatically and printed in the job logs.
- **Every job log starts with what was run**: code commit, conda env, paths, and all the
  settings of the run (section 6).
- **Code in git, data outside.** Reference data and kept results live in `hicola_data/` next
  to the repository, simulation outputs on scratch, frontend tables and logs in
  `HiCOLA/Output/` (git-ignored). Notebooks are stored without outputs (nbstripout).
- **Machine-specific values are written in full in a few places** (the scratch path in the
  run files' `output_folder`, the log folder in the `#SBATCH -o/-e` lines, the conda paths in
  `jobs/hicola_env.sh`) and changed with one `sed` on a new machine (section 11).

## 3. Where things live

```
~/                                  (any parent folder works, see section 11)
├── Hi-COLA/                        this repository                                 $HICOLA_ROOT
│   ├── CLAUDE.md                   pointer to this file for Claude Code
│   ├── PIPELINE.md                 this file
│   ├── HiCOLA/                     Hi-COLA code (frontend, FML backend), original layout
│   │   ├── FML/COLASolver/nbody    the backend binary (built, not in git)
│   │   ├── Output/                 frontend tables; logs/ = job logs (contents not in git)
│   │   └── Parameter_files/        frontend .ini files, backend .lua run files (section 5)
│   ├── jobs/                       SLURM job scripts (section 6) + hicola_env.sh, run_info.lua
│   ├── scripts/                    Python: paths, analysis helpers, halo finder, measurements
│   └── Notebooks/                  analysis notebooks and HOD parameter files (section 9)
├── HODDIES/                        HOD code (fork BielThi/HODDIES)                 $HODDIES_ROOT
├── hicola_data/                    reference data + kept measurements, README.md   $HICOLA_DATA_ROOT
└── dev-v2/                         this branch: Hi-COLA/, HODDIES/ (worktrees), hicola_data -> ../hicola_data
/scratch/biel/Hi-COLA_Output/       simulation outputs: <target>/<model>/<realization>/   run root
```

| what | where | in git? | survives? |
|---|---|---|---|
| code, run files, job scripts, notebooks (without outputs) | `Hi-COLA/` | yes, `BielThi/Hi-COLA` | yes |
| frontend tables, job logs | `Hi-COLA/HiCOLA/Output/` | no | yes (regenerable) |
| snapshots, P(k) per time step, halo catalogs | scratch run root, ~85 GB per 1024³ run | no | **no: scratch is cleaned up automatically after some weeks** |
| reference simulations, kept P(k)/HMF of every run | `hicola_data/` (6.4 GB, mostly raw Abacus catalogs) | no | yes (home) |

Environment variables (all optional; defaults in `jobs/hicola_env.sh` and
`scripts/hicola_paths.py`, which agree):

| variable | default | used for |
|---|---|---|
| `HICOLA_ROOT` | the repository the script is in | everything |
| `HICOLA_DATA_ROOT` | `<HICOLA_ROOT>/../hicola_data` | reference data, kept results |
| `HICOLA_RUN_ROOT` | `/scratch/$USER/Hi-COLA_Output` | where notebooks look for run folders; in the jobs it is replaced by the run root of the run file's `output_folder` |
| `HODDIES_ROOT` | `<HICOLA_ROOT>/../HODDIES` | HODDIES import in the HOD notebooks |
| `HICOLA_CONDA_ENV` | `hicola` | conda env activated by the jobs |

## 4. Software setup

**Conda env `hicola`** (Python 3.11; not described by a file in the repository, create it on
each machine). Packages used: `nbodykit` (FoF, P(k)), `mpi4py`, `abacusutils` (Abacus
catalogs), `cosmoprimo` + `pyclass` (CLASS), `numba`, `numpy`, `scipy`, `pandas`,
`matplotlib`, `astropy`, `fitsio`, `pynbody`, `h5py`, `hdf5plugin` (imported by HODDIES), `jax` +
`jax-power` (`pip install "jax-power @ git+https://github.com/adematti/jax-power.git"`, with
`lsstypes`; for the HOD P(k)), `sympy` + `configobj` (frontend),
`jupyter`, `nbstripout`, plus FFTW (with MPI and threads), GSL, Lua (the `lua` command is
used by the jobs) and an MPI compiler wrapper (`mpicxx`) to build the backend. The jobs
activate the env themselves.

**Backend binary** `HiCOLA/FML/COLASolver/nbody`: `cd HiCOLA/FML/COLASolver && make` with the
`hicola` env activated. The Makefile (the one file changed inside `HiCOLA/`) takes FFTW, GSL
and Lua from the active env (`$(CONDA_PREFIX)`) and FML from the clone it is in, so every
clone builds against its own code. Build options (MPI, OpenMP, FFTW threads, GSL, Lua) are
the `USE_*` lines at the top of the Makefile.

**HODDIES**: the fork `github.com/BielThi/HODDIES`, branch `hicola-dev-v2`, checked out next to
this repository (`~/dev-v2/HODDIES`, a worktree of `~/HODDIES`, whose `main` serves the `main`
pipeline). The HOD job and notebooks import it from there (`$HODDIES_ROOT`), not
from the `HODDIES` package installed in the env, which is a different copy. `hicola-dev-v2` is
the upstream branch `dev-v2` (remote `upstream` = `github.com/antoine-rocher/HODDIES`) plus two
fixes: the merge conflict left in the imports of `hod.py` resolved, and `torch` imported only
by the emulator functions that use it (the env has no torch). Since dev-v2 a HOD object is
built from a halo catalogue:
```python
from HODDIES.sim_loader import Base_catalogue
base_catalog = Base_catalogue({name: arr[name] for name in arr.dtype.names}, boxsize=boxsize, z_simu=redshift)
HOD_obj = HOD(base_catalog, param_file=PARAM_FILE)   # the HOD job also passes nthreads and seed
HOD_obj.cosmo = cosmology                           # for the RSD of get_wp / get_xiells
```
(the halo array is passed as a dict: a structured array is refused). The mock of a seed depends
on the HODDIES version: mocks made with the fork's old `main` (commit `2572a69`, numba
parallelism off) differ from dev-v2 ones and must not be mixed; the HOD job records the HODDIES
branch and commit in each output's `.json`.

**hi_class**: the env on jed has a `classy` built from `~/hi_class_public`
(`miguelzuma/hi_class_public`, commit `50f447c4`). Nothing in the pipeline uses it any more
(CLASS is used through cosmoprimo/`pyclass`); not needed on a new machine.

**nbstripout** removes notebook outputs from what git stores (the files on disk keep them).
It is declared in `.gitattributes`; each new clone runs `nbstripout --install` once (with the
env active).

**Git**: `origin` = `git@github.com:BielThi/Hi-COLA.git`, branch `main`. Tag `upstream-base`
marks the original Hi-COLA (`Hi-COLA_public`) commit this fork started from; tag
`first_release` comes from the original repository.

## 5. Runs and their names

A run is matched to a reference simulation (the **target**: same box, cosmology and initial
P(k)) and named by a **model** and a **realization**:

| | values |
|---|---|
| target | `COLA_NM4096`: 1024 Mpc/h, 1024³ particles, Euclid cosmology, Bart's initial P(k), z_init 19, outputs z = 1.65, 1.4, 1.2, 1.0 · `Abacus_Small`: 500 Mpc/h, 1024³ particles, AbacusSummit c000 cosmology, its initial P(k), z_init 49, outputs z = 1.1, 0.8, 0.5. Both 30 time steps |
| model | `LCDM` (`gravity_model = "GR"`) · `nDGP_N<rcH0>`, e.g. `nDGP_N1.2` = `gravity_model = "DGP"` with `gravity_model_dgp_rcH0overc = 1.2` (Hi-COLA's built-in DGP) |
| realization | `s<seed>` · `s<seed>-fixed` (`ic_fix_amplitude = true`) · `s<seed>-fixed-inv` (also `ic_reverse_phases = true`); the pair -fixed / -fixed-inv is averaged as "paired-fixed" |

```
HiCOLA/Parameter_files/Backend_parameter_files/
├── Hi-COLA_params.lua                          original Hi-COLA template
├── LCDM_1306dot3219_matterpower_z0.000.dat     original Hi-COLA example P(k) (σ8 ≈ 0.92, testing only)
├── Initial_power_spectra/
│   ├── Bart_InitPk.txt                         COLA_NM4096 (copy of hicola_data/reference/COLA_NM4096/raw/LCDM/AllBoxes/InitPk.txt)
│   └── abacus_cosm000_pk_cb_z0.000.dat         Abacus_Small (P_cb at z=0, from scripts/generate_abacus_cosm000_pk.py)
├── COLA_NM4096/   LCDM_s1234, LCDM_s1234-fixed, LCDM_s1234-fixed-inv,
│                  nDGP_N1.2_s1234, nDGP_N1.2_s1234-fixed, nDGP_N1.2_s1234-fixed-inv   (.lua)
├── Abacus_Small/  the same six + nDGP_N1.0_s1234
└── standalone/    runs not matched to a reference, kept as examples:
    ├── cuGal.lua                               cubic Galileon through the frontend (gravity_model = "HiCOLA")
    └── nDGP_action.lua                         nDGP from its Horndeski action through the frontend
```

Every run file is a complete Hi-COLA parameter file, like the original template: every
setting is written in it. Within a target, the files differ only in `simulation_name`,
`output_folder`, `gravity_model` (and `gravity_model_dgp_rcH0overc`), `ic_random_seed`,
`ic_fix_amplitude` and `ic_reverse_phases`. Input paths (initial P(k), frontend tables) are
relative to the Hi-COLA folder, where the jobs run `nbody`. Some values are Lua expressions
of others, as in the template (e.g. `force_nmesh = 3 * particle_Npart_1D`), so changing the
particle number also changes those.

The pipeline finds a run's outputs by its name, so a run file must follow this naming;
`jobs/run_info.lua` checks it and the jobs stop with a message otherwise:

| in the run file | must be |
|---|---|
| location | `<target>/<model>_<realization>.lua` (realization starts with `s<digits>`) |
| `simulation_name` | `"<target>_<model>_<realization>"` |
| `output_folder` | `"<run root>/<target>/<model>/<realization>/"`, run root `/scratch/biel/Hi-COLA_Output` on jed |

Standalone runs: `standalone/<name>.lua`, `output_folder = "<run root>/standalone/<name>/"`,
own `simulation_name`. They need frontend tables in `HiCOLA/Output/` (not in git).

Hi-COLA writes its usual files to the output folder:

```
pofk_<simulation_name>_cb_z<z>.txt         P_cb(k) and linear P(k) at every time step
gravitymodel_<simulation_name>_k*.txt, cosmology_<simulation_name>.txt
snapshot_<simulation_name>_z<z>/           one per output redshift
├── gadget_z<z>.*                          particles (Gadget format, one file per MPI task)
├── pofk_z<z>.txt                          P(k) of the snapshot: this is the "PS" that is kept
└── halo_catalog/                          written by the halo finder
    ├── halos200c_catalog.npy (+ .json)    FoF halos with M200c (section 8)
    ├── halos200c_hmf.txt (+ .json)        HMF: this is the "HMF" that is kept
    └── halos200c_hoddies.npy (+ .json)    the same halos in HODDIES' format
```

## 6. Running jobs

Each script in `jobs/` has a **user parameters** block at the top: edit it, then submit the
script with `sbatch`, from any folder. Paths in the parameter blocks are relative to the
Hi-COLA folder; a file given on the command line is relative to where you submit. The
`#SBATCH` lines are the production settings on jed; edit them in the script or override them
on the command line.

| script | parameters | resources (default) | run with |
|---|---|---|---|
| `backend.sh` | `RUN_FILE` | 8 nodes, 512 tasks, 6 GB/task, 5 h, exclusive | `sbatch jobs/backend.sh` |
| `halo_finder.sh` | `RUN_FILE`, `REDSHIFT`, `LINKING_LENGTH`, `NMIN`, `MDEF`, `HMF_MASS_GRID` | 8 nodes, 512 tasks, 6 GB/task, 4.5 h | `sbatch jobs/halo_finder.sh` |
| `hod_mock.sh` | `RUN_FILE`, `REDSHIFT`, `HOD_PARAM_FILE`, `TRACER`, `SEEDS`, `NTHREADS`, `LOS`, `MESHSIZE`, `ELLS`, `K_STEP`, `MARK_RCUT`, `MARK_COEFFICIENTS`, `NAME` | 1 task, 32 threads, 6 GB/thread, 1 h | `sbatch jobs/hod_mock.sh` |
| `hod_fit.sh` | `RUN_FILE`, `REDSHIFT`, `HOD_PARAM_FILE`, `TRACER`, `PRIORS`, `DATA_RPPI`, `DATA_SMU`, `STATISTICS`, `RP_REBIN`, `RP_RANGE`, `PIMAX`, `ELLS`, `S_REBIN`, `S_RANGE`, `CORRELATION`, `LOS`, `MOCK_SEED`, `MODEL_SEEDS`, `METHOD`, `MAXITER`, `POPSIZE`, `XTOL`, `START`, `START_SEED` | 1 node, 18 MPI ranks x 4 threads, 3 GB/thread, 4 h | `sbatch jobs/hod_fit.sh` |
| `power_spectrum_abacus.sh` | `REDSHIFT`, `MODEL`, `LOAD_A`, `LOAD_B`, `NMESH`, `KBINS` | 1 task, 16 threads, 40 GB, 1.5 h | `sbatch jobs/power_spectrum_abacus.sh` |
| `frontend.sh` | `MODEL_INI`, `COSMOLOGY_INI` | 1 task, 6 GB, 1 h | `sbatch jobs/frontend.sh` or `bash jobs/frontend.sh` |

```bash
sbatch ~/Hi-COLA/jobs/backend.sh                              # from anywhere
sbatch jobs/backend.sh <run.lua>                              # a run file on the command line replaces RUN_FILE
sbatch -N 2 -n 128 --time=1:00:00 jobs/backend.sh             # other resources
sbatch --job-name=LCDM_s1234 jobs/backend.sh                  # other job name (= log file name)
sbatch --dependency=afterok:<backend job id> jobs/halo_finder.sh   # start once the backend is done
```

Settings used so far, per target (set them in `halo_finder.sh` for each submission):

| | COLA_NM4096 | Abacus_Small |
|---|---|---|
| `REDSHIFT` | 1.650, 1.400, 1.200, 1.000 | 1.100, 0.800, 0.500 (reference at 0.500 only) |
| `LINKING_LENGTH` | 0.2 | 0.2 |
| `NMIN` | 30 | 50 |
| `MDEF` | 200c | 200c |
| `HMF_MASS_GRID` | `reference` (Bart's bins) | `1e11 1e16 60` |

- The halo finder processes one redshift per job; edit `REDSHIFT` and submit again for the others.
- **Logs**: `HiCOLA/Output/logs/<job name>_<job id>.out` and `.err`, at the full path written
  in the `#SBATCH -o/-e` lines (`/home/biel/dev-v2/Hi-COLA/HiCOLA/Output/logs/` on this branch). Each log starts
  with what was run: job, code commit (and whether there were uncommitted changes), conda
  env, run and data roots; for the backend and halo finder the run file with its box,
  particles, steps, z_init, gravity, cosmology, initial P(k), seed, fixed/inverted, output
  redshifts and output folder; for the halo finder also the FoF settings and HMF mass grid;
  for the Abacus P(k) its settings, input and output paths; for the frontend the `.ini` files.
- **Kept results**: at the end of the backend, halo-finder and HOD-mock jobs,
  `scripts/collect_products.py` copies the P(k), HMF and HOD galaxy P(k) and marked P(k) of each output redshift to
  `hicola_data/hicola/` (section 8). It can be rerun by hand:
  `python scripts/collect_products.py <target> <model> <realization>`.
- `lua jobs/run_info.lua <run.lua>` (env active) prints a run's identity and settings, or the
  naming error. Helpers in `jobs/hicola_env.sh`: `hicola_load_run`, `hicola_print_env`,
  `hicola_print_run`.
- Outside the job scripts and run files, the only inputs are the paths and env name in
  `jobs/hicola_env.sh`, and those are printed in every log.

## 7. Halo finder and measurements

- **FoF** (`scripts/run_halo_catalog.py`, nbodykit, MPI via `srun --mpi=pmix`): linking
  length, minimum members and mass definition from `halo_finder.sh`. The box size and the
  cosmology used for the masses (h, Ω_cdm, Ω_b, Ω_ν, Ω_Λ, T_CMB, N_eff, n_s) are read from the
  run file. Only `MDEF=200c` is accepted: the catalog fields and files are named for it.
- **HMF** (`scripts/run_halo_outputs.py`): cumulative n(>M) and dn/dlog10M with Poisson
  errors. `HMF_MASS_GRID` is `reference` (the masses of the target's Box1 LCDM reference HMF at
  the same redshift, so both compare bin by bin; stops with an error if there is none) or
  `MIN MAX N` (N log-spaced masses in Msun/h). It also writes the HODDIES halo array.
- **P(k) of the runs**: Hi-COLA's own measurement on each snapshot, `pofk_z<z>.txt`, shot
  noise subtracted, settings from the run file (`pofk_*` lines).
- **HOD galaxy P(k)** (`jobs/hod_mock.sh`: `scripts/run_hod_mock.py`,
  `scripts/galaxy_power_spectrum.py`): a HODDIES mock of the chosen tracer is made from the halo
  array, shifted to redshift space along `LOS` by v/(aH(z)) with H(z) of the run's cosmology, and
  measured with the estimator of the EMC/acm pipeline (copied from acm's
  `scripts/emc/measurements/abacus_small.py`): jaxpower, TSC, interlacing 3, compensation, shot
  noise from jaxpower; `MESHSIZE`, `ELLS`, `K_STEP` from the job (EMC values: 512, 0 2 4, 0.001),
  so it compares directly with the EMC measurements. If acm changes its estimator, the copy
  does not follow. Runs on CPU (jax).
  The same mock's **marked P(k)** (`compute_marked_spectrum`, copied from acm's
  `MarkedPowerSpectrumMultipoles`, box without randoms; not yet released in acm) uses the same
  positions and mesh settings, with the mark m = Σₙ cₙ δ_Rⁿ (δ_R: galaxy density contrast smoothed
  by a Gaussian of radius `MARK_RCUT`, cₙ = `MARK_COEFFICIENTS`), normalised by its galaxy-weighted
  mean m̄; box normalization and shot noise of the marked galaxies, as in acm. The job refuses to
  run until both mark parameters are set. With the mark m = 1 it gives the P(k) exactly.
  One job makes one mock, P(k) and marked P(k) per seed of `SEEDS` (the halo array is loaded
  once, ~1.5 min per seed at mesh 512); `NAME` contains `{seed}` (default `<TRACER>_seed{seed}`), so each seed has
  its own file. A seed gives the same mock alone or in a list, but the mock also depends on
  `NTHREADS`: keep it fixed across seeds and runs. Use the same seeds for the models compared.
- **Abacus reference P(k)** (`scripts/run_abacus_power_spectrum.py`): halo + field particles
  of the chosen subsamples (A = 3%, B = 7%), TSC, interlaced, window compensated, shot noise
  subtracted (these estimator choices are fixed in the script); Nmesh and k-bins from the job
  (so far A only, 1024, 60).
- **Abacus initial P(k)**: `scripts/generate_abacus_cosm000_pk.py` computes P_cb(k, z=0) for
  AbacusSummit c000 with CLASS (cosmoprimo) and checks σ8_cb = 0.8114.
- `scripts/euclid_reference_cosmology.py` and `scripts/abacus_reference_cosmology.py` hold the
  reference cosmologies used by the HOD notebooks; they match the run files.
- All Python scripts take every setting as a required argument (`--help` lists them); only
  paths, the Gadget particle type and thread counts have defaults.

## 8. Data and file formats

`hicola_data/README.md` describes the data in detail (sources, raw data, how each file was made).

```
hicola_data/reference/<sim>/<model>/Box<n>/<sim>_<model>_Box<n>_z<z>_{PS,HMF}.txt   (+ raw/ originals)
hicola_data/hicola/<target>/<model>/<realization>/HiCOLA_<target>_<model>_<realization>_z<z>_{PS,HMF}.txt
```

- Reference data: `COLA_NM4096` (Bart's COLA boxes; `LCDM` and `nDGP_N1.2`, Box1–10,
  z = 1.65, 1.4, 1.2, 1.0, PS and HMF; spherical-overdensity HMF, see `papers/Bart.pdf`) and
  `Abacus_Small` (AbacusSummit `small_c000_ph3000`, `LCDM`, Box1, z = 0.5, PS only).
- Redshifts are always written with 3 decimals (`z1.650`), as in Hi-COLA's own file names.
- **PS files**: text, `#` header lines, columns k [h/Mpc], P(k) [(Mpc/h)³], then optional
  columns (P_linear for Hi-COLA, N_modes for Abacus).
- **HMF files**: text, `#` header, columns M [Msun/h], n(>M) [h³/Mpc³]; Hi-COLA HMFs add
  σ_n(>M), dn/dlog10M, σ, counts(>M), counts per bin. A `.json` next to a file holds its
  metadata (run, settings, cosmology, mass grid source).
- **Halo catalog** `halos200c_catalog.npy`: structured array with `x y z` [Mpc/h], `vx vy vz`
  and `vrms` [peculiar km/s], `m200c`, `log10_m200c` [Msun/h], `r200c`, `rs` [kpc/h], `c`, `length`
  (FoF particles), `halo_id`; `.json` with box size, redshift, FoF settings, particle mass,
  cosmology and `velocity_definition: "peculiar"`. Hi-COLA's Gadget snapshots store v_pec/√a;
  the halo finder multiplies by √a (checking a against the snapshot header). Catalogs without
  `velocity_definition` were made before this fix (2026-10-14) and have velocities too large by
  √(1+z): the halo outputs, the HOD-mock job and the HOD/LRG notebooks refuse them; rerun the
  halo finder.
- **HODDIES array** `halos200c_hoddies.npy`: `x y z vx vy vz Mh log10_Mh Rh Rs c Vrms halo_id`
  (peculiar velocities).
- **HOD galaxy P(k)** `halo_catalog/hod_dev-v2/<name>_PS.txt` (+ `.h5`, `.txt.json`; `jobs/hod_mock.sh`,
  `scripts/run_hod_mock.py`), kept as `HiCOLA_<simulation_name>_z<z>_HOD-dev-v2_<name>_PS.txt`: the
  redshift-space power spectrum multipoles of a HODDIES mock (one file per seed), columns k [h/Mpc], N_modes,
  P_0, P_2, P_4 [(Mpc/h)³, shot noise subtracted]. The `.h5` is jaxpower's own format (`lsstypes`),
  the same as the EMC measurements; the `.json` holds the run, HOD settings and parameters,
  number of galaxies and density, H(z) and the estimator settings. The mock itself is not kept.
  The k-bins are 0.001 h/Mpc wide (as in EMC): the first row is the k = 0 mode and bins without
  modes (N_modes = 0) are NaN; skip both, and rebin (weighting by N_modes) for plots.
- **HOD marked P(k)** `halo_catalog/hod_dev-v2/<name>_MPS.txt` (+ `.h5`, `.txt.json`), kept as
  `HiCOLA_<simulation_name>_z<z>_HOD-dev-v2_<name>_MPS.txt`: the marked power spectrum of the same mock,
  same columns and bins. The `.json` also holds the mark (R, coefficients, m̄, min and max of m on
  the mesh), and the `.h5` has m_bar, m_min, m_max as attributes, as acm writes them. The file
  name does not contain the mark: a rerun with another mark overwrites it (check the `.json`).
- **HOD fits** `hicola_data/hicola/<target>/<model>/<realization>/HiCOLA_<simulation_name>_z<z>_HOD-dev-v2_fit_<name>/`
  (`jobs/hod_fit.sh`, `scripts/run_hod_fit.py`, `hp.hod_fit_dir`; shown by `HOD_fit_results.ipynb`,
  read with `an.load_hod_fit`), written there directly (not on scratch). The name says what is
  fitted: HODDIES' fit model name and the statistics, `<tracer>_<HOD model>_<n free>p_<statistics>`,
  e.g. `LRG_SHOD_6p_wp+xi02`; a new fit of the same thing replaces the folder (the settings of each
  fit are in its `settings.json`):
  `settings.json` (all settings of the job, run, H(z), code versions, job log), a copy of the HOD
  parameter file, `data.npz` (the fitted data: separations, values, errors, data vector, covariance,
  bin edges), `evaluations/rank_<k>.txt` (every evaluation of rank k: parameters, chi2, galaxies,
  failed flag, time; written as the fit runs), `iterations.txt` (best point of each iteration),
  and at the end `best_fit.json` (best evaluation of the whole fit, chi2 and per dof, p-value,
  galaxies, minimizer status) and `best_fit_model.npz` (model vectors at the best fit for the fit
  seed and the extra seeds, HOD occupation). `HOD_fit.ipynb` (exploration) also writes its
  `fit_info.json` and HODDIES' `best_fit_result_*.npy` into a folder of this form.

## 9. Notebooks

`Notebooks/`. Each notebook's first cell finds the repository, puts `scripts/` on the path and
imports:

- `hicola_paths as hp`, every file location:

  ```python
  hp.reference_file('COLA_NM4096', 'LCDM', box, 1.65, 'PS')         # reference measurement
  hp.hicola_file('Abacus_Small', 'nDGP_N1.2', 's1234', 0.5, 'PS')    # kept Hi-COLA measurement
  hp.pofk_file('Abacus_Small', 'nDGP_N1.2', 's1234', 0.5)            # in the run folder on scratch
  hp.hmf_file(...), hp.halo_catalog_dir(...), hp.run_file(...), hp.standalone_run_dir('cuGal')
  ```

- `analysis as an`, loading and combining measurements:

  ```python
  k, P_mean, P_std = an.load_reference('COLA_NM4096', 'LCDM', 1.65, 'PS', boxes=range(1, 11))
  k, r_mean, r_std = an.reference_ratio('COLA_NM4096', 'nDGP_N1.2', 'LCDM', 1.65, 'PS', boxes=range(1, 11))
  k, P = an.load_hicola('COLA_NM4096', 'LCDM', an.PAIRED_FIXED, 1.65, 'PS')   # mean over realizations
  k, P = an.load_hod_spectrum('Abacus_Small', 'LCDM', 's1234', 0.5, 'LRG_seed42', dk=0.01)  # HOD P(k), {ell: P_ell}
  k, P = an.load_hod_spectra('Abacus_Small', 'LCDM', 's1234', 0.5, 'LRG_seed{seed}', seeds=range(1, 11), dk=0.01)  # {ell: (n_seeds, n_k)}
  k, M = an.load_hod_spectra(..., statistic='MPS')                                     # marked P(k), same for load_hod_spectrum
  an.ratio(x_num, y_num, x_den, y_den), an.interp_to(...), an.within(...)
  ```

  `load_hicola` reads the run folder on scratch, or the kept copy in `hicola_data/hicola/` once
  scratch has been cleaned up. `boxes=None` means every `Box<n>` folder that exists.

| notebook | what it does | parameters (top cell) |
|---|---|---|
| `compare_COLA_NM4096.ipynb` | P(k) and HMF against Bart's boxes (mean ± 1σ), nDGP/LCDM ratios, (P − mean)/σ | redshift, models, runs (label → realizations averaged), reference boxes, M_min, y-range |
| `compare_Abacus_Small.ipynb` | P(k) against the Abacus box (LCDM), nDGP/LCDM, HMFs | redshift, models (e.g. nDGP_N1.2 or N1.0), runs, reference boxes, k_max, M_min, y-range |
| `LRG_Abacus_Small.ipynb` | LRG mocks (HODDIES) on the paired-fixed Abacus_Small runs: P(k), bias b²(k), centrals only, RSD multipoles | in its cells |
| `marked.ipynb` | S_X(k) = ⟨X_nDGP⟩/⟨X_GR⟩ − 1 for X = P₀, P₂ of the HOD galaxy P(k), then for the marked M₀, M₂, averaged over seeds (`hod_mock.sh`; checks all seeds and models use the same mark) | target, redshift, realization, HOD name (with `{seed}`), seeds, models, ells, dk (rebinning), k_max, show_seeds; curve = ratio of the seed means, band = std/√N of the per-seed S_X |
| `HOD.ipynb`, `HOD_Abacus.ipynb` | HODDIES mocks on one run; box size read from the halo catalog | target, model, realization, redshift; HOD in `HOD_parameter_file_z<z>.yaml` |
| `HOD_fit.ipynb` | exploration of HOD fitting with HODDIES (the fits themselves: `jobs/hod_fit.sh`) | target, model, realization, redshift, parameter file, tracer, nthreads, seed |
| `HOD_fit_results.ipynb` | results of a fit of `jobs/hod_fit.sh` (also while it runs): settings, best-fit parameters against priors, chi2 (total, per block, other seeds), data vs best-fit model with residuals, convergence, parameters of every evaluation, HOD occupation | target, model, realization, redshift, fit name (printed by the job) |
| `frontend.ipynb` | plots of the frontend tables in `HiCOLA/Output/` | table names |
| `scan.ipynb` | scan of the nDGP frontend parameters | template `Frontend_parameter_files/nDGP_scan/nDGP_scan.ini` |
| `halos.ipynb`, `overdensity.ipynb` | older analyses kept as starting points; they point to runs that no longer exist | — |

The HOD notebooks need the halo catalogs, which only exist on scratch (they are not kept).
In the comparison notebooks, a run that doesn't exist is skipped with a message.

**Jupyter on a compute node** (heavy notebooks need more memory than the login node):

```bash
Sinteract -n 72 -c 1 -m 256G -t 04:00:00 -p academic
cd ~/dev-v2/Hi-COLA        # the clone of the pipeline used (~/Hi-COLA for main)
conda activate hicola
ipnport=$(shuf -i8000-9999 -n1)
jupyter notebook --no-browser --port=${ipnport} --ip=$(hostname -i)
```

Start the server inside the clone: kernels of a remote server (e.g. VSCode connected to it)
start in the server's folder, and the notebooks find the repository (and through it the
HODDIES clone) from their working directory. From `~` the first cell fails with
`StopIteration`; from the other clone it silently uses that pipeline. One server per pipeline.

## 10. Adding runs or reference simulations

**A run**: copy a run file of the same target, rename it `<model>_<realization>.lua`, change
`simulation_name`, `output_folder` and the values that differ (gravity, seed, ...). Check it
with `lua jobs/run_info.lua <file>`. A new `nDGP_N<x>` model needs
`gravity_model = "DGP"` and `gravity_model_dgp_rcH0overc = <x>`.

**A reference simulation**:
1. Put its data in `hicola_data/reference/<sim>/<model>/Box<n>/` with the file names of
   section 8 (originals in `raw/`), and describe it in `hicola_data/README.md`.
2. Add `Backend_parameter_files/<sim>/` with run files: copy existing ones and change box,
   cosmology, initial P(k) (put it in `Initial_power_spectra/`), output redshifts, z_init and
   the names.
3. For the HOD notebooks, add a `scripts/<name>_reference_cosmology.py` if the cosmology is new.
4. Make a comparison notebook from one of the existing ones.

## 11. Setting up on another machine

1. Clone `BielThi/Hi-COLA` and `BielThi/HODDIES`, both on branch `hicola-dev-v2`, into the same folder, copy `hicola_data/`
   next to them (or export `HICOLA_DATA_ROOT` / `HODDIES_ROOT`).
2. Create the `hicola` conda env (section 4), build `nbody` with it activated, run
   `nbstripout --install` in the clone.
3. Machine-specific values written in full:
   - the scratch run root in the `output_folder` of the run files,
   - the log folder in the `#SBATCH -o/-e` lines of the job scripts,
   - in `jobs/hicola_env.sh`, the default run root (`/scratch/$USER/...`) and where conda is
     looked for (`$HOME/miniconda3`, `$HOME/anaconda3`, or a `conda` already on the PATH).

   ```bash
   sed -i 's#/scratch/biel/Hi-COLA_Output#<new scratch path>#' HiCOLA/Parameter_files/Backend_parameter_files/*/*.lua
   sed -i 's#/home/biel/dev-v2/Hi-COLA/#<new Hi-COLA path>/#' jobs/*.sh
   ```
   Export the scratch path as `HICOLA_RUN_ROOT` for the notebooks if it isn't
   `/scratch/$USER/Hi-COLA_Output`. A second clone on jed needs the log-folder `sed`.
4. Adapt the other `#SBATCH` lines (partition, cores per node, `--exclusive`, memory per
   core, time) to the cluster (section 12), and the launcher (`mpirun` in `backend.sh`,
   `srun --mpi=pmix` in `halo_finder.sh`) if needed.

   The repository keeps **jed's** `#SBATCH` lines (and log paths). On another cluster (e.g.
   Perlmutter: `--account`, `--qos`, `--constraint` instead of `--partition`, logs under its
   home), keep these edits local and don't commit them: before a `git pull`, set them aside
   with `git stash`, then `git stash pop` and resolve if a pulled change touches the same
   lines. New job scripts are committed with jed's values.
5. Test with a small run: copy a run file with `particle_Npart_1D = 192`,
   `timestep_nsteps = {10}` and a test `output_folder`/`simulation_name` (keeping the naming),
   then `sbatch -N 1 -n 16 --time=0:30:00 --mem-per-cpu=4G jobs/backend.sh <that file>` and
   the halo finder on it; to keep the test out of `hicola_data`, export `HICOLA_DATA_ROOT` to
   a test folder before submitting (the job inherits it).

## 12. The jed cluster (EPFL SCITAS)

- Partition `academic`: 72 cores and ~500 GB per node, at most 7000 MB per core
  (`--mem-per-cpu` above that is refused at submission). The production jobs use 64 tasks per
  node (8 × 64 = 512).
- Scratch `/scratch/$USER` is not backed up and old files are deleted automatically.
- SLURM creates a missing log folder; it runs a copy of the job script, which is why the jobs
  ask SLURM (`scontrol show job`) for the original script path to find the repository.
- `Sinteract` opens an interactive allocation (section 9). Jobs submitted from inside one
  inherit `SLURM_MEM_*` variables; the job scripts unset them.

## 13. Troubleshooting

| symptom | cause / fix |
|---|---|
| no log in `HiCOLA/Output/logs/` | check `sacct -j <id> -o JobID,State,ExitCode,WorkDir` and `scontrol show job <id>` (`StdOut=`); the log goes where the `#SBATCH -o` line says |
| job stops with "simulation_name is …, expected …" or "output_folder … should end with …" | the run file doesn't follow the naming (section 5); usually a copied file whose names weren't updated |
| `Failed to create snapshot folder` | `nbody` run by hand without creating the nested `output_folder` first (`mkdir -p`); the job script does it |
| halo-finder job `DependencyNeverSatisfied` | the backend job it waits for failed; read the backend log |
| "halo velocities are not marked as peculiar velocities" | the halo catalog was made before the velocity fix (2026-10-14); rerun `jobs/halo_finder.sh` on that run and redshift |
| `--mass-grid reference: there is no reference HMF` | `HMF_MASS_GRID="reference"` for a target/redshift without reference HMF; use `MIN MAX N` |
| sbatch: `exceeds MaxMemPerCPU` | more than 7000 MB per core requested; lower `--mem-per-cpu`/`--mem` or add cores |
| notebook: "Skipping … No PS for …" | the run (or that redshift) doesn't exist on scratch nor in `hicola_data/hicola/` |
| HOD notebook: missing `halos200c_hoddies.npy` | run the halo finder on that run and redshift (the run must still be on scratch) |
| `CondaError: Run 'conda init'` | conda not found as a shell function; `jobs/hicola_env.sh` looks in `$CONDA_EXE`, the PATH, `~/miniconda3`, `~/anaconda3` |

## 14. Known issues

- Hi-COLA creates only the last level of `output_folder` (FML `create_folder`);
  `jobs/backend.sh` creates the nested run folder before starting `nbody`.
- `Frontend_parameter_files/nDGP/nDGP.ini` does not run (`Hubble0 = '1.0/Hubble_dS'` is not a
  number). The nDGP tables used by `standalone/nDGP_action.lua` exist only in
  `HiCOLA/Output/` on jed.
- The COLA_NM4096 reference HMF is spherical-overdensity with an unstated Δ; Hi-COLA HMFs are
  FoF + M200c.
- Abacus_Small has no reference HMF yet.
- `halos.ipynb` and `overdensity.ipynb` point to runs that were removed.

## 15. History

Reorganised in October 2026 (commits `d3fb7a1` … `99bd72a`); `57ee534` is the last commit
with the old layout, and older run files can be recovered from it (`git show 57ee534:<path>`).

Old names → current names, for older notes, plots and logs:

| old | now |
|---|---|
| `LCDM_Bart`, `LCDM_Bart_a`, `LCDM_Bart_ma` | `COLA_NM4096/LCDM_s1234`, `_s1234-fixed`, `_s1234-fixed-inv` |
| `nDGP_Bart`, `nDGP_Bart_a`, `nDGP_Bart_ma` | `COLA_NM4096/nDGP_N1.2_s1234`, `-fixed`, `-fixed-inv` |
| `LCDM_Abacus`, `_a`, `_ma` / `nDGP_Abacus`, `_a`, `_ma` | `Abacus_Small/LCDM_…` / `Abacus_Small/nDGP_N1.2_…` |
| `nDGP1dot0_Abacus` | `Abacus_Small/nDGP_N1.0_s1234` |
| `/scratch/biel/<old name>/` | `/scratch/biel/Hi-COLA_Output/<target>/<model>/<realization>/` |
| `Comparison/Box<n>/Bart_{LCDM,nDGP}_{snap,hmf}_<s>.txt` | `hicola_data/reference/COLA_NM4096/{LCDM,nDGP_N1.2}/Box<n>/…_z<z>_{PS,HMF}.txt` (snap 0–3 = z 1.65, 1.4, 1.2, 1.0) |
| `Comparison/Abacus1/` | `hicola_data/reference/Abacus_Small/LCDM/Box1/` |
| `Hi-COLA_Output/Frontend/` | `HiCOLA/Output/` |
| `backend*.sh`, `halo_finder512*.sh`, `power_spectrum_abacus.sh` at the root | `jobs/` |
| `power_spectrum*.ipynb` (8 versions) | `compare_COLA_NM4096.ipynb`, `compare_Abacus_Small.ipynb`, `LRG_Abacus_Small.ipynb` |

Velocity fix (2026-10-14): halo catalogs made before it carry Gadget velocities (v_pec/√a)
labelled as km/s, too large by √(1+z) (1.22 at z = 0.5). Real-space results (positions, masses,
HMF, P(k)) are unaffected; anything in redshift space made from those catalogs (HOD mocks, RSD
multipoles, HODDIES clustering with RSD) should be redone after rerunning the halo finder.

Fixed during the reorganisation: the `_a`/`_ma` Bart runs (`LCDM_Bart_a/_ma`) wrote to the
same output folder as `LCDM_Bart`, and the nDGP `_a/_ma` Bart runs used seed 4321 and
z_init 49 instead of 1234 and 19 like their LCDM partners; runs made before should be redone.
