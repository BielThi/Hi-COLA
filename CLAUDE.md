# Hi-COLA comparison pipeline: notes for Claude

Read `PIPELINE.md` first: it describes the whole pipeline (layout, runs, jobs, data, notebooks,
setup on a new machine, troubleshooting, history). The data is described in
`../hicola_data/README.md`.

## How this repository is meant to stay

- `HiCOLA/` is the original Hi-COLA code with its original layout and output file names. Don't
  restructure it; the only local change there is `HiCOLA/FML/COLASolver/Makefile`.
- Run files (`HiCOLA/Parameter_files/Backend_parameter_files/<target>/<model>_<realization>.lua`)
  are complete, explicit Hi-COLA parameter files. No `dofile`, shared base files or generated
  values. A new run is a copied and edited file following the naming in PIPELINE.md section 5.
- Parameters the user chooses are visible where they are used: run files, the parameter block
  at the top of each `jobs/*.sh`, the parameter cell of each notebook. No hidden defaults in
  Python (settings are required arguments). Paths and other dependencies are automatic and
  printed in the job logs.
- Jobs are submitted directly: `sbatch jobs/<script>.sh`, no wrapper scripts.
- Prefer the simplest change that keeps these properties; ask before adding automation.

## Practical rules

- Never commit data or outputs: `hicola_data/` is outside the repo, `HiCOLA/Output/*` is
  ignored, notebooks are stripped by nbstripout (`nbstripout --install` once per clone).
- Scratch (`/scratch/$USER`) is cleaned up automatically: only `hicola_data/hicola/` keeps
  results. Don't overwrite reference files in `hicola_data/reference/`.
- Test changes to jobs or scripts with a small run, never by touching production outputs:
  a copied run file with `particle_Npart_1D = 192`, `timestep_nsteps = {10}`, a test
  `output_folder` (keeping the naming), `HICOLA_DATA_ROOT` exported to a test folder, and
  `sbatch -N 1 -n 16 --time=0:30:00 --mem-per-cpu=4G` (PIPELINE.md section 11).
- When changing run files, check that every parameter is unchanged except the intended ones
  (evaluate the old and new files with `lua` and diff the globals).
- On jed (EPFL SCITAS): partition `academic`, 72 cores per node, max 7000 MB per core.
- Commit with clear messages; push to GitHub only when asked.
- Related repository: `../HODDIES` (fork `BielThi/HODDIES`), used by the HOD notebooks.
