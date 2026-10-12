# Environment shared by the Hi-COLA pipeline jobs (sourced by jobs/*.sh).
# Every variable can be set beforehand to override its default, e.g. on another machine:
#   export HICOLA_RUN_ROOT=/scratch/other/Hi-COLA_Output

HICOLA_ROOT="${HICOLA_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
export HICOLA_ROOT
# Simulation outputs (large, regenerable): snapshots, P(k), halo catalogs. Used by the notebooks;
# in the jobs it is replaced by the run root of the run file's output_folder (hicola_load_run).
export HICOLA_RUN_ROOT="${HICOLA_RUN_ROOT:-/scratch/$USER/Hi-COLA_Output}"
# Reference data and the kept Hi-COLA measurements, next to the repository
export HICOLA_DATA_ROOT="${HICOLA_DATA_ROOT:-$(dirname "$HICOLA_ROOT")/hicola_data}"
export HICOLA_CONDA_ENV="${HICOLA_CONDA_ENV:-hicola}"

export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"

if [ "${CONDA_DEFAULT_ENV:-}" != "$HICOLA_CONDA_ENV" ]; then
    # conda's activation scripts are not written for `set -u`
    case $- in *u*) hicola_nounset=1; set +u ;; *) hicola_nounset= ;; esac
    if [ "$(type -t conda)" != function ]; then
        for conda_exe in "${CONDA_EXE:-}" "$(command -v conda)" "$HOME/miniconda3/bin/conda" "$HOME/anaconda3/bin/conda"; do
            if [ -n "$conda_exe" ] && [ -x "$conda_exe" ]; then
                eval "$("$conda_exe" shell.bash hook)"
                break
            fi
        done
    fi
    conda activate "$HICOLA_CONDA_ENV"
    if [ -n "$hicola_nounset" ]; then set -u; fi
    unset hicola_nounset
fi

# Load the identity and settings of a run file into shell variables (see jobs/run_info.lua),
# and point HICOLA_RUN_ROOT at the run root its output_folder is in:  hicola_load_run <run.lua>
hicola_load_run() {
    local info
    info=$(lua "$HICOLA_ROOT/jobs/run_info.lua" "$1") || return 1
    eval "$info"
    export HICOLA_RUN_ROOT="$RUN_ROOT"
}

# Header of the job logs: code version and where everything is read from and written to
hicola_print_env() {
    local commit
    commit=$(git -C "$HICOLA_ROOT" rev-parse --short HEAD 2>/dev/null || echo unknown)
    if [ -n "$(git -C "$HICOLA_ROOT" status --porcelain --untracked-files=no 2>/dev/null)" ]; then
        commit="$commit + uncommitted changes"
    fi
    echo "================================================================"
    echo "Job:            ${SLURM_JOB_NAME:-interactive} ${SLURM_JOB_ID:-}  ($(date '+%Y-%m-%d %H:%M'))"
    echo "Hi-COLA:        $HICOLA_ROOT  (commit $commit)"
    echo "Conda env:      ${CONDA_DEFAULT_ENV:-none}"
    echo "Run root:       $HICOLA_RUN_ROOT"
    echo "Data root:      $HICOLA_DATA_ROOT"
}

# Settings of the run loaded with hicola_load_run
hicola_print_run() {
    echo "Run file:       $1"
    echo "Run:            $RUN_LABEL  (simulation_name = $SIM_NAME)"
    echo "Box:            $BOXSIZE Mpc/h, ${NPART}^3 particles, $NSTEPS time steps from z = $Z_INIT"
    echo "Gravity:        $GRAVITY"
    echo "Cosmology:      $COSMOLOGY"
    echo "Initial P(k):   $IC_FILE"
    echo "Realization:    seed $IC_SEED, fixed amplitude $IC_FIXED, reversed phases $IC_INVERTED"
    echo "Output z:       $REDSHIFTS"
    echo "Output folder:  $OUTPUT_FOLDER"
}
