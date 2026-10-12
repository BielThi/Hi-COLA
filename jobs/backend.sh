#!/bin/bash
#SBATCH --job-name=backend
#SBATCH -o HiCOLA/Output/logs/%x_%j.out
#SBATCH -e HiCOLA/Output/logs/%x_%j.err
#SBATCH -N 8
#SBATCH -n 512
#SBATCH --mem-per-cpu=6G
#SBATCH --cpus-per-task=1
#SBATCH --time=5:00:00
#SBATCH --partition=academic
#SBATCH --exclusive
# Run one Hi-COLA backend simulation, then copy its P(k) to $HICOLA_DATA_ROOT/hicola/.
#
# Edit the parameters below, then submit from the Hi-COLA folder:
#   cd ~/Hi-COLA && sbatch jobs/backend.sh
# A run file given on the command line replaces RUN_FILE:  sbatch jobs/backend.sh <run.lua>

# ---------------------------- user parameters ----------------------------
# Run file, relative to the Hi-COLA folder (or absolute). It sets the model, realization,
# output folder and all simulation settings, see PIPELINE.md.
RUN_FILE="HiCOLA/Parameter_files/Backend_parameter_files/Abacus_Small/LCDM_s1234.lua"
# -------------------------------------------------------------------------

set -euo pipefail
RUN_FILE=${1:-$RUN_FILE}
HICOLA_ROOT=${HICOLA_ROOT:-$SLURM_SUBMIT_DIR}
if [ ! -f "$HICOLA_ROOT/jobs/hicola_env.sh" ]; then
    echo "Submit this job from the Hi-COLA folder: cd <Hi-COLA> && sbatch jobs/backend.sh" >&2
    exit 1
fi
source "$HICOLA_ROOT/jobs/hicola_env.sh"
cd "$HICOLA_ROOT"
RUN_FILE=$(readlink -f "$RUN_FILE")

eval "$(hicola_run_info "$RUN_FILE")"
echo "Run:    $RUN_LABEL ($SIM_NAME)"
echo "Config: $RUN_FILE"
echo "Output: $OUTPUT_FOLDER"

# If this job is submitted from inside an interactive allocation, inherited
# SLURM_MEM_* variables can conflict with this script's --mem-per-cpu request.
unset SLURM_MEM_PER_NODE SLURM_MEM_PER_GPU SLURM_CPU_BIND SLURM_CPU_BIND_LIST SLURM_CPU_BIND_TYPE

export TMPDIR=/tmp/${USER}_${SLURM_JOB_ID}
mkdir -p "$TMPDIR"
# Hi-COLA creates only the last level of output_folder, so make the nested run folder here
mkdir -p "$OUTPUT_FOLDER"

time mpirun -np "${SLURM_NTASKS}" "$HICOLA_ROOT/HiCOLA/FML/COLASolver/nbody" "$RUN_FILE"

if [ "$TARGET" != standalone ]; then
    python "$HICOLA_ROOT/scripts/collect_products.py" "$TARGET" "$MODEL" "$REALIZATION"
fi
