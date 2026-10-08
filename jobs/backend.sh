#!/bin/bash
#SBATCH -N 8
#SBATCH -n 512
#SBATCH --mem-per-cpu=6G
#SBATCH --cpus-per-task=1
#SBATCH --time=5:00:00
#SBATCH --partition=academic
#SBATCH --exclusive
# Run one Hi-COLA backend simulation, then copy its P(k) to $HICOLA_DATA_ROOT/hicola/.
# Submit it with jobs/submit_backend.sh <run.lua> [sbatch options].

set -euo pipefail
: "${HICOLA_ROOT:?submit this job with jobs/submit_backend.sh}"
source "$HICOLA_ROOT/jobs/hicola_env.sh"

RUN_FILE=$1
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
