#!/bin/bash
#SBATCH -N 8
#SBATCH -n 512
#SBATCH --mem-per-cpu=6G
#SBATCH --cpus-per-task=1
#SBATCH --time=4:30:00
#SBATCH --partition=academic
# FoF halo finding on one snapshot of a Hi-COLA run, then the HMF and HODDIES outputs, then
# copy the HMF to $HICOLA_DATA_ROOT/hicola/.
# Submit it with jobs/submit_halo_finder.sh <run.lua> [redshift] [sbatch options].

set -euo pipefail
: "${HICOLA_ROOT:?submit this job with jobs/submit_halo_finder.sh}"
source "$HICOLA_ROOT/jobs/hicola_env.sh"

TARGET=$1
MODEL=$2
REALIZATION=$3
REDSHIFT=$4
RUN=(--target "$TARGET" --model "$MODEL" --realization "$REALIZATION" --redshift "$REDSHIFT")
echo "Run: $TARGET/$MODEL/$REALIZATION at z=$REDSHIFT"

# If this job is submitted from inside an interactive allocation, inherited
# SLURM_MEM_* variables can conflict with this script's --mem-per-cpu request.
unset SLURM_MEM_PER_NODE SLURM_MEM_PER_GPU SLURM_CPU_BIND SLURM_CPU_BIND_LIST SLURM_CPU_BIND_TYPE

time srun --mpi=pmix --cpu-bind=none -n "${SLURM_NTASKS}" \
    python "$HICOLA_ROOT/scripts/run_halo_catalog.py" "${RUN[@]}"

time python "$HICOLA_ROOT/scripts/run_halo_outputs.py" "${RUN[@]}"

python "$HICOLA_ROOT/scripts/collect_products.py" "$TARGET" "$MODEL" "$REALIZATION"
