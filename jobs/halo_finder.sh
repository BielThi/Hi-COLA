#!/bin/bash
#SBATCH --job-name=halo_finder
#SBATCH -o HiCOLA/Output/logs/%x_%j.out
#SBATCH -e HiCOLA/Output/logs/%x_%j.err
#SBATCH -N 8
#SBATCH -n 512
#SBATCH --mem-per-cpu=6G
#SBATCH --cpus-per-task=1
#SBATCH --time=4:30:00
#SBATCH --partition=academic
# FoF halo finding on one snapshot of a Hi-COLA run, then the HMF and the HODDIES halo array,
# then copy the HMF to $HICOLA_DATA_ROOT/hicola/.
#
# Edit the parameters below, then submit from the Hi-COLA folder:
#   cd ~/Hi-COLA && sbatch jobs/halo_finder.sh
# To start after the backend job: sbatch --dependency=afterok:<backend job id> jobs/halo_finder.sh

# ---------------------------- user parameters ----------------------------
# Run file of the simulation to process, relative to the Hi-COLA folder (or absolute)
RUN_FILE="HiCOLA/Parameter_files/Backend_parameter_files/Abacus_Small/LCDM_s1234.lua"
# Output redshift of the snapshot to process
REDSHIFT="0.500"

# FoF / halo definition. Leave NMIN empty for the target's default (TARGETS in scripts/hicola_paths.py:
# 30 for COLA_NM4096, 50 for Abacus_Small). Box size and cosmology come from the target.
LINKING_LENGTH="0.2"
NMIN=""
MDEF="200c"
# -------------------------------------------------------------------------

set -euo pipefail
HICOLA_ROOT=${HICOLA_ROOT:-$SLURM_SUBMIT_DIR}
if [ ! -f "$HICOLA_ROOT/jobs/hicola_env.sh" ]; then
    echo "Submit this job from the Hi-COLA folder: cd <Hi-COLA> && sbatch jobs/halo_finder.sh" >&2
    exit 1
fi
source "$HICOLA_ROOT/jobs/hicola_env.sh"
cd "$HICOLA_ROOT"
RUN_FILE=$(readlink -f "$RUN_FILE")

eval "$(hicola_run_info "$RUN_FILE")"
if [ "$TARGET" = standalone ]; then
    echo "The halo finder needs a run matched to a reference simulation, not $RUN_LABEL" >&2
    exit 1
fi
RUN=(--target "$TARGET" --model "$MODEL" --realization "$REALIZATION" --redshift "$REDSHIFT")
FOF=(--linking-length "$LINKING_LENGTH" --mdef "$MDEF")
if [ -n "$NMIN" ]; then FOF+=(--nmin "$NMIN"); fi
echo "Run: $RUN_LABEL ($SIM_NAME) at z=$REDSHIFT"

# If this job is submitted from inside an interactive allocation, inherited
# SLURM_MEM_* variables can conflict with this script's --mem-per-cpu request.
unset SLURM_MEM_PER_NODE SLURM_MEM_PER_GPU SLURM_CPU_BIND SLURM_CPU_BIND_LIST SLURM_CPU_BIND_TYPE

time srun --mpi=pmix --cpu-bind=none -n "${SLURM_NTASKS}" \
    python "$HICOLA_ROOT/scripts/run_halo_catalog.py" "${RUN[@]}" "${FOF[@]}"

time python "$HICOLA_ROOT/scripts/run_halo_outputs.py" "${RUN[@]}"

python "$HICOLA_ROOT/scripts/collect_products.py" "$TARGET" "$MODEL" "$REALIZATION"
