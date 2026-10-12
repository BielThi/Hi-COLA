#!/bin/bash
#SBATCH --job-name=backend
#SBATCH -o /home/biel/Hi-COLA/HiCOLA/Output/logs/%x_%j.out
#SBATCH -e /home/biel/Hi-COLA/HiCOLA/Output/logs/%x_%j.err
#SBATCH -N 8
#SBATCH -n 512
#SBATCH --mem-per-cpu=6G
#SBATCH --cpus-per-task=1
#SBATCH --time=5:00:00
#SBATCH --partition=academic
#SBATCH --exclusive
# Run one Hi-COLA backend simulation, then copy its P(k) to $HICOLA_DATA_ROOT/hicola/.
#
# Edit the parameters below, then submit it (from any folder):
#   sbatch ~/Hi-COLA/jobs/backend.sh          (or, from the Hi-COLA folder: sbatch jobs/backend.sh)
# A run file given on the command line replaces RUN_FILE:  sbatch jobs/backend.sh <run.lua>
# Logs: the #SBATCH -o/-e lines (HiCOLA/Output/logs/ of this clone; update them in another clone).

# ---------------------------- user parameters ----------------------------
# Run file, relative to the Hi-COLA folder (or absolute). It sets the model, realization,
# output folder and all simulation settings, see PIPELINE.md.
RUN_FILE="HiCOLA/Parameter_files/Backend_parameter_files/Abacus_Small/LCDM_s1234.lua"
# -------------------------------------------------------------------------

set -euo pipefail
# The Hi-COLA folder this script belongs to. SLURM runs a copy of the script, so ask SLURM
# where the original is; outside SLURM, use this file's location.
if [ -n "${SLURM_JOB_ID:-}" ]; then
    SCRIPT=$(scontrol show job "$SLURM_JOB_ID" 2>/dev/null | sed -n 's/^ *Command=\([^ ]*\).*/\1/p')
else
    SCRIPT=${BASH_SOURCE[0]}
fi
HICOLA_ROOT=${HICOLA_ROOT:-$(cd "$(dirname "${SCRIPT:-$SLURM_SUBMIT_DIR/jobs/x}")/.." && pwd)}
if [ ! -f "$HICOLA_ROOT/jobs/hicola_env.sh" ]; then
    echo "Cannot find the Hi-COLA folder (looked in $HICOLA_ROOT)" >&2
    exit 1
fi
# A run file given on the command line is relative to where you submitted; RUN_FILE above is
# relative to the Hi-COLA folder
if [ $# -ge 1 ]; then RUN_FILE=$(cd "${SLURM_SUBMIT_DIR:-$PWD}" && readlink -f "$1"); fi
source "$HICOLA_ROOT/jobs/hicola_env.sh"
cd "$HICOLA_ROOT"
RUN_FILE=$(readlink -f "$RUN_FILE")

hicola_load_run "$RUN_FILE"
hicola_print_env
hicola_print_run "$RUN_FILE"
echo "Resources:      ${SLURM_JOB_NUM_NODES} nodes, ${SLURM_NTASKS} MPI tasks"
echo "================================================================"

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
