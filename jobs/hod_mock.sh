#!/bin/bash
#SBATCH --job-name=hod_mock
#SBATCH -o /global/homes/b/bielthi/Hi-COLA/HiCOLA/Output/logs/%x_%j.out
#SBATCH -e /global/homes/b/bielthi/Hi-COLA/HiCOLA/Output/logs/%x_%j.err
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=32
#SBATCH --time=0:30:00
#SBATCH --constraint=cpu
#SBATCH --qos=shared
#SBATCH --account=desi
# HODDIES galaxy mock on the halo array of one snapshot of a Hi-COLA run, written in the format
# of the EMC AbacusHOD mocks (FITS, see scripts/run_hod_mock.py), then copied to
# $HICOLA_DATA_ROOT/hicola/.
#
# Edit the parameters below, then submit it (from any folder):
#   sbatch ~/Hi-COLA/jobs/hod_mock.sh          (or, from the Hi-COLA folder: sbatch jobs/hod_mock.sh)
# To start after the halo finder: sbatch --dependency=afterok:<halo_finder job id> jobs/hod_mock.sh

# ---------------------------- user parameters ----------------------------
# Run file of the simulation, relative to the Hi-COLA folder (or absolute)
RUN_FILE="HiCOLA/Parameter_files/Backend_parameter_files/Abacus_Small/nDGP_N1.2_s1234.lua"
# Output redshift of the snapshot (its halo catalog must exist: jobs/halo_finder.sh)
REDSHIFT="0.500"

# HODDIES parameter file (relative to the Hi-COLA folder) and the tracer of it to populate.
# If the file sets the tracer's `density`, HODDIES rescales Ac to reach it.
HOD_PARAM_FILE="Notebooks/HOD_parameter_file_z0.500.yaml"
TRACER="LRG"
SEED="42"
NTHREADS="32"   # HODDIES threads; the mock of a seed depends on it (at most --cpus-per-task)

# Output: halo_catalog/hod/<NAME>.fits in the run folder
NAME="${TRACER}_seed${SEED}"
# -------------------------------------------------------------------------
# Box size is taken from the run file.

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
source "$HICOLA_ROOT/jobs/hicola_env.sh"
cd "$HICOLA_ROOT"
RUN_FILE=$(readlink -f "$RUN_FILE")
HOD_PARAM_FILE=$(readlink -f "$HOD_PARAM_FILE")

hicola_load_run "$RUN_FILE"
if [ "$TARGET" = standalone ]; then
    echo "The HOD mock needs a run matched to a reference simulation, not $RUN_LABEL" >&2
    exit 1
fi
RUN=(--target "$TARGET" --model "$MODEL" --realization "$REALIZATION" --redshift "$REDSHIFT" --boxsize "$BOXSIZE")
MOCK=(--hod-param-file "$HOD_PARAM_FILE" --tracer "$TRACER" --seed "$SEED" --name "$NAME" --nthreads "$NTHREADS")

hicola_print_env
hicola_print_run "$RUN_FILE"
echo "Snapshot z:     $REDSHIFT"
echo "HOD:            $TRACER from $HOD_PARAM_FILE, seed $SEED, $NTHREADS threads"
echo "Output name:    $NAME"
echo "Resources:      ${SLURM_CPUS_PER_TASK:-?} CPUs"
echo "================================================================"

export NUMBA_NUM_THREADS="$NTHREADS"
time python "$HICOLA_ROOT/scripts/run_hod_mock.py" "${RUN[@]}" "${MOCK[@]}"

python "$HICOLA_ROOT/scripts/collect_products.py" "$TARGET" "$MODEL" "$REALIZATION"
