#!/bin/bash
# Submit the halo finder for a Hi-COLA run, one job per redshift.
#
#   jobs/submit_halo_finder.sh <run.lua> [redshift] [sbatch options]
#
# Without a redshift, all output redshifts of the run are submitted.
# e.g.  jobs/submit_halo_finder.sh HiCOLA/Parameter_files/Backend_parameter_files/Abacus_Small/LCDM_s1234.lua 0.5
#
# The FoF settings (box size, reference cosmology, minimum halo size) come from the run's
# target, see TARGETS in scripts/hicola_paths.py.

set -euo pipefail
if [ $# -lt 1 ]; then
    sed -n '2,10p' "$0"
    exit 1
fi
source "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/hicola_env.sh"

RUN_FILE=$(readlink -f "$1")
shift
eval "$(hicola_run_info "$RUN_FILE")"
if [ "$TARGET" = standalone ]; then
    echo "The halo finder needs a run matched to a reference simulation, not $RUN_LABEL" >&2
    exit 1
fi

if [ $# -gt 0 ] && [[ $1 =~ ^[0-9.]+$ ]]; then
    REDSHIFTS=$1
    shift
fi

LOG_DIR="$HICOLA_ROOT/HiCOLA/Output/COLA/$RUN_LABEL/logs"
mkdir -p "$LOG_DIR"
for z in $REDSHIFTS; do
    z=$(printf '%.3f' "$z")
    sbatch --job-name="halos_${SIM_NAME}_z$z" --chdir="$LOG_DIR" \
        --output="$LOG_DIR/%x_%j.out" --error="$LOG_DIR/%x_%j.err" \
        "$@" "$HICOLA_ROOT/jobs/halo_finder.sh" "$TARGET" "$MODEL" "$REALIZATION" "$z"
done
