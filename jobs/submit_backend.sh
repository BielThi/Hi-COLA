#!/bin/bash
# Submit a Hi-COLA backend run.
#
#   jobs/submit_backend.sh <run.lua> [sbatch options]
#
# e.g.  jobs/submit_backend.sh HiCOLA/Parameter_files/Backend_parameter_files/Abacus_Small/LCDM_s1234.lua
#       jobs/submit_backend.sh <run.lua> -N 2 -n 128 --time=1:00:00
#
# Logs go to HiCOLA/Output/COLA/<target>/<model>/<realization>/logs/ (standalone/<name>/logs/
# for the standalone runs).

set -euo pipefail
if [ $# -lt 1 ]; then
    sed -n '2,10p' "$0"
    exit 1
fi
source "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/hicola_env.sh"

RUN_FILE=$(readlink -f "$1")
shift
eval "$(hicola_run_info "$RUN_FILE")"

LOG_DIR="$HICOLA_ROOT/HiCOLA/Output/COLA/$RUN_LABEL/logs"
mkdir -p "$LOG_DIR"
sbatch --job-name="$SIM_NAME" --chdir="$LOG_DIR" \
    --output="$LOG_DIR/%x_%j.out" --error="$LOG_DIR/%x_%j.err" \
    "$@" "$HICOLA_ROOT/jobs/backend.sh" "$RUN_FILE"
