#!/bin/bash
# Submit the Abacus_Small reference P(k) measurement (jobs/power_spectrum_abacus.sh).
#
#   jobs/submit_power_spectrum_abacus.sh [sbatch options]
#
# Logs go to HiCOLA/Output/logs/.

set -euo pipefail
source "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/hicola_env.sh"

LOG_DIR="$HICOLA_ROOT/HiCOLA/Output/logs"
mkdir -p "$LOG_DIR"
sbatch --job-name=power_spectrum_abacus --chdir="$LOG_DIR" \
    --output="$LOG_DIR/%x_%j.out" --error="$LOG_DIR/%x_%j.err" \
    "$@" "$HICOLA_ROOT/jobs/power_spectrum_abacus.sh"
