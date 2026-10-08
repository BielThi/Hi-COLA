#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=16
#SBATCH --mem=40G
#SBATCH --time=1:30:00
#SBATCH --partition=academic
# Matter P(k) of the Abacus_Small reference box, written to
# $HICOLA_DATA_ROOT/reference/Abacus_Small/LCDM/Box1/Abacus_Small_LCDM_Box1_z<z>_PS.txt.
# Submit it with jobs/submit_power_spectrum_abacus.sh [sbatch options].

set -euo pipefail
: "${HICOLA_ROOT:?submit this job with jobs/submit_power_spectrum_abacus.sh}"
source "$HICOLA_ROOT/jobs/hicola_env.sh"

# ---------------------------- user parameters ----------------------------
REDSHIFT="0.500"
MODEL="LCDM"

LOAD_A=1  # load subsample A (3%)
LOAD_B=0  # load subsample B (7%); combine with A for the full 10% subsample

NMESH="1024"  # was 512; k_max = pi*NMESH/Lbox, raised to push the Nyquist cutoff further out
KBINS="60"    # independent of NMESH: fixes the mostly-empty/noisy low-k bins from before
# -------------------------------------------------------------------------

time python "$HICOLA_ROOT/scripts/run_abacus_power_spectrum.py" \
  --redshift "${REDSHIFT}" \
  --load-A "${LOAD_A}" \
  --load-B "${LOAD_B}" \
  --nmesh "${NMESH}" \
  --kbins "${KBINS}" \
  --nthread "${SLURM_CPUS_PER_TASK}" \
  --model "${MODEL}"
