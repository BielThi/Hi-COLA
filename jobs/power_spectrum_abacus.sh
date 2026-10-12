#!/bin/bash
#SBATCH --job-name=power_spectrum_abacus
#SBATCH -o HiCOLA/Output/logs/%x_%j.out
#SBATCH -e HiCOLA/Output/logs/%x_%j.err
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=16
#SBATCH --mem=40G
#SBATCH --time=1:30:00
#SBATCH --partition=academic
# Matter P(k) of the Abacus_Small reference box, written to
# $HICOLA_DATA_ROOT/reference/Abacus_Small/LCDM/Box1/Abacus_Small_LCDM_Box1_z<z>_PS.txt.
#
# Edit the parameters below, then submit from the Hi-COLA folder:
#   cd ~/Hi-COLA && sbatch jobs/power_spectrum_abacus.sh

# ---------------------------- user parameters ----------------------------
# The estimator itself (TSC mass assignment, interlacing, shot-noise subtraction) is fixed in
# scripts/run_abacus_power_spectrum.py.
REDSHIFT="0.500"
MODEL="LCDM"

LOAD_A=1  # load subsample A (3%)
LOAD_B=0  # load subsample B (7%); combine with A for the full 10% subsample

NMESH="1024"  # was 512; k_max = pi*NMESH/Lbox, raised to push the Nyquist cutoff further out
KBINS="60"    # independent of NMESH: fixes the mostly-empty/noisy low-k bins from before
# -------------------------------------------------------------------------

set -euo pipefail
HICOLA_ROOT=${HICOLA_ROOT:-$SLURM_SUBMIT_DIR}
if [ ! -f "$HICOLA_ROOT/jobs/hicola_env.sh" ]; then
    echo "Submit this job from the Hi-COLA folder: cd <Hi-COLA> && sbatch jobs/power_spectrum_abacus.sh" >&2
    exit 1
fi
source "$HICOLA_ROOT/jobs/hicola_env.sh"

hicola_print_env
echo "Measurement:    Abacus_Small $MODEL P(k) at z = $REDSHIFT, subsamples A=$LOAD_A B=$LOAD_B, Nmesh $NMESH, $KBINS k-bins"
echo "Estimator:      halo + field particles, TSC, interlaced, window compensated, shot noise subtracted"
echo "                (fixed in scripts/run_abacus_power_spectrum.py)"
PYTHONPATH="$HICOLA_ROOT/scripts" python - "$MODEL" "$REDSHIFT" <<'EOF'
import sys
import hicola_paths as hp
model, redshift = sys.argv[1], float(sys.argv[2])
print(f"Reads:          {hp.abacus_halo_dir(redshift)}")
print(f"Writes:         {hp.reference_file('Abacus_Small', model, 1, redshift, 'PS')}")
EOF
echo "Resources:      ${SLURM_CPUS_PER_TASK} threads"
echo "================================================================"

time python "$HICOLA_ROOT/scripts/run_abacus_power_spectrum.py" \
  --redshift "${REDSHIFT}" \
  --load-A "${LOAD_A}" \
  --load-B "${LOAD_B}" \
  --nmesh "${NMESH}" \
  --kbins "${KBINS}" \
  --nthread "${SLURM_CPUS_PER_TASK}" \
  --model "${MODEL}"
