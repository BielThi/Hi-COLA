#!/bin/bash
#SBATCH --job-name=hmf_abacus
#SBATCH -o /home/biel/Hi-COLA/HiCOLA/Output/logs/%x_%j.out
#SBATCH -e /home/biel/Hi-COLA/HiCOLA/Output/logs/%x_%j.err
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=0:30:00
#SBATCH --partition=academic
# Halo mass function of the Abacus_Small reference box from its CompaSO halo catalog
# (scripts/run_abacus_hmf.py), written to
# $HICOLA_DATA_ROOT/reference/Abacus_Small/LCDM/Box1/Abacus_Small_LCDM_Box1_z<z>_HMF.txt.
#
# Edit the parameters below, then submit it (from any folder):
#   sbatch ~/Hi-COLA/jobs/hmf_abacus.sh          (or, from the Hi-COLA folder: sbatch jobs/hmf_abacus.sh)

# ---------------------------- user parameters ----------------------------
REDSHIFT="0.500"
MODEL="LCDM"     # only names the output file; the Abacus_Small box is ΛCDM only
CLEANED="1"      # 1: merger-tree cleaned CompaSO catalog (as in the DESI HOD fits), 0: raw
MASS_GRID="1e11 1e16 60"   # MIN MAX N [Msun/h]: the HMF_MASS_GRID of the Abacus_Small runs in halo_finder.sh
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
source "$HICOLA_ROOT/jobs/hicola_env.sh"
read -r -a MASS_GRID_LIST <<< "$MASS_GRID"

hicola_print_env
echo "Measurement:    Abacus_Small $MODEL HMF at z = $REDSHIFT, CompaSO masses (cleaned $CLEANED), mass grid $MASS_GRID"
PYTHONPATH="$HICOLA_ROOT/scripts" python - "$MODEL" "$REDSHIFT" <<'EOF'
import sys
import hicola_paths as hp
model, redshift = sys.argv[1], float(sys.argv[2])
print(f"Reads:          {hp.abacus_halo_dir(redshift)}")
print(f"Writes:         {hp.reference_file('Abacus_Small', model, 1, redshift, 'HMF')}")
EOF
echo "================================================================"

time python "$HICOLA_ROOT/scripts/run_abacus_hmf.py" --redshift "$REDSHIFT" --model "$MODEL" --cleaned "$CLEANED" \
    --mass-grid "${MASS_GRID_LIST[@]}"
