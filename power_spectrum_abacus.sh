#!/bin/bash
#SBATCH --job-name=power_spectrum_abacus
#SBATCH -o %x.out
#SBATCH -e %x.err
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=16
#SBATCH --mem=40G
#SBATCH --time=1:30:00
#SBATCH --partition=academic
#SBATCH -D "/home/biel/Hi-COLA/Hi-COLA_Output/logs/"

set -euo pipefail

source /home/biel/miniconda3/etc/profile.d/conda.sh
conda activate hicola

# ---------------------------- user parameters ----------------------------
# Edit this block before submitting the job.

BOX="Abacus1"
SIM_NAME="AbacusSummit_small_c000_ph3000"
REDSHIFT="z0.500"
MODEL="LCDM"

LOAD_A=1  # load subsample A (3%)
LOAD_B=0  # load subsample B (7%); combine with A for the full 10% subsample

NMESH="1024"  # was 512; k_max = pi*NMESH/Lbox, raised to push the Nyquist cutoff further out
KBINS="60"    # independent of NMESH: fixes the mostly-empty/noisy low-k bins from before

PATH_TO_SIM="/home/biel/Hi-COLA/Comparison/${BOX}/small/${SIM_NAME}/halos/${REDSHIFT}"
OUTPUT_DIR="/home/biel/Hi-COLA/Comparison/${BOX}"

# -------------------------------------------------------------------------

echo "Box: ${BOX}"
echo "Sim: ${SIM_NAME} at ${REDSHIFT}"
echo "Path to sim: ${PATH_TO_SIM}"
echo "Output dir: ${OUTPUT_DIR}"

time /home/biel/miniconda3/envs/hicola/bin/python /home/biel/Hi-COLA/scripts/run_abacus_power_spectrum.py \
  --path-to-sim "${PATH_TO_SIM}" \
  --load-A "${LOAD_A}" \
  --load-B "${LOAD_B}" \
  --nmesh "${NMESH}" \
  --kbins "${KBINS}" \
  --nthread "${SLURM_CPUS_PER_TASK}" \
  --model "${MODEL}" \
  --output-dir "${OUTPUT_DIR}"
