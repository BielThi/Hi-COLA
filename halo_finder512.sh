#!/bin/bash
#SBATCH --job-name=halo_finder
#SBATCH -o %x.out
#SBATCH -e %x.err
#SBATCH -N 8
#SBATCH -n 512
#SBATCH --mem-per-cpu=6G
#SBATCH --cpus-per-task=1
#SBATCH --time=4:30:00
#SBATCH --partition=academic
#SBATCH -D "/home/biel/Hi-COLA/Hi-COLA_Output/logs/"

set -euo pipefail

source /home/biel/miniconda3/etc/profile.d/conda.sh
conda activate hicola

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

# ---------------------------- user parameters ----------------------------
# Edit this block before submitting the job.

MODEL="nDGP_Bart"
REDSHIFT="1.000"

# Hi-COLA Gadget output to process.
SNAPSHOT_DIR="/scratch/biel/${MODEL}/snapshot_${MODEL}_z${REDSHIFT}"
SNAPSHOT_GLOB="${SNAPSHOT_DIR}/gadget_z${REDSHIFT}.*"

# Folder where all halo-finder products will be written.
OUTPUT_DIR="${SNAPSHOT_DIR}/halo_catalog"

# FoF / halo definition settings.
BOXSIZE="1024.0"
LINKING_LENGTH="0.2"
NMIN="30"
MDEF="200c"

# Output filenames.
HALO_CATALOG="${OUTPUT_DIR}/halos200c_catalog.npy"
HMF_OUTPUT="${OUTPUT_DIR}/halos200c_hmf.txt"
HODDIES_OUTPUT="${OUTPUT_DIR}/halos200c_hoddies.npy"

# -------------------------------------------------------------------------

# If this job is submitted from inside an interactive allocation, inherited
# SLURM_MEM_* variables can conflict with this script's --mem-per-cpu request.
unset SLURM_MEM_PER_NODE
unset SLURM_MEM_PER_GPU
unset SLURM_CPU_BIND
unset SLURM_CPU_BIND_LIST
unset SLURM_CPU_BIND_TYPE

echo "Model: ${MODEL}"
echo "Redshift: ${REDSHIFT}"
echo "Snapshot: ${SNAPSHOT_GLOB}"
echo "Halo catalog: ${HALO_CATALOG}"
echo "HMF output: ${HMF_OUTPUT}"
echo "HODDIES output: ${HODDIES_OUTPUT}"

time srun --mpi=pmix --cpu-bind=none -n ${SLURM_NTASKS} /home/biel/miniconda3/envs/hicola/bin/python /home/biel/Hi-COLA/scripts/run_halo_catalog.py \
  --snapshot "${SNAPSHOT_GLOB}" \
  --output "${HALO_CATALOG}" \
  --boxsize "${BOXSIZE}" \
  --redshift "${REDSHIFT}" \
  --linking-length "${LINKING_LENGTH}" \
  --nmin "${NMIN}" \
  --mdef "${MDEF}"

time /home/biel/miniconda3/envs/hicola/bin/python /home/biel/Hi-COLA/scripts/run_halo_outputs.py \
  --halo-catalog "${HALO_CATALOG}" \
  --hmf-output "${HMF_OUTPUT}" \
  --hoddies-output "${HODDIES_OUTPUT}" \
  --boxsize "${BOXSIZE}" \
  --redshift "${REDSHIFT}" \
  --model "${MODEL}"
