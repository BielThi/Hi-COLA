#!/bin/bash
#SBATCH --job-name=hod_fof
#SBATCH -o %x.out
#SBATCH -e %x.err
#SBATCH -N 8
#SBATCH -n 512
#SBATCH --mem-per-cpu=6G
#SBATCH --cpus-per-task=1
#SBATCH --time=4:00:00
#SBATCH --partition=academic
#SBATCH -D "/home/biel/Hi-COLA/Hi-COLA_Output/logs/"

source /home/biel/miniconda3/etc/profile.d/conda.sh
conda activate hicola

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

# If this job is submitted from inside an interactive allocation, inherited
# SLURM_MEM_* variables can conflict with this script's --mem-per-cpu request.
unset SLURM_MEM_PER_NODE
unset SLURM_MEM_PER_GPU
unset SLURM_CPU_BIND
unset SLURM_CPU_BIND_LIST
unset SLURM_CPU_BIND_TYPE

time srun --mpi=pmix --cpu-bind=none -n ${SLURM_NTASKS} /home/biel/miniconda3/envs/hicola/bin/python /home/biel/Hi-COLA/scripts/run_hod_fof.py \
  --snapshot "/scratch/biel/LCDM/snapshot_LCDM_z0.506/gadget_z0.506.*" \
  --output "/scratch/biel/LCDM/snapshot_LCDM_z0.506/hod_fof/halos200c_hod_arr.npy" \
  --boxsize 1024.0 \
  --redshift 0.506 \
  --linking-length 0.2 \
  --nmin 30 \
  --mdef 200c
