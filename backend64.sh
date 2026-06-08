#!/bin/bash
#SBATCH --job-name=cuGal64
#SBATCH -o %x.out
#SBATCH -e %x.err
#SBATCH -N 1
#SBATCH -n 64
#SBATCH --mem-per-cpu=6G
#SBATCH --cpus-per-task=1
#SBATCH --time=1:00:00
#SBATCH --partition=academic
#SBATCH -D "/home/biel/Hi-COLA/Hi-COLA_Output/logs/"

source /home/biel/miniconda3/etc/profile.d/conda.sh
conda activate hicola

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

time mpirun -np ${SLURM_NTASKS} /home/biel/Hi-COLA/HiCOLA/FML/COLASolver/nbody /home/biel/Hi-COLA/HiCOLA/Parameter_files/Backend_parameter_files/cuGal.lua