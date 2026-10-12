#!/bin/bash
#SBATCH --job-name=frontend
#SBATCH -o HiCOLA/Output/logs/%x_%j.out
#SBATCH -e HiCOLA/Output/logs/%x_%j.err
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=6G
#SBATCH --time=1:00:00
#SBATCH --partition=academic
# Run the Hi-COLA frontend: writes <model>_<cosmology>_{expansion,force}.txt to the
# output_directory of the model .ini (HiCOLA/Output for the files in this repository).
#
# Edit the parameters below, then submit from the Hi-COLA folder:
#   cd ~/Hi-COLA && sbatch jobs/frontend.sh
# It is quick, so it can also run directly:  bash jobs/frontend.sh
# Two .ini files given on the command line replace the ones below.

# ---------------------------- user parameters ----------------------------
# Paths relative to the Hi-COLA folder (or absolute)
MODEL_INI="HiCOLA/Parameter_files/Frontend_parameter_files/examples/cubic_Galileon/cuGal_f0p2.ini"
COSMOLOGY_INI="HiCOLA/Parameter_files/Frontend_parameter_files/examples/WMAP_cosmological_parameters.ini"
# -------------------------------------------------------------------------

set -euo pipefail
if [ -n "${SLURM_JOB_ID:-}" ]; then
    HICOLA_ROOT=${HICOLA_ROOT:-$SLURM_SUBMIT_DIR}   # under SLURM this script runs from a copy
else
    HICOLA_ROOT=${HICOLA_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}
fi
if [ ! -f "$HICOLA_ROOT/jobs/hicola_env.sh" ]; then
    echo "Submit this job from the Hi-COLA folder: cd <Hi-COLA> && sbatch jobs/frontend.sh" >&2
    exit 1
fi
source "$HICOLA_ROOT/jobs/hicola_env.sh"
cd "$HICOLA_ROOT"
MODEL_INI=$(readlink -f "${1:-$MODEL_INI}")
COSMOLOGY_INI=$(readlink -f "${2:-$COSMOLOGY_INI}")

hicola_print_env
echo "Model .ini:     $MODEL_INI"
echo "Cosmology .ini: $COSMOLOGY_INI"
echo "Writes to:      $(sed -n "s/^output_directory *= *['\"]\{0,1\}\([^'\"]*\).*/\1/p" "$MODEL_INI")  (relative to $HICOLA_ROOT)"
echo "================================================================"

# output_directory in the .ini files is relative to the Hi-COLA folder
python3 -m HiCOLA.Frontend.generate_simulation_input "$MODEL_INI" "$COSMOLOGY_INI"
