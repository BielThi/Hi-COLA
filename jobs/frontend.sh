#!/bin/bash
#SBATCH --job-name=frontend
#SBATCH -o /home/biel/dev-v2/Hi-COLA/HiCOLA/Output/logs/%x_%j.out
#SBATCH -e /home/biel/dev-v2/Hi-COLA/HiCOLA/Output/logs/%x_%j.err
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=6G
#SBATCH --time=1:00:00
#SBATCH --partition=academic
# Run the Hi-COLA frontend: writes <model>_<cosmology>_{expansion,force}.txt to the
# output_directory of the model .ini (HiCOLA/Output for the files in this repository).
#
# Edit the parameters below, then submit it (from any folder):
#   sbatch ~/Hi-COLA/jobs/frontend.sh          (or, from the Hi-COLA folder: sbatch jobs/frontend.sh)
# It is quick, so it can also run directly:  bash jobs/frontend.sh
# Two .ini files given on the command line replace the ones below.

# ---------------------------- user parameters ----------------------------
# Paths relative to the Hi-COLA folder (or absolute)
MODEL_INI="HiCOLA/Parameter_files/Frontend_parameter_files/examples/cubic_Galileon/cuGal_f0p2.ini"
COSMOLOGY_INI="HiCOLA/Parameter_files/Frontend_parameter_files/examples/WMAP_cosmological_parameters.ini"
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
# .ini files given on the command line are relative to where you ran/submitted it; the ones
# above are relative to the Hi-COLA folder
if [ $# -ge 2 ]; then
    MODEL_INI=$(cd "${SLURM_SUBMIT_DIR:-$PWD}" && readlink -f "$1")
    COSMOLOGY_INI=$(cd "${SLURM_SUBMIT_DIR:-$PWD}" && readlink -f "$2")
fi
source "$HICOLA_ROOT/jobs/hicola_env.sh"
cd "$HICOLA_ROOT"
MODEL_INI=$(readlink -f "$MODEL_INI")
COSMOLOGY_INI=$(readlink -f "$COSMOLOGY_INI")

hicola_print_env
echo "Model .ini:     $MODEL_INI"
echo "Cosmology .ini: $COSMOLOGY_INI"
echo "Writes to:      $(sed -n "s/^output_directory *= *['\"]\{0,1\}\([^'\"]*\).*/\1/p" "$MODEL_INI")  (relative to $HICOLA_ROOT)"
echo "================================================================"

# output_directory in the .ini files is relative to the Hi-COLA folder
python3 -m HiCOLA.Frontend.generate_simulation_input "$MODEL_INI" "$COSMOLOGY_INI"
