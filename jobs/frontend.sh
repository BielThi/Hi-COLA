#!/bin/bash
# Run the Hi-COLA frontend: writes <model>_<cosmology>_{expansion,force}.txt to the
# output_directory of the model .ini (HiCOLA/Output for the files in this repository).
# It is quick, so it runs directly rather than as a SLURM job.
#
# Edit the parameters below, then run:   bash jobs/frontend.sh
# Two .ini files given on the command line replace the ones below.

# ---------------------------- user parameters ----------------------------
# Paths relative to the Hi-COLA folder (or absolute)
MODEL_INI="HiCOLA/Parameter_files/Frontend_parameter_files/examples/cubic_Galileon/cuGal_f0p2.ini"
COSMOLOGY_INI="HiCOLA/Parameter_files/Frontend_parameter_files/examples/WMAP_cosmological_parameters.ini"
# -------------------------------------------------------------------------

set -euo pipefail
source "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/hicola_env.sh"
cd "$HICOLA_ROOT"
MODEL_INI=$(readlink -f "${1:-$MODEL_INI}")
COSMOLOGY_INI=$(readlink -f "${2:-$COSMOLOGY_INI}")

# output_directory in the .ini files is relative to the Hi-COLA folder
python3 -m HiCOLA.Frontend.generate_simulation_input "$MODEL_INI" "$COSMOLOGY_INI"
