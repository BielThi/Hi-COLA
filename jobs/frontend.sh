#!/bin/bash
# Run the Hi-COLA frontend: writes <model>_<cosmology>_{expansion,force}.txt to the
# output_directory of the model .ini (HiCOLA/Output for the files in this repository).
#
#   jobs/frontend.sh <model.ini> <cosmology.ini>
#
# e.g.  jobs/frontend.sh HiCOLA/Parameter_files/Frontend_parameter_files/nDGP/nDGP.ini \
#           HiCOLA/Parameter_files/Frontend_parameter_files/examples/Planck18_cosmological_parameters.ini

set -euo pipefail
if [ $# -ne 2 ]; then
    sed -n '2,9p' "$0"
    exit 1
fi
MODEL_INI=$(readlink -f "$1")
COSMOLOGY_INI=$(readlink -f "$2")
source "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/hicola_env.sh"

# output_directory in the .ini files is relative to the repository root
cd "$HICOLA_ROOT"
python3 -m HiCOLA.Frontend.generate_simulation_input "$MODEL_INI" "$COSMOLOGY_INI"
