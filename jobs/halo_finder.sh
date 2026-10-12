#!/bin/bash
#SBATCH --job-name=halo_finder
#SBATCH -o HiCOLA/Output/logs/%x_%j.out
#SBATCH -e HiCOLA/Output/logs/%x_%j.err
#SBATCH -N 8
#SBATCH -n 512
#SBATCH --mem-per-cpu=6G
#SBATCH --cpus-per-task=1
#SBATCH --time=4:30:00
#SBATCH --partition=academic
# FoF halo finding on one snapshot of a Hi-COLA run, then the HMF and the HODDIES halo array,
# then copy the HMF to $HICOLA_DATA_ROOT/hicola/.
#
# Edit the parameters below, then submit from the Hi-COLA folder:
#   cd ~/Hi-COLA && sbatch jobs/halo_finder.sh
# To start after the backend job: sbatch --dependency=afterok:<backend job id> jobs/halo_finder.sh

# ---------------------------- user parameters ----------------------------
# Run file of the simulation to process, relative to the Hi-COLA folder (or absolute)
RUN_FILE="HiCOLA/Parameter_files/Backend_parameter_files/Abacus_Small/LCDM_s1234.lua"
# Output redshift of the snapshot to process
REDSHIFT="0.500"

# FoF / halo definition
LINKING_LENGTH="0.2"   # in units of the mean particle separation
NMIN="50"              # minimum particles per halo; so far 30 for COLA_NM4096 runs, 50 for Abacus_Small runs
MDEF="200c"            # mass definition; only 200c for now: the catalog fields (m200c, r200c) and files (halos200c_*) are named for it

# HMF mass grid: "reference" = the masses of the target's reference HMF at this redshift, to
# compare bin by bin (COLA_NM4096 only, so far); or "MIN MAX N" = N log-spaced masses in Msun/h
HMF_MASS_GRID="1e11 1e16 60"
# -------------------------------------------------------------------------
# Box size and cosmology (for the halo masses) are taken from the run file.

set -euo pipefail
HICOLA_ROOT=${HICOLA_ROOT:-$SLURM_SUBMIT_DIR}
if [ ! -f "$HICOLA_ROOT/jobs/hicola_env.sh" ]; then
    echo "Submit this job from the Hi-COLA folder: cd <Hi-COLA> && sbatch jobs/halo_finder.sh" >&2
    exit 1
fi
source "$HICOLA_ROOT/jobs/hicola_env.sh"
cd "$HICOLA_ROOT"
RUN_FILE=$(readlink -f "$RUN_FILE")

hicola_load_run "$RUN_FILE"
if [ "$TARGET" = standalone ]; then
    echo "The halo finder needs a run matched to a reference simulation, not $RUN_LABEL" >&2
    exit 1
fi
RUN=(--target "$TARGET" --model "$MODEL" --realization "$REALIZATION" --redshift "$REDSHIFT" --boxsize "$BOXSIZE")
FOF=(--linking-length "$LINKING_LENGTH" --nmin "$NMIN" --mdef "$MDEF")
COSMOLOGY=(--h "$COSMO_H" --omega-cdm "$COSMO_OMEGA_CDM" --omega-b "$COSMO_OMEGA_B" --omega-mnu "$COSMO_OMEGA_MNU"
           --omega-lambda "$COSMO_OMEGA_LAMBDA" --tcmb "$COSMO_TCMB" --neff "$COSMO_NEFF" --n-s "$COSMO_NS")
read -r -a GRID <<< "$HMF_MASS_GRID"

hicola_print_env
hicola_print_run "$RUN_FILE"
echo "Snapshot z:     $REDSHIFT"
echo "FoF:            linking length $LINKING_LENGTH, at least $NMIN particles, mass definition $MDEF"
echo "Halo masses:    computed with the run's cosmology (above)"
PYTHONPATH="$HICOLA_ROOT/scripts" python - "$TARGET" "$REDSHIFT" "${GRID[@]}" <<'EOF'
import sys
from run_halo_outputs import mass_grid
print(f"HMF mass grid:  {mass_grid(sys.argv[3:], sys.argv[1], float(sys.argv[2]))[1]}")
EOF
echo "Resources:      ${SLURM_JOB_NUM_NODES} nodes, ${SLURM_NTASKS} MPI tasks"
echo "================================================================"

# If this job is submitted from inside an interactive allocation, inherited
# SLURM_MEM_* variables can conflict with this script's --mem-per-cpu request.
unset SLURM_MEM_PER_NODE SLURM_MEM_PER_GPU SLURM_CPU_BIND SLURM_CPU_BIND_LIST SLURM_CPU_BIND_TYPE

time srun --mpi=pmix --cpu-bind=none -n "${SLURM_NTASKS}" \
    python "$HICOLA_ROOT/scripts/run_halo_catalog.py" "${RUN[@]}" "${FOF[@]}" "${COSMOLOGY[@]}"

time python "$HICOLA_ROOT/scripts/run_halo_outputs.py" "${RUN[@]}" --mass-grid "${GRID[@]}"

python "$HICOLA_ROOT/scripts/collect_products.py" "$TARGET" "$MODEL" "$REALIZATION"
