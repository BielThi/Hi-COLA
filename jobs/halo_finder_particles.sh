#!/bin/bash
#SBATCH --job-name=halo_finder_particles
#SBATCH -o /home/biel/dev-v2/Hi-COLA/HiCOLA/Output/logs/%x_%j.out
#SBATCH -e /home/biel/dev-v2/Hi-COLA/HiCOLA/Output/logs/%x_%j.err
#SBATCH -N 8
#SBATCH -n 512
#SBATCH --mem-per-cpu=6G
#SBATCH --cpus-per-task=1
#SBATCH --time=4:30:00
#SBATCH --partition=academic
# Tryout: FoF halo finding on one snapshot of a Hi-COLA run (as jobs/halo_finder.sh) that also
# keeps a random subsample of the halo member particles, for HODDIES' use_particles (satellites on
# particles instead of an NFW profile; jobs/hod_fit_particles.sh). Writes, next to the main halo
# catalog and without touching it, halo_catalog/halos200c_ptcl_{catalog,hoddies,particles}.npy
# (same halo_id in the three; Rh comoving in the HODDIES array). scripts/run_halo_particles.py.
#
# Edit the parameters below, then submit it (from any folder):
#   sbatch ~/dev-v2/Hi-COLA/jobs/halo_finder_particles.sh

# ---------------------------- user parameters ----------------------------
# Run file of the simulation to process, relative to the Hi-COLA folder (or absolute)
RUN_FILE="HiCOLA/Parameter_files/Backend_parameter_files/Abacus_Small/LCDM_s1234.lua"
# Output redshift of the snapshot to process
REDSHIFT="0.500"

# FoF / halo definition
LINKING_LENGTH="0.2"   # in units of the mean particle separation
NMIN="50"              # minimum particles per halo; so far 30 for COLA_NM4096 runs, 50 for Abacus_Small runs
MDEF="200c"            # mass definition; only 200c for now: the catalog fields (m200c, r200c) and files (halos200c_*) are named for it

# Member particles kept: random fraction of every halo's particles, and the seed of the draw
# (10% keeps ~100 particles in a 1e13 Msun/h halo for 1024^3 particles in 500 Mpc/h)
PARTICLE_FRACTION="0.1"
PARTICLE_SEED="1"
# -------------------------------------------------------------------------
# Box size and cosmology (for the halo masses) are taken from the run file.

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
cd "$HICOLA_ROOT"
RUN_FILE=$(readlink -f "$RUN_FILE")

hicola_load_run "$RUN_FILE"
if [ "$TARGET" = standalone ]; then
    echo "The halo finder needs a run matched to a reference simulation, not $RUN_LABEL" >&2
    exit 1
fi
RUN=(--target "$TARGET" --model "$MODEL" --realization "$REALIZATION" --redshift "$REDSHIFT" --boxsize "$BOXSIZE")
FOF=(--linking-length "$LINKING_LENGTH" --nmin "$NMIN" --mdef "$MDEF")
COSMO_ARGS=(--h "$COSMO_H" --omega-cdm "$COSMO_OMEGA_CDM" --omega-b "$COSMO_OMEGA_B" --omega-mnu "$COSMO_OMEGA_MNU"
           --omega-lambda "$COSMO_OMEGA_LAMBDA" --tcmb "$COSMO_TCMB" --neff "$COSMO_NEFF" --n-s "$COSMO_NS")

hicola_print_env
hicola_print_run "$RUN_FILE"
echo "Snapshot z:     $REDSHIFT"
echo "FoF:            linking length $LINKING_LENGTH, at least $NMIN particles, mass definition $MDEF"
echo "Particles:      fraction $PARTICLE_FRACTION of the member particles (seed $PARTICLE_SEED)"
echo "Resources:      ${SLURM_JOB_NUM_NODES} nodes, ${SLURM_NTASKS} MPI tasks"
echo "================================================================"

# If this job is submitted from inside an interactive allocation, inherited
# SLURM_MEM_* variables can conflict with this script's --mem-per-cpu request.
unset SLURM_MEM_PER_NODE SLURM_MEM_PER_GPU SLURM_CPU_BIND SLURM_CPU_BIND_LIST SLURM_CPU_BIND_TYPE

time srun --mpi=pmix --cpu-bind=none -n "${SLURM_NTASKS}" \
    python "$HICOLA_ROOT/scripts/run_halo_particles.py" "${RUN[@]}" "${FOF[@]}" "${COSMO_ARGS[@]}" \
    --particle-fraction "$PARTICLE_FRACTION" --particle-seed "$PARTICLE_SEED"
