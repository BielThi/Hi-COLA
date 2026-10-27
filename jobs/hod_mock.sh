#!/bin/bash
#SBATCH --job-name=hod_mock
#SBATCH -o /home/biel/Hi-COLA/HiCOLA/Output/logs/%x_%j.out
#SBATCH -e /home/biel/Hi-COLA/HiCOLA/Output/logs/%x_%j.err
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=6G
#SBATCH --time=1:00:00
#SBATCH --partition=academic
# HODDIES galaxy mock on the halo array of one snapshot of a Hi-COLA run, and its redshift-space
# power spectrum and marked power spectrum multipoles measured as in the EMC/acm pipeline
# (scripts/run_hod_mock.py, scripts/galaxy_power_spectrum.py), then copied to $HICOLA_DATA_ROOT/hicola/.
#
# Edit the parameters below, then submit it (from any folder):
#   sbatch ~/Hi-COLA/jobs/hod_mock.sh          (or, from the Hi-COLA folder: sbatch jobs/hod_mock.sh)
# To start after the halo finder: sbatch --dependency=afterok:<halo_finder job id> jobs/hod_mock.sh

# ---------------------------- user parameters ----------------------------
# Run file of the simulation, relative to the Hi-COLA folder (or absolute)
RUN_FILE="HiCOLA/Parameter_files/Backend_parameter_files/Abacus_Small/nDGP_N1.2_s1234.lua"
# Output redshift of the snapshot (its halo catalog must exist: jobs/halo_finder.sh)
REDSHIFT="0.500"

# HOD: the best fit of a fit of jobs/hod_fit.sh (its name, line 'Fit:' of its log), made on the run of
# HOD_FIT_MODEL at the same target, realization and redshift: "own" = this run's model (each model
# with its own fit), or a model, e.g. "LCDM" (the LCDM fit for every model). The fits fix the galaxy
# density: HODDIES rescales Ac to reach it.
HOD_FIT="LRG_SHOD_6p_ptcl_r1-max"
HOD_FIT_MODEL="own"
TRACER="LRG"       # tracer of the fit
# Satellites: "particles" (random member particles of the halo: jobs/halo_finder.sh with PARTICLE_FRACTION > 0)
# or "nfw" (HODDIES' NFW profile); must be those of the fit
SATELLITES="particles"
SEEDS="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20"      # HOD seeds, one mock, P(k) and marked P(k) each, e.g. "1 2 3 4 5"; ~1.5 min per seed at mesh 512 (raise --time above for more than ~35 seeds)
NTHREADS="32"   # HODDIES threads; the mock of a seed depends on it: keep it fixed across seeds

# Power spectrum (EMC/acm settings so far: z, 512, 0 2 4, 0.001)
LOS="z"          # line of sight (box axis) of the redshift-space shift
MESHSIZE="512"
ELLS="0 2 4"     # multipoles
K_STEP="0.001"   # k-bin width [h/Mpc]

# Marked power spectrum (same settings as the P(k) above): mark m = sum_n c_n delta_R^n, with delta_R
# the galaxy density contrast smoothed by a Gaussian of radius R, normalised by its mean over the galaxies
MARK_RCUT="25"          # R [Mpc/h], e.g. "10"
MARK_COEFFICIENTS="1 -1"  # c_0 c_1 c_2 ..., e.g. "1 -0.5"

# Output: halo_catalog/hod/<HOD_FIT>_fit-<model of the fit>_seed<seed>_PS.txt and _MPS.txt in the run folder
# -------------------------------------------------------------------------
# Box size and cosmology (for H(z) in the redshift-space shift) are taken from the run file.

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
    echo "The HOD mock needs a run matched to a reference simulation, not $RUN_LABEL" >&2
    exit 1
fi
FIT_MODEL=$([ "$HOD_FIT_MODEL" = own ] && echo "$MODEL" || echo "$HOD_FIT_MODEL")
NAME="${HOD_FIT}_fit-${FIT_MODEL}_seed{seed}"
RUN=(--target "$TARGET" --model "$MODEL" --realization "$REALIZATION" --redshift "$REDSHIFT" --boxsize "$BOXSIZE")
COSMO_ARGS=(--h "$COSMO_H" --omega-cdm "$COSMO_OMEGA_CDM" --omega-b "$COSMO_OMEGA_B" --omega-mnu "$COSMO_OMEGA_MNU"
           --tcmb "$COSMO_TCMB" --neff "$COSMO_NEFF")
read -r -a SEED_LIST <<< "$SEEDS"
MOCK=(--hod-fit "$HOD_FIT" --hod-fit-model "$FIT_MODEL" --tracer "$TRACER" --seeds "${SEED_LIST[@]}" --name "$NAME" --nthreads "$NTHREADS"
      --satellites "$SATELLITES")
read -r -a ELL_LIST <<< "$ELLS"
PK=(--los "$LOS" --meshsize "$MESHSIZE" --ells "${ELL_LIST[@]}" --k-step "$K_STEP")
if [ -z "$MARK_RCUT" ] || [ -z "$MARK_COEFFICIENTS" ]; then
    echo "Set MARK_RCUT and MARK_COEFFICIENTS in the user parameters of $SCRIPT" >&2
    exit 1
fi
read -r -a MARK_COEFFICIENT_LIST <<< "$MARK_COEFFICIENTS"
MARK=(--mark-rcut "$MARK_RCUT" --mark-coefficients "${MARK_COEFFICIENT_LIST[@]}")

hicola_print_env
hicola_print_run "$RUN_FILE"
echo "Snapshot z:     $REDSHIFT"
echo "HOD:            $TRACER, best fit of $HOD_FIT on $FIT_MODEL, seeds $SEEDS, $NTHREADS threads"
echo "Satellites:     $SATELLITES"
echo "P(k):           redshift space along $LOS, mesh ${MESHSIZE}^3, ells $ELLS, k-step $K_STEP (jaxpower, TSC, interlacing 3)"
echo "Marked P(k):    mark sum_n c_n delta_R^n, R = $MARK_RCUT Mpc/h, coefficients $MARK_COEFFICIENTS"
echo "Output names:   ${NAME}_PS.txt and ${NAME}_MPS.txt for each seed"
echo "Resources:      ${SLURM_CPUS_PER_TASK:-?} CPUs"
echo "================================================================"

export NUMBA_NUM_THREADS="$NTHREADS"
time python "$HICOLA_ROOT/scripts/run_hod_mock.py" "${RUN[@]}" "${COSMO_ARGS[@]}" "${MOCK[@]}" "${PK[@]}" "${MARK[@]}"

python "$HICOLA_ROOT/scripts/collect_products.py" "$TARGET" "$MODEL" "$REALIZATION"
