#!/bin/bash
#SBATCH --job-name=hod_fit
#SBATCH -o /home/biel/dev-v2/Hi-COLA/HiCOLA/Output/logs/%x_%j.out
#SBATCH -e /home/biel/dev-v2/Hi-COLA/HiCOLA/Output/logs/%x_%j.err
# -n: MPI ranks, each evaluates POPSIZE/ranks points per iteration (POPSIZE a multiple of it);
# --cpus-per-task: threads per rank (HODDIES and pycorr), the mocks depend on it
#SBATCH -N 1
#SBATCH -n 18
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=3G
#SBATCH --time=4:00:00
#SBATCH --partition=academic
# Fit of the HOD parameters on the halos of one snapshot of a Hi-COLA run to measured clustering
# (wp and/or xi_ell), with HODDIES mocks and a stochopy minimizer (scripts/run_hod_fit.py). Uses the
# halos of jobs/halo_finder.sh (halos200c_hoddies_dev-v2.npy) and, for SATELLITES="particles", its
# member particles (halos200c_particles.npy, PARTICLE_FRACTION > 0 in the halo finder).
# Results in $HICOLA_DATA_ROOT/hicola/<target>/<model>/<realization>/HiCOLA_<sim>_z<z>_HOD-dev-v2_fit_<name>/,
# named after what is fitted: <tracer>_<HOD model>_<n free>p_<ptcl|nfw>[_<statistics>], the
# statistics only when not all of wp + xi_0 + xi_2, and _r<min>-<max> when the fit range is
# restricted, e.g. LRG_SHOD_6p_ptcl, LRG_SHOD_6p_nfw_wp, LRG_SHOD_6p_ptcl_r0.1-10. A new fit of the
# same thing replaces the earlier results. Shown by Notebooks/HOD_fit_results.ipynb.
#
# Edit the parameters below, then submit it (from any folder):
#   sbatch ~/dev-v2/Hi-COLA/jobs/hod_fit.sh          (or, from the Hi-COLA folder: sbatch jobs/hod_fit.sh)

# ---------------------------- user parameters ----------------------------
# Run file of the simulation, relative to the Hi-COLA folder (or absolute), and the snapshot
# (its halos, and particles for SATELLITES="particles", must exist: jobs/halo_finder.sh)
RUN_FILE="HiCOLA/Parameter_files/Backend_parameter_files/Abacus_Small/LCDM_s1234.lua"
REDSHIFT="0.500"

# HOD model and the values of the fixed parameters: HODDIES parameter file (relative to the
# Hi-COLA folder); its fit settings (fit_param) are not used, they are set here.
HOD_PARAM_FILE="Notebooks/HOD_parameter_file_optimizer.yaml"
TRACER="LRG"
# Free parameters and their (flat) priors, name:min:max
PRIORS="M_0:12.5:13.5 M_1:13:14.5 alpha:0.5:1.5 f_sigv:0.5:1.5 log_Mcent:12.4:13.5 sigma_M:0.05:1"
# Satellites: "particles" = a random member particle of the halo, velocity
# v_halo + f_sigv (v_particle - v_halo) (the NFW profile for halos with too few stored particles);
# "nfw" = HODDIES' NFW profile (concentration c, radius Rh), velocity dispersion Vrms x f_sigv
SATELLITES="particles"

# Data: pycorr estimators, relative to $HICOLA_DATA_ROOT (or absolute)
DATA_RPPI="DESI_Y1_LRG/allcounts_LRG_GCcomb_0.4_0.6_default_angular_bitwise_log_njack128_nran4_split20_rppi.npy"
DATA_SMU="DESI_Y1_LRG/allcounts_LRG_GCcomb_0.4_0.6_default_angular_bitwise_log_njack128_nran4_split20_smu.npy"
STATISTICS="wp xi_ells"   # fitted statistics: "wp xi_ells" (all), "wp" or "xi_ells"; data vector: wp, then xi_ell per ell
# wp(rp): rebinning of the rp bins, rp range (bins entirely inside), pi_max [Mpc/h]
RP_REBIN="2"
RP_RANGE="0.04 32"        # with the files above: 17 bins, 0.046 to 31.6 Mpc/h
PIMAX="40"
# xi_ell(s): multipoles, rebinning of the s bins, s range [Mpc/h]
ELLS="0 2"
S_REBIN="1"
S_RANGE="0.2 32"          # with the files above: 26 bins, 0.215 to 31.6 Mpc/h
# Covariance = (jackknife errors of the data) x (a correlation matrix, e.g. from mocks; relative to
# $HICOLA_DATA_ROOT or absolute), or "none" for a diagonal covariance. CORRELATION_STATISTICS: what
# the matrix covers, with the bins above (wp, then xi_ell per ell); the blocks of the fitted
# statistics are used. The wp and xi_ell settings above are needed for all of them, fitted or not.
CORRELATION="DESI_Y1_LRG/M_corr_lrg_z0.5.txt"
CORRELATION_STATISTICS="wp xi_ells"
# Fit range [Mpc/h], the same for rp and s: only the bins (of the ranges above) entirely between
# FIT_R_MIN and FIT_R_MAX are fitted, e.g. "1e-1" and "1e1"; "None" = no limit on that side. With a
# limit the results are named ..._r<min>-<max> (e.g. _r0.1-10, _rmin-10), without one as usual.
FIT_R_MIN="1.0"
FIT_R_MAX="None"

# Model: line of sight of the redshift-space distortions (box axis); HOD seed of every mock of the
# fit; extra seeds of the best-fit model (for its mock-to-mock scatter)
LOS="z"
MOCK_SEED="1"
MODEL_SEEDS="2 3 4 5 6 7 8 9 10 11"
# HOD modelled at the end for comparison with the best fit (all the tracer's values of this
# HODDIES parameter file, same seeds), relative to the Hi-COLA folder; or "none"
REFERENCE_PARAM_FILE="Notebooks/HOD_parameter_file_z0.500_reference.yaml"

# Minimizer (stochopy): method (cmaes, de, pso, cpso, vdcma, na), iterations, population per
# iteration (evaluated in parallel by the MPI ranks), tolerance on the parameters; starting point
# "random" (drawn in the priors) or "param_file" (the values of HOD_PARAM_FILE), and the seed of
# the start and of the minimizer
METHOD="cmaes"
MAXITER="500"
POPSIZE="72"
XTOL="1e-2"
START="random"
START_SEED="1"
# -------------------------------------------------------------------------
# Box size and cosmology (for H(z) of the redshift-space distortions) are taken from the run file.

set -euo pipefail
# The Hi-COLA folder this script belongs to. SLURM runs a copy of the script, so ask SLURM
# where the original is; outside SLURM, use this file's location.
if [ -n "${SLURM_JOB_ID:-}" ]; then
    SCRIPT=$(scontrol show job "$SLURM_JOB_ID" 2>/dev/null | sed -n 's/^ *Command=\([^ ]*\).*/\1/p')
    JOB_LOG=$(scontrol show job "$SLURM_JOB_ID" 2>/dev/null | sed -n 's/^ *StdOut=\(.*\)$/\1/p')
else
    SCRIPT=${BASH_SOURCE[0]}
    JOB_LOG=""
fi
HICOLA_ROOT=${HICOLA_ROOT:-$(cd "$(dirname "${SCRIPT:-$SLURM_SUBMIT_DIR/jobs/x}")/.." && pwd)}
if [ ! -f "$HICOLA_ROOT/jobs/hicola_env.sh" ]; then
    echo "Cannot find the Hi-COLA folder (looked in $HICOLA_ROOT)" >&2
    exit 1
fi
source "$HICOLA_ROOT/jobs/hicola_env.sh"
cd "$HICOLA_ROOT"
RUN_FILE=$(readlink -f "$RUN_FILE")
HOD_PARAM_FILE=$(readlink -f "$HOD_PARAM_FILE")
if [ "$REFERENCE_PARAM_FILE" != none ]; then REFERENCE_PARAM_FILE=$(readlink -f "$REFERENCE_PARAM_FILE"); fi
data_path() { case "$1" in /*|none) echo "$1" ;; *) echo "$HICOLA_DATA_ROOT/$1" ;; esac; }
DATA_RPPI=$(data_path "$DATA_RPPI")
DATA_SMU=$(data_path "$DATA_SMU")
CORRELATION=$(data_path "$CORRELATION")

hicola_load_run "$RUN_FILE"
if [ "$TARGET" = standalone ]; then
    echo "The HOD fit needs a run matched to a reference simulation, not $RUN_LABEL" >&2
    exit 1
fi
NTHREADS="${SLURM_CPUS_PER_TASK:-1}"
NRANKS="${SLURM_NTASKS:-1}"
if (( POPSIZE % NRANKS != 0 )); then
    echo "Warning: POPSIZE $POPSIZE is not a multiple of the $NRANKS MPI ranks: some ranks wait at each iteration"
fi
RUN=(--target "$TARGET" --model "$MODEL" --realization "$REALIZATION" --redshift "$REDSHIFT" --boxsize "$BOXSIZE")
COSMO_ARGS=(--h "$COSMO_H" --omega-cdm "$COSMO_OMEGA_CDM" --omega-b "$COSMO_OMEGA_B" --omega-mnu "$COSMO_OMEGA_MNU"
           --tcmb "$COSMO_TCMB" --neff "$COSMO_NEFF" --n-s "$COSMO_NS")
read -r -a PRIOR_LIST <<< "$PRIORS"
HOD=(--hod-param-file "$HOD_PARAM_FILE" --tracer "$TRACER" --priors "${PRIOR_LIST[@]}" --nthreads "$NTHREADS")
read -r -a STATISTIC_LIST <<< "$STATISTICS"
read -r -a CORRELATION_STATISTIC_LIST <<< "$CORRELATION_STATISTICS"
read -r -a RP_LIST <<< "$RP_RANGE"
read -r -a S_LIST <<< "$S_RANGE"
read -r -a ELL_LIST <<< "$ELLS"
DATA=(--statistics "${STATISTIC_LIST[@]}" --data-rppi "$DATA_RPPI" --data-smu "$DATA_SMU" --correlation "$CORRELATION" --correlation-statistics "${CORRELATION_STATISTIC_LIST[@]}"
      --pimax "$PIMAX" --rp-rebin "$RP_REBIN" --rp-range "${RP_LIST[@]}"
      --ells "${ELL_LIST[@]}" --s-rebin "$S_REBIN" --s-range "${S_LIST[@]}" --fit-r-min "$FIT_R_MIN" --fit-r-max "$FIT_R_MAX" --los "$LOS"
      --satellites "$SATELLITES")
read -r -a MODEL_SEED_LIST <<< "$MODEL_SEEDS"
FIT=(--method "$METHOD" --maxiter "$MAXITER" --popsize "$POPSIZE" --xtol "$XTOL" --start "$START"
     --start-seed "$START_SEED" --mock-seed "$MOCK_SEED" --model-seeds "${MODEL_SEED_LIST[@]}" --reference-param-file "$REFERENCE_PARAM_FILE"
     --job-log "$JOB_LOG")

hicola_print_env
hicola_print_run "$RUN_FILE"
echo "Snapshot z:     $REDSHIFT"
echo "HOD:            $TRACER, $HOD_PARAM_FILE"
echo "Priors:         $PRIORS"
echo "Satellites:     $SATELLITES"
echo "Data:           fitted: $STATISTICS; wp: rp in [$RP_RANGE] Mpc/h, rebin $RP_REBIN, pi_max $PIMAX; xi_ell: ells $ELLS, s in [$S_RANGE] Mpc/h, rebin $S_REBIN"
echo "                $DATA_RPPI"
echo "                $DATA_SMU"
echo "Fit range:      $FIT_R_MIN to $FIT_R_MAX Mpc/h (rp and s)"
echo "Covariance:     jackknife errors x $CORRELATION (covering $CORRELATION_STATISTICS)"
echo "Model:          RSD along $LOS, mock seed $MOCK_SEED, extra seeds of the best fit $MODEL_SEEDS"
echo "Reference HOD:  $REFERENCE_PARAM_FILE"
echo "Minimizer:      $METHOD, maxiter $MAXITER, popsize $POPSIZE, xtol $XTOL, start $START (seed $START_SEED)"
echo "Resources:      ${SLURM_JOB_NUM_NODES:-1} node(s), $NRANKS MPI ranks x $NTHREADS threads"
echo "================================================================"

# If this job is submitted from inside an interactive allocation, inherited
# SLURM_MEM_* variables can conflict with this script's --mem-per-cpu request.
unset SLURM_MEM_PER_NODE SLURM_MEM_PER_GPU SLURM_CPU_BIND SLURM_CPU_BIND_LIST SLURM_CPU_BIND_TYPE

export NUMBA_NUM_THREADS="$NTHREADS"
time srun --mpi=pmix --cpu-bind=cores -n "$NRANKS" -c "$NTHREADS" \
    python "$HICOLA_ROOT/scripts/run_hod_fit.py" "${RUN[@]}" "${COSMO_ARGS[@]}" "${HOD[@]}" "${DATA[@]}" "${FIT[@]}"
