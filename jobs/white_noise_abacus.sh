#!/bin/bash
#SBATCH --job-name=white_noise_abacus
#SBATCH -o /home/biel/Hi-COLA/HiCOLA/Output/logs/%x_%j.out
#SBATCH -e /home/biel/Hi-COLA/HiCOLA/Output/logs/%x_%j.err
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=72
#SBATCH --mem-per-cpu=3G
#SBATCH --time=2:00:00
#SBATCH --partition=academic
# White noise W(x) of an AbacusSummit simulation for Hi-COLA's ic_random_field_type = "read_white_noise",
# regenerated with zeldovich-PLT (the AbacusSummit IC code) from the simulation's seed (PIPELINE.md, section 6):
#   A  N_LOW^3,  Abacus's own parameters and P(k)      -> compared with Abacus's stored ic_dens_N<N_LOW>.asdf
#   B  N_LOW^3,  flat P(k) (white noise), same seed    -> same phases as A
#   C  N_HIGH^3, flat P(k), same seed                  -> same modes as B (resolution invariance);
#      turned into $OUTPUT_DIR/white_noise_N<N_HIGH>.f32 by scripts/make_white_noise.py
# Checks: scripts/validate_abacus_white_noise.py (reference values in its docstring).
# Needs zeldovich-PLT built in its own conda env (PIPELINE.md, section 6), and in $HICOLA_DATA_ROOT the
# Abacus IC density and CLASS_power (hicola_data/README.md).
# Memory: about 3 x 4 N_HIGH^3 bytes for zeldovich-PLT and ~10 x 4 N_HIGH^3 for the conversion (N_HIGH = 1024: ~45 GB).
#
# Edit the parameters below, then submit it (from any folder):
#   sbatch ~/Hi-COLA/jobs/white_noise_abacus.sh

# ---------------------------- user parameters ----------------------------
SIM_NAME="AbacusSummit_small_c000_ph3000"   # reference simulation (the IC header values below are its own)
IC_DIR="reference/Abacus_Small/LCDM/Box1/raw/ic"   # ic_dens_N<N_LOW>.asdf and CLASS_power, relative to $HICOLA_DATA_ROOT
N_LOW=576               # grid of Abacus's stored IC density (checks A, AB)
N_HIGH=1024             # grid of the white noise for Hi-COLA (= ic_nmesh of the run files)
# IC header of SIM_NAME (ic_dens_N576.asdf, header), passed to zeldovich-PLT unchanged
BOXSIZE=500.0
ZD_SEED=42321
INITIAL_REDSHIFT=99.0
OMEGA_M=0.315192
N_S=0.9649
ZD_PK_NORM=8.0
ZD_PK_SIGMA_RATIO=0.02113950155394198
ZD_F_CLUSTER=0.9954948095129318
ZD_PLT_TARGET_Z=12.0
ZELDOVICH_DIR="$HOME/zeldovich-PLT"         # clone of github.com/abacusorg/zeldovich-PLT, built in build/
ZELDOVICH_ENV="zeldovich"                   # conda env it was built in
WORK_DIR="/scratch/biel/Hi-COLA_IC/Abacus_Small/work"      # zeldovich-PLT runs A, B, C (large, deletable)
OUTPUT_DIR="/scratch/biel/Hi-COLA_IC/Abacus_Small/ph3000"  # white_noise_N<N_HIGH>.f32
THREADS=72
# -------------------------------------------------------------------------

set -euo pipefail
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

ZELDOVICH="$ZELDOVICH_DIR/build/zeldovich"
STORED="$HICOLA_DATA_ROOT/$IC_DIR/ic_dens_N$N_LOW.asdf"
CLASS_POWER="$HICOLA_DATA_ROOT/$IC_DIR/CLASS_power"
OUTPUT="$OUTPUT_DIR/white_noise_N$N_HIGH.f32"
for f in "$ZELDOVICH" "$ZELDOVICH_DIR/build/eigmodes128" "$STORED" "$CLASS_POWER"; do
    [ -e "$f" ] || { echo "Missing: $f" >&2; exit 1; }
done
if [ -e "$OUTPUT" ]; then
    echo "$OUTPUT exists (delete it to remake it)" >&2
    exit 1
fi

hicola_print_env
echo "Simulation:     $SIM_NAME (seed $ZD_SEED, $BOXSIZE Mpc/h, z_init $INITIAL_REDSHIFT)"
echo "zeldovich-PLT:  $ZELDOVICH_DIR (commit $(git -C "$ZELDOVICH_DIR" rev-parse --short HEAD 2>/dev/null || echo unknown)), env $ZELDOVICH_ENV"
echo "Reads:          $STORED, $CLASS_POWER"
echo "Writes:         $WORK_DIR/{A_N${N_LOW}_abacus,B_N${N_LOW}_white,C_N${N_HIGH}_white}/, $OUTPUT"
echo "================================================================"

mkdir -p "$WORK_DIR" "$OUTPUT_DIR"
cd "$WORK_DIR"
printf '%s\n' "1e-10 1" "1e-8 1" "1e-6 1" "1e-4 1" "1e-2 1" "1e0 1" "1e1 1" "1e2 1" "1e3 1" "1e10 1" > white_noise.pow

# zeldovich-PLT parameter file: <run> <N> <abacus|white>
write_par() {
    local run=$1 n=$2 kind=$3
    cat > "$run.par" <<EOF
BoxSize = $BOXSIZE
CPD = 1
ICFormat = "RVZel"
InitialRedshift = $INITIAL_REDSHIFT
ZD_NumBlock = 2
ZD_Seed = $ZD_SEED
ZD_Version = 2
ZD_k_cutoff = 1.0
ZD_CornerModes = 1
ZD_Pk_scale = 1.0
ZD_Pk_smooth = 0.0
ZD_qdensity = 2
ZD_density_filename = "density{:d}"
Omega_M = $OMEGA_M
ZD_n_s = $N_S
NP = $((n * n * n))
InitialConditionsDirectory = "$WORK_DIR/$run"
ZD_f_cluster = $ZD_F_CLUSTER
EOF
    if [ "$kind" = abacus ]; then
        cat >> "$run.par" <<EOF
ZD_Pk_filename = "$CLASS_POWER"
ZD_Pk_norm = $ZD_PK_NORM
ZD_Pk_sigma_ratio = $ZD_PK_SIGMA_RATIO
ZD_qPLT = 0
ZD_qPLT_rescale = 1
ZD_PLT_filename = "$ZELDOVICH_DIR/build/eigmodes128"
ZD_PLT_target_z = $ZD_PLT_TARGET_Z
EOF
    else
        cat >> "$run.par" <<EOF
ZD_Pk_filename = "$WORK_DIR/white_noise.pow"
ZD_Pk_norm = 0.0
ZD_qPLT = 0
ZD_qPLT_rescale = 0
ZD_Pk_sigma_ratio = 1.0
EOF
    fi
}

A=A_N${N_LOW}_abacus; B=B_N${N_LOW}_white; C=C_N${N_HIGH}_white
write_par "$A" "$N_LOW" abacus
write_par "$B" "$N_LOW" white
write_par "$C" "$N_HIGH" white
for run in "$A" "$B" "$C"; do
    rm -rf "$run"; mkdir -p "$run"
    echo "---------------- zeldovich-PLT $run ----------------"
    OMP_NUM_THREADS=$THREADS conda run --no-capture-output -n "$ZELDOVICH_ENV" "$ZELDOVICH" "$run.par" > "$run.log" 2>&1 \
        || { cat "$run.log"; exit 1; }
    tail -n 3 "$run.log"
done

echo "---------------- checks ----------------"
python "$HICOLA_ROOT/scripts/validate_abacus_white_noise.py" \
    --stored "$STORED" --abacus "$A/density$N_LOW" \
    --white-low "$B/density$N_LOW" --white-high "$C/density$N_HIGH" \
    --n-low "$N_LOW" --n-high "$N_HIGH" --threads "$THREADS"

echo "---------------- white noise ----------------"
python "$HICOLA_ROOT/scripts/make_white_noise.py" \
    --input "$C/density$N_HIGH" --nmesh "$N_HIGH" --input-order xyz \
    --output "$OUTPUT.part" --threads "$THREADS"
mv "$OUTPUT.part" "$OUTPUT"
echo "Done: $OUTPUT"
