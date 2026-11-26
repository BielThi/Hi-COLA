#!/bin/bash
#SBATCH --job-name=white_noise_degrace_dgp
#SBATCH -o /home/biel/Hi-COLA/HiCOLA/Output/logs/%x_%j.out
#SBATCH -e /home/biel/Hi-COLA/HiCOLA/Output/logs/%x_%j.err
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=72
#SBATCH --mem-per-cpu=6G
#SBATCH --time=3:00:00
#SBATCH --partition=academic
# White noise W(x) of DEGRACE-DGP boxes for Hi-COLA's ic_random_field_type = "read_white_noise",
# from GLAM's regenerated ICs (scripts/make_white_noise_degrace_dgp.py), one file per box:
#   $OUTPUT_ROOT/Box<n>/white_noise_N2048.f32 (+ .json)
# Whole node: about 170 GB of memory per box (psi on the 2048^3 lattice + FFTs).
#
# Edit the parameters below, then submit it (from any folder):
#   sbatch ~/Hi-COLA/jobs/white_noise_degrace_dgp.sh

# ---------------------------- user parameters ----------------------------
BOXES="2 3 4 5 6"                                  # boxes with ICs (DEGRACE-DGP boxNNNN/ICs)
IC_ROOT="/pscratch/sd/b/bielthi/DEGRACE-DGP"       # boxNNNN/ICs/ (scratch copy of the drop-off)
LOADER="$HOME/hicola_data/DEGRACE-DGP/DEGRACE-DGP-loader"
PK_FILE="HiCOLA/Parameter_files/Backend_parameter_files/Initial_power_spectra/degrace_dgp_PkTable_z0.000.dat"  # shape of P(k), relative to the Hi-COLA folder
LATTICE_OFFSET="0.125"   # Mpc/h, GLAM's lattice offset on each axis (checked against the particles)
CUTOFF_N="1024"          # GLAM's IC cutoff sphere |k| < CUTOFF_N k_f (checked: the script stops if there is power outside)
NORM_KMAX="5.2052"       # h/Mpc, normalise <|W|^2> = 1 below this k (end of GLAM's own PkTable.dat; the P(k) tail above is extrapolated)
VALIDATE_BOXES="2"       # boxes compared with the CIC density of their particles
VALIDATE_NMESH="495"     # CIC grid of that check (odd: incommensurate with the particle lattice, see the script)
OUTPUT_ROOT="/pscratch/sd/b/bielthi/Hi-COLA_IC/DEGRACE_DGP"
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
export OMP_NUM_THREADS=1

hicola_print_env
echo "Boxes:          $BOXES (validated: $VALIDATE_BOXES, CIC $VALIDATE_NMESH^3)"
echo "Reads:          $IC_ROOT/boxNNNN/ICs, loader $LOADER, P(k) shape $HICOLA_ROOT/$PK_FILE"
echo "Writes:         $OUTPUT_ROOT/Box<n>/white_noise_N2048.f32"
echo "Lattice offset: $LATTICE_OFFSET Mpc/h; IC cutoff |n| < $CUTOFF_N; normalised over k < $NORM_KMAX h/Mpc"
echo "================================================================"

for box in $BOXES; do
    out="$OUTPUT_ROOT/Box$box/white_noise_N2048.f32"
    if [ -e "$out" ]; then
        echo "$out exists: skipping box $box (delete it to remake it)"
        continue
    fi
    mkdir -p "$OUTPUT_ROOT/Box$box"
    nmesh_val=0
    for vb in $VALIDATE_BOXES; do [ "$vb" = "$box" ] && nmesh_val=$VALIDATE_NMESH; done
    echo "---------------- box $box ----------------"
    time python "$HICOLA_ROOT/scripts/make_white_noise_degrace_dgp.py" \
        --ic-dir "$IC_ROOT/$(printf 'box%04d' "$box")/ICs" \
        --loader "$LOADER" \
        --pk "$HICOLA_ROOT/$PK_FILE" \
        --boxsize 1024 --npart-1d 2048 \
        --lattice-offset "$LATTICE_OFFSET" \
        --cutoff-n "$CUTOFF_N" --norm-kmax "$NORM_KMAX" \
        --threads 128 --read-threads 16 \
        --validate-nmesh "$nmesh_val" --validate-procs 28 \
        --output "$out.part"
    mv "$out.part" "$out"
    mv "$out.part.json" "$out.json"
done
