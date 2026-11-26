#!/bin/bash
#SBATCH --job-name=pofk_snapshot
#SBATCH -o /home/biel/Hi-COLA/HiCOLA/Output/logs/%x_%j.out
#SBATCH -e /home/biel/Hi-COLA/HiCOLA/Output/logs/%x_%j.err
#SBATCH -N 4
#SBATCH -n 256
#SBATCH --mem-per-cpu=6G
#SBATCH --cpus-per-task=1
#SBATCH --time=2:00:00
#SBATCH --partition=academic
# P(k) of one snapshot of a Hi-COLA run with other density assignment methods, measured afterwards on the
# Gadget particles with Hi-COLA's own estimator (scripts/pofk_snapshot), then copied to
# $HICOLA_DATA_ROOT/hicola/ as kind PS_<method>.
# The mesh, interlacing and shot-noise subtraction are the run's own (pofk_nmesh, pofk_interlacing,
# pofk_subtract_shotnoise); only the method changes. The run's own method is measured too and compared
# with Hi-COLA's pofk_z<z>.txt (they must agree; the snapshot positions are single precision).
# A run file given on the command line replaces RUN_FILE:  sbatch jobs/pofk_snapshot.sh <run.lua>
#
# Build the program once (and after changes to it or to HiCOLA/FML), in the hicola env:
#   make -C scripts/pofk_snapshot
# Edit the parameters below, then submit it (from any folder):
#   sbatch ~/Hi-COLA/jobs/pofk_snapshot.sh          (or, from the Hi-COLA folder: sbatch jobs/pofk_snapshot.sh)
# To start after the backend job: sbatch --dependency=afterok:<backend job id> jobs/pofk_snapshot.sh <run.lua>
# Memory: two nmesh^3 double-precision grids (interlacing) + the particles, ~0.5 TB for nmesh = 3072.

# ---------------------------- user parameters ----------------------------
# Run file of the simulation to process, relative to the Hi-COLA folder (or absolute)
RUN_FILE="HiCOLA/Parameter_files/Backend_parameter_files/Abacus_Small/LCDM_steps100_zini19_mesh3_s1234.lua"
# Output redshift of the snapshot to process
REDSHIFT="0.500"
# Density assignment methods to measure (NGP CIC TSC PCS PQS)
METHODS="CIC TSC PCS"
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
# A run file given on the command line is relative to where you submitted; RUN_FILE above is
# relative to the Hi-COLA folder
if [ $# -ge 1 ]; then RUN_FILE=$(cd "${SLURM_SUBMIT_DIR:-$PWD}" && readlink -f "$1"); fi
source "$HICOLA_ROOT/jobs/hicola_env.sh"
cd "$HICOLA_ROOT"
RUN_FILE=$(readlink -f "$RUN_FILE")

hicola_load_run "$RUN_FILE"
if [ "$TARGET" = standalone ]; then
    echo "The P(k) job needs a run matched to a reference simulation, not $RUN_LABEL" >&2
    exit 1
fi
PROGRAM="$HICOLA_ROOT/scripts/pofk_snapshot/pofk_snapshot"
if [ ! -x "$PROGRAM" ]; then
    echo "$PROGRAM is not built: make -C scripts/pofk_snapshot (hicola env)" >&2
    exit 1
fi
SNAPSHOT_DIR="$OUTPUT_FOLDER/snapshot_${SIM_NAME}_z${REDSHIFT}"
if ! ls "$SNAPSHOT_DIR/gadget_z${REDSHIFT}".* > /dev/null 2>&1; then
    echo "No Gadget snapshot $SNAPSHOT_DIR/gadget_z${REDSHIFT}.*" >&2
    exit 1
fi
flag() { [ "$1" = true ] && echo 1 || echo 0; }

hicola_print_env
hicola_print_run "$RUN_FILE"
echo "Snapshot:       $SNAPSHOT_DIR/gadget_z${REDSHIFT}.*"
echo "P(k):           methods $METHODS; the run's nmesh $POFK_NMESH, interlacing $POFK_INTERLACING,"
echo "                shot noise subtracted $POFK_SUBTRACT_SHOTNOISE (Hi-COLA's own: $POFK_METHOD)"
echo "Program:        $PROGRAM"
echo "Resources:      ${SLURM_JOB_NUM_NODES:-1} nodes, ${SLURM_NTASKS:-1} MPI tasks"
echo "================================================================"

for METHOD in $METHODS; do
    OUT="$SNAPSHOT_DIR/pofk_${METHOD}_z${REDSHIFT}.txt"
    echo "---- $METHOD"
    time mpirun -np "${SLURM_NTASKS:-1}" "$PROGRAM" "$SNAPSHOT_DIR/gadget_z${REDSHIFT}" "$BOXSIZE" "$((NPART ** 3))" \
        "$POFK_NMESH" "$METHOD" "$(flag "$POFK_INTERLACING")" "$(flag "$POFK_SUBTRACT_SHOTNOISE")" "$OUT"
done

# The run's own method must reproduce Hi-COLA's P(k)
if [[ " $METHODS " == *" $POFK_METHOD "* ]] && [ -f "$SNAPSHOT_DIR/pofk_z${REDSHIFT}.txt" ]; then
    python - "$SNAPSHOT_DIR/pofk_z${REDSHIFT}.txt" "$SNAPSHOT_DIR/pofk_${POFK_METHOD}_z${REDSHIFT}.txt" <<'EOF'
import sys
import numpy as np
k0, p0 = np.loadtxt(sys.argv[1], usecols=(0, 1), unpack=True)
k1, p1 = np.loadtxt(sys.argv[2], usecols=(0, 1), unpack=True)
assert np.allclose(k0, k1, rtol=1e-5), "k bins differ from Hi-COLA's"
r = np.abs(p1 / p0 - 1)
print(f"Check against Hi-COLA's own P(k): max |P/P_HiCOLA - 1| = {r.max():.2e} (median {np.median(r):.2e}) over {len(k0)} bins")
EOF
fi

python "$HICOLA_ROOT/scripts/collect_products.py" "$TARGET" "$MODEL" "$REALIZATION"
