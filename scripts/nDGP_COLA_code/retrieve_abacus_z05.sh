#!/bin/bash
#SBATCH --qos=xfer
#SBATCH --time=12:00:00
#SBATCH --job-name=abacus_z05
#SBATCH --output=%x.%j.log
#SBATCH --licenses=scratch
#
# Retrieve the AbacusSummit_highbase_c000_ph100 full particle slice at z = 0.5
# from NERSC tape (HPSS) into Thibaud's scratch, verify it, and share it.
#
# Default destination: /pscratch/sd/b/bielthi/abacus_highbase (Thibaud, user bielthi)
# Default sharing:     read access for Rafaela (user gsponerr)
#
# Usage (from a Perlmutter login node):
#   Thibaud runs:   sbatch retrieve_abacus_z05.sh
#   anyone else:    DEST=$PSCRATCH/abacus_highbase SHARE_WITH=bielthi sbatch retrieve_abacus_z05.sh
#   no sharing:     SHARE_WITH="" sbatch retrieve_abacus_z05.sh
# Only the owner of a scratch directory can write into it, so with the default
# destination the job must be submitted by Thibaud (the script checks this).
#
# Related public data (already on disk, no retrieval needed):
#   IC density/displacement grids (ic_dens_N576/N1152, ic_disp_N576/N1152):
#   /global/cfs/cdirs/desi/public/cosmosim/AbacusSummit/ic/AbacusSummit_highbase_c000_ph100/
#   Halo catalogues and 10% particle subsamples:
#   /global/cfs/cdirs/desi/public/cosmosim/AbacusSummit/AbacusSummit_highbase_c000_ph100/halos/z0.500/
#
# What it does:
#   0. records quota and the archive listing (index only, no tape read)
#   1. extracts the htar archive straight to disk (no intermediate .tar copy)
#   2. checks every file against the checksums.crc32 files
#   3. gives read access to the users in SHARE_WITH (setfacl)
#   4. writes README_PROVENANCE.txt
# Re-running is safe: the extraction is skipped once it has completed.
#
# If sbatch rejects the xfer QOS or the license line, check NERSC's
# "xfer queue" documentation for the current Perlmutter options.

set -euo pipefail

TAPE_DIR=/nersc/projects/desi/cosmosim/Abacus/AbacusSummit_highbase_c000_ph100/Abacus_AbacusSummit_highbase_c000_ph100_slices
ARCHIVE=Abacus_AbacusSummit_highbase_c000_ph100_z0.500.tar
EXPECTED_BYTES=443766582272
DEST=${DEST:-/pscratch/sd/b/bielthi/abacus_highbase}   # Thibaud's scratch
SHARE_WITH=${SHARE_WITH-gsponerr}         # NERSC usernames, space separated
NPROC_CKSUM=${NPROC_CKSUM:-4}

log() { echo "[$(date '+%F %T')] $*"; }

PARENT=$(dirname "$DEST")
if [[ ! -w "$PARENT" ]]; then
    echo "Cannot write to $PARENT as user ${USER}." >&2
    echo "Either Thibaud (bielthi) submits this job, or set DEST to your own scratch:" >&2
    echo "  DEST=\$PSCRATCH/abacus_highbase SHARE_WITH=bielthi sbatch retrieve_abacus_z05.sh" >&2
    exit 1
fi
mkdir -p "$DEST"
cd "$DEST"
# Lustre striping suited to files of a few hundred MB (NERSC helper, if available)
if command -v stripe_medium >/dev/null 2>&1; then stripe_medium "$DEST" || true; fi

# ---------------------------------------------------------------- 0. checks
log "Quota before retrieval:"
showquota 2>/dev/null || log "(showquota not available)"
log "Archive size on tape: $((EXPECTED_BYTES / 1000000000)) GB; listing contents (index only)"
htar -tvf "${TAPE_DIR}/${ARCHIVE}" > listing_z0.500.txt
log "Archive lists $(grep -c . listing_z0.500.txt) entries -> listing_z0.500.txt"

# ---------------------------------------------------------------- 1. extract
if [[ -f .extraction_complete ]]; then
    log "Extraction already completed earlier, skipping."
else
    log "Extracting ${ARCHIVE} into ${DEST}"
    htar -xvf "${TAPE_DIR}/${ARCHIVE}"
    date '+%F %T' > .extraction_complete
    log "Extraction finished."
fi
log "Size on disk: $(du -sh --exclude='*.txt' . | cut -f1)"

# ---------------------------------------------------------------- 2. checksums
# checksums.crc32 lines are expected in GNU cksum format: CRC SIZE FILENAME
log "Verifying checksums (${NPROC_CKSUM} parallel processes)"
n_dirs=0; n_bad=0
while IFS= read -r ck; do
    d=$(dirname "$ck")
    n_dirs=$((n_dirs + 1))
    computed=$(cd "$d" && awk '{print $3}' checksums.crc32 \
               | xargs -r -P "$NPROC_CKSUM" -n 8 cksum | sort -k3)
    expected=$(sort -k3 "$ck")
    if [[ "$computed" == "$expected" ]]; then
        log "  OK   $d"
    else
        n_bad=$((n_bad + 1))
        log "  FAIL $d  (differences below; also check the checksums.crc32 format)"
        diff <(echo "$expected") <(echo "$computed") | head -20 || true
    fi
done < <(find "$DEST" -name checksums.crc32 | sort)
if (( n_dirs == 0 )); then
    CKSUM_RESULT="no checksums.crc32 files found"
elif (( n_bad == 0 )); then
    CKSUM_RESULT="all ${n_dirs} directories OK"
else
    CKSUM_RESULT="${n_bad} of ${n_dirs} directories FAILED"
fi
log "Checksums: ${CKSUM_RESULT}"

# ---------------------------------------------------------------- 3. share
if [[ -n "$SHARE_WITH" ]]; then
    for u in $SHARE_WITH; do
        log "Granting read access to ${u}"
        setfacl -m "u:${u}:X" "$PARENT"                # traverse only, no listing of the rest
        setfacl -R -m "u:${u}:rX" "$DEST"
        setfacl -R -d -m "u:${u}:rX" "$DEST"           # also for files created later
    done
fi

# ---------------------------------------------------------------- 4. provenance
cat > README_PROVENANCE.txt <<EOF
AbacusSummit_highbase_c000_ph100, full particle slice at z = 0.500
Source (HPSS): ${TAPE_DIR}/${ARCHIVE}
Archive size:  ${EXPECTED_BYTES} bytes
Retrieved:     $(cat .extraction_complete) by ${USER}, Slurm job ${SLURM_JOB_ID:-interactive}
Checksums:     ${CKSUM_RESULT}
Shared with:   ${SHARE_WITH:-nobody}
Note: scratch (${PARENT}) is purged after ~8 weeks without access. Re-run this script to
fetch the slice again; keep derived products on CFS (sync_products_to_cfs.sh).
EOF
log "Done. See ${DEST}/README_PROVENANCE.txt"
