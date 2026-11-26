#!/bin/bash
#
# Copy the small, long-lived analysis products from scratch to CFS, so the
# scratch purge only ever removes data that can be regenerated.
#
# Usage:
#   ./sync_products_to_cfs.sh [SRC] [DST]
# Defaults:
#   SRC = /pscratch/sd/b/bielthi/ndgp_cola_validation/products   (Thibaud's scratch)
#   DST = /global/cfs/cdirs/desi/users/$USER/ndgp_cola_validation/products
#         (change DST if the products should live in a shared DESI project folder)
#
# Copied:   power spectra, halo catalogues, matched-halo tables, HOD mocks, logs,
#           parameter files - everything under SRC except particle data.
# Skipped:  particle snapshots and IC grids (regenerate or re-fetch them instead).
#
# Suggested layout on Thibaud's scratch (/pscratch/sd/b/bielthi):
#   abacus_highbase/                 Abacus slice (retrieve_abacus_z05.sh)
#   ndgp_cola_validation/ics/        IC fields (e.g. delta_..._N1728.f32)
#   ndgp_cola_validation/runs/       COLA run directories with snapshots
#   ndgp_cola_validation/products/   P(k), halo catalogues, tables, mocks
# Anyone with read access can run this script; DST must be writable by the user who runs it.

set -euo pipefail

SCRATCH_ROOT=/pscratch/sd/b/bielthi
SRC=${1:-${SCRATCH_ROOT}/ndgp_cola_validation/products}
DST=${2:-/global/cfs/cdirs/desi/users/${USER}/ndgp_cola_validation/products}

if [[ ! -d "$SRC" ]]; then
    echo "Source directory $SRC does not exist" >&2
    exit 1
fi
mkdir -p "$DST"

rsync -av --partial \
      --chmod=Dug+rwx,Do+rx,Fug+rw,Fo+r \
      --exclude='*.f32' --exclude='*.pack9.asdf' --exclude='snap*' \
      --exclude='*.gadget*' --exclude='gadget*' --exclude='*particles*' \
      "$SRC/" "$DST/"

echo
echo "Synced $(du -sh "$SRC" | cut -f1) on scratch -> $(du -sh "$DST" | cut -f1) on CFS"
echo "Destination: $DST"
