#!/usr/bin/env bash
#
# Build refs/<version>/recruit/ from a local GetOrganelleDB checkout.
# See tasks/3_refdata.md §2.1.
#
# Runs everything inside neoformit/daff-wf5-recruit — the same image used
# by the RECRUIT process at pipeline runtime (minimap2/seqtk/samtools),
# so the .mmi indices are built with the same minimap2 version that will
# consume them.
#
# Usage:
#   scripts/refdata/build_recruit.sh \
#       --getorganelle-db reference-material/getorganelle-db/0.0.1 \
#       --out refs/v2026.07/recruit

set -euo pipefail

IMAGE="neoformit/daff-wf5-recruit:latest"

GETORGANELLE_DB=""
OUT_DIR=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --getorganelle-db) GETORGANELLE_DB="$2"; shift 2 ;;
        --out)             OUT_DIR="$2";         shift 2 ;;
        *) echo "ERROR: unknown arg: $1" >&2; exit 1 ;;
    esac
done

[[ -n "$GETORGANELLE_DB" ]] || { echo "ERROR: --getorganelle-db required" >&2; exit 1; }
[[ -n "$OUT_DIR" ]]         || { echo "ERROR: --out required" >&2; exit 1; }

SEED_DIR="${GETORGANELLE_DB}/SeedDatabase"
for f in embplant_pt.fasta embplant_mt.fasta animal_mt.fasta; do
    [[ -f "${SEED_DIR}/${f}" ]] || {
        echo "ERROR: missing ${SEED_DIR}/${f}" >&2
        exit 1
    }
done

# Stage into a sibling dir, atomically rename at the end (task spec §1).
STAGING="${OUT_DIR}.staging"
rm -rf "$STAGING"
mkdir -p "$STAGING"

echo "[recruit] copying seed FASTAs → ${STAGING}/"
cp "${SEED_DIR}/embplant_pt.fasta" "${STAGING}/plant_pt.fa"
cp "${SEED_DIR}/embplant_mt.fasta" "${STAGING}/plant_mt.fa"
cp "${SEED_DIR}/animal_mt.fasta"   "${STAGING}/animal_mt.fa"

echo "[recruit] building minimap2 indices in ${IMAGE}"
docker run --rm \
    --user "$(id -u):$(id -g)" \
    -v "$(realpath "$STAGING"):/work" \
    -w /work \
    "$IMAGE" \
    bash -c '
        set -euo pipefail
        for name in plant_pt plant_mt animal_mt; do
            echo "  → ${name}.mmi"
            minimap2 -d "${name}.mmi" "${name}.fa"
        done
    '

# Atomic promote
rm -rf "$OUT_DIR"
mkdir -p "$(dirname "$OUT_DIR")"
mv "$STAGING" "$OUT_DIR"

echo "[recruit] done: ${OUT_DIR}"
ls -lh "$OUT_DIR"
