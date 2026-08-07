#!/usr/bin/env bash
#
# Build refs/<version>/annotate/mitos/refseq89m/ — MITOS2 reference data.
# See tasks/completed/31_annotate_merge.md §2.
#
# Downloads the metazoan MITOS2 reference set ("Reference data for
# MITOS2", Zenodo 10.5281/zenodo.4284483, refseq89m.tar.bz2) and
# extracts it into the bundle. Not refseq89f (fungal) or refseq89o
# (other) — animal_mt is the only target with a MITOS2 branch
# (CONSTITUTION constraint 1).
#
# Usage:
#   scripts/refdata/build_annotate.sh \
#       --out     refs/v2026.09/annotate \
#       --scratch /tmp/refdata-annotate

set -euo pipefail

ZENODO_URL="https://zenodo.org/records/4284483/files/refseq89m.tar.bz2?download=1"
EXPECTED_MD5="bb6325e27612e61a6995e88e8ffcecf8"
EXPECTED_BYTES=19629319

OUT_DIR=""
SCRATCH=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --out)      OUT_DIR="$2";  shift 2 ;;
        --scratch)  SCRATCH="$2";  shift 2 ;;
        *) echo "ERROR: unknown arg: $1" >&2; exit 1 ;;
    esac
done

[[ -n "$OUT_DIR" ]] || { echo "ERROR: --out required" >&2; exit 1; }
[[ -n "$SCRATCH" ]] || SCRATCH="${OUT_DIR}.scratch"

mkdir -p "$SCRATCH"
STAGING="${OUT_DIR}.staging"
rm -rf "$STAGING"
mkdir -p "$STAGING/mitos"

# ── 1. Download + verify ──────────────────────────────────────────────────
ARCHIVE="${SCRATCH}/refseq89m.tar.bz2"
if [[ -f "$ARCHIVE" ]]; then
    echo "[annotate] refseq89m.tar.bz2 (cached)"
else
    echo "[annotate] downloading refseq89m.tar.bz2 ..."
    curl -fsSL "$ZENODO_URL" -o "$ARCHIVE"
fi

ACTUAL_BYTES=$(stat -c%s "$ARCHIVE")
ACTUAL_MD5=$(md5sum "$ARCHIVE" | cut -d' ' -f1)

if [[ "$ACTUAL_BYTES" != "$EXPECTED_BYTES" ]]; then
    echo "ERROR: refseq89m.tar.bz2 size mismatch:" \
        "got ${ACTUAL_BYTES}, expected ${EXPECTED_BYTES}" >&2
    exit 1
fi
if [[ "$ACTUAL_MD5" != "$EXPECTED_MD5" ]]; then
    echo "ERROR: refseq89m.tar.bz2 md5 mismatch:" \
        "got ${ACTUAL_MD5}, expected ${EXPECTED_MD5}" >&2
    exit 1
fi
echo "[annotate] md5 verified: ${ACTUAL_MD5}"

# ── 2. Extract ─────────────────────────────────────────────────────────────
echo "[annotate] extracting refseq89m/ ..."
tar -xjf "$ARCHIVE" -C "${STAGING}/mitos"

# ── 3. Provenance sidecar ─────────────────────────────────────────────────
DOWNLOAD_DATE=$(date -u +%Y-%m-%d)
SOURCE_SHA256=$(sha256sum "$ARCHIVE" | cut -d' ' -f1)

cat > "${STAGING}/provenance.json" <<EOF
{
  "source": "Reference data for MITOS2",
  "source_url": "https://doi.org/10.5281/zenodo.4284483",
  "archive": "refseq89m.tar.bz2",
  "release": "refseq89",
  "kingdom": "metazoan",
  "downloaded_at": "${DOWNLOAD_DATE}",
  "archive_bytes": ${EXPECTED_BYTES},
  "archive_md5": "${EXPECTED_MD5}",
  "archive_sha256": "${SOURCE_SHA256}"
}
EOF

# ── 4. Atomic promote ──────────────────────────────────────────────────────
rm -rf "$OUT_DIR"
mkdir -p "$(dirname "$OUT_DIR")"
mv "$STAGING" "$OUT_DIR"

echo "[annotate] done: ${OUT_DIR}"
ls -lh "${OUT_DIR}/mitos/refseq89m" | head
