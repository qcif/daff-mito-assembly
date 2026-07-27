#!/usr/bin/env bash
#
# Build refs/<version>/validate/ — RefSeq organelle BLAST DBs.
# See tasks/3_refdata.md §2.2.
#
# Downloads RefSeq plastid + mitochondrion genomic FASTA files, splits the
# mitochondrion by kingdom (taxonkit on host), and builds three makeblastdb
# databases inside a blast biocontainer.
#
# Tools required on host: curl, python3, taxonkit (at TAXONKIT_BIN)
# Docker image for makeblastdb: BLAST_IMAGE
#
# Usage:
#   scripts/refdata/build_validate.sh \
#       --out         refs/v2026.07/validate \
#       --scratch     /tmp/refdata-validate \
#       [--ncbi-api-key KEY]

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

REFSEQ_FTP="https://ftp.ncbi.nlm.nih.gov/refseq/release"
BLAST_IMAGE="quay.io/biocontainers/blast:2.17.0--h66d330f_0"
TAXONKIT_BIN="/home/cameron/.local/bin/taxonkit"
TAXDUMP_DIR="/home/cameron/.taxonkit/new_v0.20"

OUT_DIR=""
SCRATCH=""
NCBI_API_KEY="${NCBI_API_KEY:-}"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --out)           OUT_DIR="$2";       shift 2 ;;
        --scratch)       SCRATCH="$2";       shift 2 ;;
        --ncbi-api-key)  NCBI_API_KEY="$2";  shift 2 ;;
        *) echo "ERROR: unknown arg: $1" >&2; exit 1 ;;
    esac
done

[[ -n "$OUT_DIR" ]]  || { echo "ERROR: --out required" >&2; exit 1; }
[[ -n "$SCRATCH" ]]  || SCRATCH="${OUT_DIR}.scratch"

mkdir -p "$SCRATCH"
STAGING="${OUT_DIR}.staging"
rm -rf "$STAGING"
mkdir -p "$STAGING"

# ── 1. Pull blast image (idempotent) ──────────────────────────────────────────
echo "[validate] pulling ${BLAST_IMAGE} ..."
docker pull "$BLAST_IMAGE"

# ── 2. RefSeq release number ──────────────────────────────────────────────────
echo "[validate] fetching RefSeq release number ..."
REFSEQ_RELEASE=$(curl -fsSL "${REFSEQ_FTP}/RELEASE_NUMBER")
echo "[validate] RefSeq release: ${REFSEQ_RELEASE}"
echo "$REFSEQ_RELEASE" > "${SCRATCH}/REFSEQ_RELEASE"

DOWNLOAD_DATE=$(date -u +%Y-%m-%d)

# ── helper: download all *.genomic.fna.gz from a RefSeq release subdirectory ──
download_refseq_dir() {
    local subdir="$1"   # e.g. "plastid"
    local dest="$2"     # local directory to save into
    mkdir -p "$dest"
    local base_url="${REFSEQ_FTP}/${subdir}"
    # Scrape the index for all genomic.fna.gz filenames
    local files
    files=$(curl -fsSL "${base_url}/" \
        | grep -o 'href="[^"]*genomic\.fna\.gz"' \
        | sed 's/href="//;s/"//')
    if [[ -z "$files" ]]; then
        echo "ERROR: no genomic.fna.gz found at ${base_url}/" >&2
        exit 1
    fi
    while IFS= read -r fname; do
        local out="${dest}/${fname}"
        if [[ -f "$out" ]]; then
            echo "  → ${fname} (cached)"
        else
            echo "  → ${fname}"
            curl -fsSL "${base_url}/${fname}" -o "$out"
        fi
    done <<< "$files"
}

# ── 3. Download plastid ───────────────────────────────────────────────────────
PLASTID_RAW="${SCRATCH}/plastid_raw"
echo "[validate] downloading plastid genomic FNA ..."
download_refseq_dir "plastid" "$PLASTID_RAW"

echo "[validate] concatenating + decompressing plastid ..."
zcat "${PLASTID_RAW}"/plastid.*.genomic.fna.gz > "${SCRATCH}/refseq_pt.fa"

# ── 4. Download mitochondrion ─────────────────────────────────────────────────
MT_RAW="${SCRATCH}/mt_raw"
echo "[validate] downloading mitochondrion genomic FNA ..."
download_refseq_dir "mitochondrion" "$MT_RAW"

echo "[validate] concatenating + decompressing mitochondrion ..."
zcat "${MT_RAW}"/mitochondrion.*.genomic.fna.gz > "${SCRATCH}/refseq_mt_all.fa"

# ── 5. Kingdom-split mitochondrion ────────────────────────────────────────────
echo "[validate] splitting mitochondrion by kingdom ..."
API_KEY_ARG=""
[[ -n "$NCBI_API_KEY" ]] && API_KEY_ARG="--ncbi-api-key ${NCBI_API_KEY}"

# shellcheck disable=SC2086
python3 "${SCRIPT_DIR}/split_refseq_mt.py" \
    --fasta              "${SCRATCH}/refseq_mt_all.fa" \
    --out-metazoa        "${SCRATCH}/refseq_mt_metazoa.fa" \
    --out-viridiplantae  "${SCRATCH}/refseq_mt_viridiplantae.fa" \
    --taxonkit           "$TAXONKIT_BIN" \
    --taxdump            "$TAXDUMP_DIR" \
    $API_KEY_ARG

# ── 6. Build BLAST databases ──────────────────────────────────────────────────
echo "[validate] building BLAST databases ..."

build_blastdb() {
    local fasta="$1"
    local out_prefix="$2"
    local label="$3"
    echo "  → ${label}"
    docker run --rm \
        --user "$(id -u):$(id -g)" \
        -v "$(realpath "$SCRATCH"):/scratch" \
        -v "$(realpath "$STAGING"):/out" \
        -w /scratch \
        "$BLAST_IMAGE" \
        makeblastdb \
            -in    "/scratch/$(basename "$fasta")" \
            -dbtype nucl \
            -parse_seqids \
            -out   "/out/${out_prefix}" \
            -title "$label"
}

build_blastdb "${SCRATCH}/refseq_pt.fa"            "refseq_pt"            "RefSeq plastid"
build_blastdb "${SCRATCH}/refseq_mt_metazoa.fa"    "refseq_mt_metazoa"    "RefSeq mitochondrion Metazoa"
build_blastdb "${SCRATCH}/refseq_mt_viridiplantae.fa" "refseq_mt_viridiplantae" "RefSeq mitochondrion Viridiplantae"

# ── 7. Write provenance sidecar ───────────────────────────────────────────────
cat > "${STAGING}/provenance.json" <<EOF
{
  "refseq_release": "${REFSEQ_RELEASE}",
  "downloaded_at":  "${DOWNLOAD_DATE}",
  "taxdump_dir":    "${TAXDUMP_DIR}",
  "plastid_source": "${REFSEQ_FTP}/plastid/",
  "mt_source":      "${REFSEQ_FTP}/mitochondrion/"
}
EOF

# ── 8. Atomic promote ─────────────────────────────────────────────────────────
rm -rf "$OUT_DIR"
mkdir -p "$(dirname "$OUT_DIR")"
mv "$STAGING" "$OUT_DIR"

echo "[validate] done: ${OUT_DIR}"
ls -lh "$OUT_DIR"
