#!/usr/bin/env bash
#
# Build the protein panel with the label check (task 64): select, run an
# all-vs-all blastp per origin in the pinned BLAST image, exclude
# representatives that look like another gene, re-select, and repeat
# until a pass flags nothing. Then write the final panel.
#
# Usage:
#   scripts/refdata/build_proteins_checked.sh --work /tmp/proteins-work \
#       -- --gene-sets assets/organelle_gene_sets.json \
#          --cds-jsonl plant_pt=pt_cds.jsonl ... \
#          --out refs/v2026.10_1/proteins --refseq-release 237 ...
#
# Everything after `--` is passed to `build_proteins.py full`.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BLAST_IMAGE="quay.io/biocontainers/blast:2.17.0--h66d330f_0"
MAX_PASSES=3
BLAST_EVALUE="1e-3"
BLAST_MAX_TARGETS=500
BLAST_THREADS=8

WORK=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --work) WORK="$2"; shift 2 ;;
        --) shift; break ;;
        *) echo "ERROR: unknown arg: $1" >&2; exit 1 ;;
    esac
done
[[ -n "$WORK" ]] || { echo "ERROR: --work required" >&2; exit 1; }
[[ $# -gt 0 ]]   || { echo "ERROR: no build_proteins.py args" >&2; exit 1; }

mkdir -p "$WORK"
WORK="$(cd "$WORK" && pwd)"
EXCLUSIONS="${WORK}/exclusions.json"
SEL="${WORK}/sel"
rm -f "$EXCLUSIONS"

for pass in $(seq 1 "$MAX_PASSES"); do
    echo "[label-check] pass ${pass}: select" >&2
    exclude_args=()
    [[ -f "$EXCLUSIONS" ]] && exclude_args=(--exclude "$EXCLUSIONS")
    python3 "${SCRIPT_DIR}/build_proteins.py" full "$@" \
        "${exclude_args[@]}" --emit-selection "$SEL"

    for faa in "$SEL"/*.faa; do
        origin="$(basename "$faa" .faa)"
        echo "[label-check] pass ${pass}: blastp ${origin}" >&2
        docker run --rm --user "$(id -u):$(id -g)" \
            -v "${SEL}:/sel" "$BLAST_IMAGE" bash -c "
            makeblastdb -dbtype prot -in /sel/${origin}.faa \
                -out /sel/${origin}.db -logfile /dev/null >/dev/null &&
            blastp -query /sel/${origin}.faa -db /sel/${origin}.db \
                -outfmt '6 qseqid sseqid bitscore' -evalue ${BLAST_EVALUE} \
                -max_target_seqs ${BLAST_MAX_TARGETS} \
                -num_threads ${BLAST_THREADS} \
                -out /sel/${origin}.hits.tsv"
    done

    status=0
    python3 "${SCRIPT_DIR}/build_proteins.py" label_check \
        --selection "$SEL" --exclusions "$EXCLUSIONS" \
        --pass-number "$pass" --max-passes "$MAX_PASSES" || status=$?
    [[ $status -eq 0 ]] && break
    [[ $status -eq 2 ]] || exit "$status"
done

echo "[label-check] clean; writing final panel" >&2
python3 "${SCRIPT_DIR}/build_proteins.py" full "$@" --exclude "$EXCLUSIONS"
