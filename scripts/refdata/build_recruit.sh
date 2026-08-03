#!/usr/bin/env bash
#
# Build refs/<version>/recruit/ from a local GetOrganelleDB checkout.
# See tasks/3_refdata.md §2.1 and tasks/completed/28_plastid_masked_mt_panel.md.
#
# Runs everything inside neoformit/daff-wf5-recruit — the same image used
# by the RECRUIT process at pipeline runtime (minimap2/seqtk/samtools),
# so the .mmi indices are built with the same minimap2 version that will
# consume them.
#
# `plant_pt` and `animal_mt` are copied straight from the GetOrganelle
# seed database. `plant_mt` is not: GetOrganelleDB ships exactly one plant
# mitochondrion seed, which is far too thin for the three jobs we ask the
# panel to do (recruitment, the spec §2.1.5 coverage split, and the spec
# §3.7.1 sibling discrimination). It is instead derived from the RefSeq
# Viridiplantae mitochondrion set, with plastid-derived regions (NUPTs)
# masked out — see mask_panel.py for why an unmasked broad panel is
# actively dangerous.
#
# Usage:
#   scripts/refdata/build_recruit.sh \
#       --getorganelle-db reference-material/getorganelle-db/0.0.1 \
#       --refseq-mt       refs/v2026.08/validate/refseq_mt_viridiplantae.fa \
#       --out             refs/v2026.08/recruit

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

IMAGE="neoformit/daff-wf5-recruit:latest"

# Masking knobs — exposed rather than hard-coded (CONSTITUTION rule 18).
# See task 28 §3.2 for what each one buys.
#   asm20 is the divergent genome-to-genome preset; asm10 would mask less
#   and map-ont is wrong for assembly-to-assembly comparison.
MASK_PRESET="asm20"
#   Below 200 bp a hit is a conserved gene fragment shared by both
#   organelles rather than a genuine NUPT insertion.
MASK_MIN_ALIGN_LEN=200
#   A mitogenome measuring >60 % plastid is more likely a misannotated
#   RefSeq record than a real genome: skip it, do not emit mostly-N.
MASK_MAX_MASKED_FRAC=0.60

GETORGANELLE_DB=""
REFSEQ_MT=""
OUT_DIR=""
THREADS="$(nproc)"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --getorganelle-db)  GETORGANELLE_DB="$2";       shift 2 ;;
        --refseq-mt)        REFSEQ_MT="$2";             shift 2 ;;
        --out)              OUT_DIR="$2";               shift 2 ;;
        --mask-preset)      MASK_PRESET="$2";           shift 2 ;;
        --mask-min-align-len)  MASK_MIN_ALIGN_LEN="$2"; shift 2 ;;
        --mask-max-masked-frac) MASK_MAX_MASKED_FRAC="$2"; shift 2 ;;
        --threads)          THREADS="$2";               shift 2 ;;
        *) echo "ERROR: unknown arg: $1" >&2; exit 1 ;;
    esac
done

[[ -n "$GETORGANELLE_DB" ]] || { echo "ERROR: --getorganelle-db required" >&2; exit 1; }
[[ -n "$REFSEQ_MT" ]]       || { echo "ERROR: --refseq-mt required (build_validate.sh emits it)" >&2; exit 1; }
[[ -n "$OUT_DIR" ]]         || { echo "ERROR: --out required" >&2; exit 1; }
[[ -f "$REFSEQ_MT" ]]       || { echo "ERROR: no such file: ${REFSEQ_MT}" >&2; exit 1; }

SEED_DIR="${GETORGANELLE_DB}/SeedDatabase"
for f in embplant_pt.fasta animal_mt.fasta; do
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
cp "${SEED_DIR}/animal_mt.fasta"   "${STAGING}/animal_mt.fa"

# ── plant_mt: align RefSeq mitogenomes against the plastid panel ──────────────
echo "[recruit] aligning RefSeq Viridiplantae mt → plant_pt (${MASK_PRESET})"
cp "$REFSEQ_MT" "${STAGING}/refseq_mt_viridiplantae.fa"

docker run --rm \
    --user "$(id -u):$(id -g)" \
    -v "$(realpath "$STAGING"):/work" \
    -w /work \
    "$IMAGE" \
    bash -c "
        set -euo pipefail
        minimap2 -x ${MASK_PRESET} -t ${THREADS} --secondary=no \
            plant_pt.fa refseq_mt_viridiplantae.fa > mt_vs_pt.paf
    "

echo "[recruit] masking plastid-derived regions out of plant_mt"
python3 "${SCRIPT_DIR}/mask_panel.py" \
    --fasta            "${STAGING}/refseq_mt_viridiplantae.fa" \
    --paf              "${STAGING}/mt_vs_pt.paf" \
    --out              "${STAGING}/plant_mt.fa" \
    --report           "${STAGING}/plant_mt.masking.json" \
    --min-align-len    "$MASK_MIN_ALIGN_LEN" \
    --max-masked-frac  "$MASK_MAX_MASKED_FRAC"

# Intermediates: the source FASTA is retained in validate/, the PAF is
# reproducible from it and the recorded preset.
rm -f "${STAGING}/refseq_mt_viridiplantae.fa" "${STAGING}/mt_vs_pt.paf"

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
