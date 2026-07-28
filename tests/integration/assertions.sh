#!/usr/bin/env bash
# Post-run assertions for integration tests.
# Run after: nextflow run . -profile integration
#
# Assertions are progressive — biology checks are commented out and
# uncommented as downstream stages become real.
#
# Progressive uncomment plan:
#   Stage lands (task)            | Uncomment block
#   ------------------------------|--------------------------------------------------
#   RECRUIT (real minimap2)       | Coverage gate status per sample
#   METAFLYE (real assembler)     | Assembly-length bounds (expected/*/assembly_bounds.json)
#   ANNOTATE (real annotator)     | Barcode loci presence (expected/*/expected_loci.txt)
#   REPORT (real Jinja render)    | run-report.html size > 100 KB, metadata.json schema
#   VALIDATE (real BLAST)         | Taxonomic-identification consistency checks
#
# Each commented block carries a: # TODO(task-N): uncomment when <stage> lands
# Run: grep TODO tests/integration/assertions.sh  — to see outstanding work.
#
# Usage: bash tests/integration/assertions.sh [outdir]

set -euo pipefail

OUTDIR="${1:-tests/integration/output}"
FAILED=0

# --- Per-sample structural checks (always applicable) ---
for sample in INT-ANIMAL-01 INT-PLANT-01-pt INT-PLANT-01-mt; do
    for f in metadata.json report.html; do
        if [[ ! -s "$OUTDIR/$sample/$f" ]]; then
            echo "FAIL: $sample/$f missing or empty"
            FAILED=1
        else
            echo "OK:   $sample/$f"
        fi
    done
done

# --- Run-level checks ---
for f in run_manifest.json run-report.html; do
    if [[ ! -s "$OUTDIR/$f" ]]; then
        echo "FAIL: $f missing or empty"
        FAILED=1
    else
        echo "OK:   $f"
    fi
done

# --- Biology checks (uncomment as stages land) ---
#
# Uncomment the RECRUIT block when COVERAGE_GATE is real (task 10-COV):
# for sample in INT-ANIMAL-01 INT-PLANT-01-pt INT-PLANT-01-mt; do
#     status_file="$OUTDIR/$sample/sample_status.json"
#     if [[ ! -s "$status_file" ]]; then
#         echo "FAIL: $sample/sample_status.json missing"
#         FAILED=1
#         continue
#     fi
#     status=$(jq -r .status "$status_file")
#     if [[ "$status" != "ok" ]]; then
#         echo "FAIL: $sample coverage gate status=$status (expected ok)"
#         FAILED=1
#     fi
# done
#
# Uncomment when METAFLYE is real:
# declare -A SAMPLE_TARGET=(
#     [INT-ANIMAL-01]=animal_mt
#     [INT-PLANT-01-pt]=plant_pt
#     [INT-PLANT-01-mt]=plant_mt
# )
# for sample in INT-ANIMAL-01 INT-PLANT-01-pt INT-PLANT-01-mt; do
#     asm="$OUTDIR/$sample/organelle_assembly.fasta"
#     if [[ ! -s "$asm" ]]; then
#         echo "FAIL: $sample assembly missing or empty"
#         FAILED=1
#         continue
#     fi
#
#     asm_bp=$(grep -v '^>' "$asm" | tr -d '\n' | wc -c)
#     target="${SAMPLE_TARGET[$sample]}"
#     bounds="tests/integration/expected/${target}/assembly_bounds.json"
#     min_bp=$(jq .min_bp "$bounds")
#     max_bp=$(jq .max_bp "$bounds")
#     if (( asm_bp < min_bp || asm_bp > max_bp )); then
#         echo "FAIL: $sample assembly ${asm_bp} bp outside [${min_bp}, ${max_bp}]"
#         FAILED=1
#     else
#         echo "OK:   $sample assembly ${asm_bp} bp in [${min_bp}, ${max_bp}]"
#     fi
# done
#
# Uncomment when MINIPROT_EXTRACT is real:
# for sample in INT-ANIMAL-01 INT-PLANT-01-pt INT-PLANT-01-mt; do
#     barcodes="$OUTDIR/$sample/barcodes.fasta"
#     if [[ ! -s "$barcodes" ]]; then
#         echo "FAIL: $sample/barcodes.fasta missing or empty"
#         FAILED=1
#         continue
#     fi
#     target="${SAMPLE_TARGET[$sample]}"
#     expected_loci="tests/integration/expected/${target}/expected_loci.txt"
#     found=0; missing=0
#     while IFS= read -r locus; do
#         if grep -q ">${locus}" "$barcodes"; then
#             (( found++ )) || true
#         else
#             echo "WARN: $sample missing barcode $locus (not a hard fail)"
#             (( missing++ )) || true
#         fi
#     done < "$expected_loci"
#     if (( found < 3 )); then
#         echo "FAIL: $sample only $found expected loci recovered (need >= 3)"
#         FAILED=1
#     else
#         echo "OK:   $sample $found loci recovered"
#     fi
# done

if [[ "$FAILED" -eq 0 ]]; then
    echo "All assertions passed."
fi
exit "$FAILED"
