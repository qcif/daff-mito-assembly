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
#   METAFLYE (real assembler)     | Assembly-length bounds (expected/*/assembly_bounds.json) [done]
#   BANDAGE_NG (real renderer)    | PNG magic bytes per sample [done]
#   BIN_TARGET (real C3)          | Contig bp bounds + circularity (expected/*/bin_bounds.json) [done]
#   Plastid canonicalisation (C4) | Canonicalisation branch + isoform files (plant_pt) [done]
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
# Use -e (exists) not -s (non-empty): stub script blocks produce empty placeholder
# files until real tools land. Content checks are in the biology blocks below.
for sample in INT-ANIMAL-01 INT-PLANT-01-pt INT-PLANT-01-mt; do
    for f in metadata.json report.html; do
        if [[ ! -e "$OUTDIR/$sample/$f" ]]; then
            echo "FAIL: $sample/$f missing"
            FAILED=1
        else
            echo "OK:   $sample/$f"
        fi
    done
done

# --- Run-level checks ---
for f in run_manifest.json run-report.html; do
    if [[ ! -e "$OUTDIR/$f" ]]; then
        echo "FAIL: $f missing"
        FAILED=1
    else
        echo "OK:   $f"
    fi
done

# --- Biology checks (uncomment as stages land) ---

# COVERAGE_GATE is real (task 15):
for sample in INT-ANIMAL-01 INT-PLANT-01-pt INT-PLANT-01-mt; do
    status_file="$OUTDIR/$sample/coverage_gate/sample_status.json"
    if [[ ! -s "$status_file" ]]; then
        echo "FAIL: $sample/coverage_gate/sample_status.json missing"
        FAILED=1
        continue
    fi
    status=$(jq -r .status "$status_file")
    if [[ "$status" != "ok" ]]; then
        echo "FAIL: $sample coverage gate status=$status (expected ok)"
        FAILED=1
    else
        echo "OK:   $sample coverage gate status=ok"
    fi
done


# METAFLYE is real (task 16):
declare -A SAMPLE_TARGET=(
    [INT-ANIMAL-01]=animal_mt
    [INT-PLANT-01-pt]=plant_pt
    [INT-PLANT-01-mt]=plant_mt
)
for sample in INT-ANIMAL-01 INT-PLANT-01-pt INT-PLANT-01-mt; do
    asm="$OUTDIR/$sample/assembly/assembly.fasta"
    if [[ ! -s "$asm" ]]; then
        echo "FAIL: $sample assembly missing or empty"
        FAILED=1
        continue
    fi

    asm_bp=$(grep -v '^>' "$asm" | tr -d '\n' | wc -c)
    target="${SAMPLE_TARGET[$sample]}"
    bounds="tests/integration/expected/${target}/assembly_bounds.json"
    min_bp=$(jq .min_bp "$bounds")
    max_bp=$(jq .max_bp "$bounds")
    if (( asm_bp < min_bp || asm_bp > max_bp )); then
        echo "FAIL: $sample assembly ${asm_bp} bp outside [${min_bp}, ${max_bp}]"
        FAILED=1
    else
        echo "OK:   $sample assembly ${asm_bp} bp in [${min_bp}, ${max_bp}]"
    fi
done

# BANDAGE_NG is real (task 17):
for sample in INT-ANIMAL-01 INT-PLANT-01-pt INT-PLANT-01-mt; do
    png="$OUTDIR/$sample/assembly/${sample}.graph.png"
    if [[ ! -s "$png" ]]; then
        echo "FAIL: $sample graph PNG missing or empty"
        FAILED=1
        continue
    fi
    # PNG magic bytes: 89 50 4E 47 0D 0A 1A 0A
    magic=$(head -c 8 "$png" | xxd -p)
    if [[ "$magic" != "89504e470d0a1a0a" ]]; then
        echo "FAIL: $sample graph PNG magic bytes bad ($magic)"
        FAILED=1
    else
        echo "OK:   $sample graph PNG ($(wc -c < "$png") bytes)"
    fi
done


# BIN_TARGET is real (task 18):
for sample in INT-ANIMAL-01 INT-PLANT-01-pt INT-PLANT-01-mt; do
    tgt="$OUTDIR/$sample/bin_target/target.fasta"
    meta="$OUTDIR/$sample/bin_target/bin_metadata.json"
    if [[ ! -s "$tgt" ]]; then
        echo "FAIL: $sample target.fasta missing or empty"
        FAILED=1
        continue
    fi
    n_contigs=$(grep -c '^>' "$tgt")
    tgt_bp=$(grep -v '^>' "$tgt" | tr -d '\n' | wc -c)
    target="${SAMPLE_TARGET[$sample]}"
    bounds="tests/integration/expected/${target}/bin_bounds.json"
    min_bp=$(jq .target_min_bp "$bounds")
    max_bp=$(jq .target_max_bp "$bounds")
    max_contigs=$(jq .target_max_contigs "$bounds")
    if (( tgt_bp < min_bp || tgt_bp > max_bp )); then
        echo "FAIL: $sample target ${tgt_bp} bp outside [${min_bp}, ${max_bp}]"
        FAILED=1
    elif (( n_contigs > max_contigs )); then
        echo "FAIL: $sample $n_contigs contigs > max ${max_contigs}"
        FAILED=1
    else
        echo "OK:   $sample target ${tgt_bp} bp / ${n_contigs} contigs"
    fi
    # animal_mt: additionally assert circularity.
    if [[ "$target" == "animal_mt" ]]; then
        circular=$(jq -r .circular "$meta")
        if [[ "$circular" != "true" ]]; then
            echo "FAIL: $sample animal_mt target not detected circular"
            FAILED=1
        fi
    fi
done


# Plastid canonicalisation (C4) is real (task 20):
for sample in INT-ANIMAL-01 INT-PLANT-01-pt INT-PLANT-01-mt; do
    meta="$OUTDIR/$sample/bin_target/bin_metadata.json"
    iso_dir="$OUTDIR/$sample/bin_target/plastid_isoforms"

    if [[ "$sample" != "INT-PLANT-01-pt" ]]; then
        if [[ -d "$iso_dir" ]]; then
            echo "FAIL: $sample unexpected plastid_isoforms/ directory"
            FAILED=1
        else
            echo "OK:   $sample no plastid_isoforms/ (not plant_pt)"
        fi
        continue
    fi

    if [[ ! -s "$meta" ]]; then
        echo "FAIL: $sample bin_metadata.json missing"
        FAILED=1
        continue
    fi

    branch=$(jq -r '.plastid_canonicalisation.branch // empty' "$meta")
    case "$branch" in
        canonical)
            tgt="$OUTDIR/$sample/bin_target/target.fasta"
            if [[ -s "$iso_dir/path1.fasta" && -s "$iso_dir/path2.fasta" ]]; then
                echo "OK:   $sample plastid_isoforms/path{1,2}.fasta present"
            else
                echo "FAIL: $sample plastid_isoforms/path1.fasta or path2.fasta missing"
                FAILED=1
            fi
            if cmp -s "$iso_dir/path1.fasta" "$tgt"; then
                echo "OK:   $sample target.fasta == plastid_isoforms/path1.fasta"
            else
                echo "FAIL: $sample target.fasta differs from plastid_isoforms/path1.fasta"
                FAILED=1
            fi
            ;;
        resolved_circle)
            if [[ -d "$iso_dir" ]]; then
                echo "FAIL: $sample resolved_circle but plastid_isoforms/ exists"
                FAILED=1
            else
                echo "OK:   $sample resolved_circle, no plastid_isoforms/"
            fi
            ;;
        non_canonical)
            n=$(jq -r '.plastid_canonicalisation.edge_count // "?"' "$meta")
            echo "WARN: $sample plastid graph non_canonical (edge_count=$n)"
            ;;
        *)
            echo "FAIL: $sample plastid_canonicalisation branch missing or unrecognised ($branch)"
            FAILED=1
            ;;
    esac
done

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
