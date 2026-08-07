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
#   RECRUIT (real minimap2)       | Coverage gate status per sample [done — task 25 makes
#                                 |   the gate sibling-aware and asserts per-sample expected
#                                 |   status + the recruited-pool split (expected/*/coverage_bounds.json)]
#   METAFLYE (real assembler)     | Assembly-length bounds (expected/*/assembly_bounds.json) [done]
#   BANDAGE_NG (real renderer)    | PNG magic bytes per sample [done]
#   BIN_TARGET (real C3)          | Contig bp bounds + circularity (expected/*/bin_bounds.json) [done —
#                                 |   recalibrated by task 23: adds n_target_selected >= 1, no
#                                 |   sibling_organelle emitted, circular_method == flye_circ. Task 25
#                                 |   adds the sibling_carryover cross-check and drops the
#                                 |   INT-PLANT-01-mt plastid-contig check — that sample no longer
#                                 |   reaches BIN_TARGET (see ASSEMBLING_SAMPLES below)]
#   Plastid canonicalisation (C4) | Canonicalisation branch + isoform files (plant_pt) [done —
#                                 |   task 24 adds substitution_applied && n_target_selected >= 1
#                                 |   and target_source == c4_plastid_path1 on the canonical branch]
#   MINIPROT_CDS + EXTRACT_BARCODES | Barcode loci presence + cds.gff coherence
#                                 |   (expected/*/expected_loci.txt) [done — task 30]
#   ANNOTATE (real annotator)     | annotation_summary.json status/counts/crosscheck
#                                 |   (expected/animal_mt/annotation_bounds.json) [done — task 31]
#   REPORT (real Jinja render)    | run-report.html size > 100 KB, metadata.json schema
#   VALIDATE (real BLAST)         | Taxonomic-identification consistency checks [done — task 21]
#
# Each commented block carries a: # TODO(task-N): uncomment when <stage> lands
# Run: grep TODO tests/integration/assertions.sh  — to see outstanding work.
#
# Usage: bash tests/integration/assertions.sh [outdir]
#
# NOTE: publishDir never deletes. If a sample stopped producing an
# output (e.g. INT-PLANT-01-mt no longer assembles — see
# ASSEMBLING_SAMPLES below), a stale file from an earlier run will still
# be sitting in the outdir and the absence checks below will report on
# it rather than on this run. Clear the outdir and re-run with -resume
# (which republishes from the work cache) before trusting a result.

set -euo pipefail

OUTDIR="${1:-tests/integration/output}"
FAILED=0

SAMPLES=(INT-ANIMAL-01 INT-PLANT-01-pt INT-PLANT-01-mt)

declare -A SAMPLE_TARGET=(
    [INT-ANIMAL-01]=animal_mt
    [INT-PLANT-01-pt]=plant_pt
    [INT-PLANT-01-mt]=plant_mt
)

# Samples expected to clear the coverage gate and reach assembly.
# INT-PLANT-01-mt soft-fails at C2: task 25 §2 measured its recruited
# pool at 70.5 % plastid, leaving 27.29× mitochondrial depth against a
# 30× MIN. Everything downstream of the gate is therefore asserted over
# this list, not over SAMPLES. The fixture is under-sequenced for its
# declared target and cannot support mitogenome assertions — replacing
# it is tracked in tasks/todo.md.
ASSEMBLING_SAMPLES=(INT-ANIMAL-01 INT-PLANT-01-pt)

# --- Per-sample structural checks (always applicable) ---
# Use -e (exists) not -s (non-empty): stub script blocks produce empty placeholder
# files until real tools land. Content checks are in the biology blocks below.
for sample in "${SAMPLES[@]}"; do
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

# COVERAGE_GATE is real (task 15; made sibling-aware by task 25):
for sample in "${SAMPLES[@]}"; do
    status_file="$OUTDIR/$sample/coverage_gate/sample_status.json"
    if [[ ! -s "$status_file" ]]; then
        echo "FAIL: $sample/coverage_gate/sample_status.json missing"
        FAILED=1
        continue
    fi
    target="${SAMPLE_TARGET[$sample]}"
    bounds="tests/integration/expected/${target}/coverage_bounds.json"

    status=$(jq -r .status "$status_file")
    want_status=$(jq -r .expected_status "$bounds")
    if [[ "$status" != "$want_status" ]]; then
        echo "FAIL: $sample coverage gate status=$status (expected $want_status)"
        FAILED=1
    else
        echo "OK:   $sample coverage gate status=$status"
    fi

    # The estimate must rest on target-assigned bases wherever the
    # target has a sibling organelle panel (spec §2.1.5). A silent
    # fallback to the whole recruited pool is the defect task 25 fixed.
    basis=$(jq -r '.coverage_basis // "missing"' "$status_file")
    want_basis=$(jq -r .coverage_basis "$bounds")
    if [[ "$basis" != "$want_basis" ]]; then
        echo "FAIL: $sample coverage_basis=$basis (expected $want_basis)"
        FAILED=1
    else
        echo "OK:   $sample coverage_basis=$basis"
    fi

    est=$(jq -r .estimated_cov "$status_file")
    min_cov=$(jq -r .min_estimated_cov "$bounds")
    max_cov=$(jq -r .max_estimated_cov "$bounds")
    if awk "BEGIN{exit !($est >= $min_cov && $est <= $max_cov)}"; then
        echo "OK:   $sample estimated_cov=${est}x in [${min_cov}, ${max_cov}]"
    else
        echo "FAIL: $sample estimated_cov=${est}x outside [${min_cov}, ${max_cov}]"
        FAILED=1
    fi

    # Sibling-organelle carry-over in the recruited pool. plant_mt is
    # plastid-dominated and plant_pt is not; asserting both directions
    # catches a regression that disabled the split as well as one that
    # mis-assigned reads.
    want_max_sib=$(jq -r '.max_sibling_organelle_fraction // empty' "$bounds")
    want_min_sib=$(jq -r '.min_sibling_organelle_fraction // empty' "$bounds")
    if [[ -n "$want_max_sib" || -n "$want_min_sib" ]]; then
        sib=$(jq -r '.sibling_organelle_fraction // "missing"' "$status_file")
        if [[ "$sib" == "missing" ]]; then
            echo "FAIL: $sample sibling_organelle_fraction not recorded"
            FAILED=1
        elif [[ -n "$want_max_sib" ]] \
             && awk "BEGIN{exit !($sib > $want_max_sib)}"; then
            echo "FAIL: $sample sibling fraction ${sib} > ${want_max_sib}"
            FAILED=1
        elif [[ -n "$want_min_sib" ]] \
             && awk "BEGIN{exit !($sib < $want_min_sib)}"; then
            echo "FAIL: $sample sibling fraction ${sib} < ${want_min_sib}"
            FAILED=1
        else
            echo "OK:   $sample sibling_organelle_fraction=${sib}"
        fi
    fi
done

# A soft-failed sample must not have run the assembler (spec §2.1.4).
for sample in "${SAMPLES[@]}"; do
    if [[ " ${ASSEMBLING_SAMPLES[*]} " == *" $sample "* ]]; then
        continue
    fi
    if [[ -s "$OUTDIR/$sample/assembly/assembly.fasta" ]]; then
        echo "FAIL: $sample soft-failed the gate but produced an assembly"
        FAILED=1
    else
        echo "OK:   $sample soft-failed, no assembly attempted"
    fi
done


# METAFLYE is real (task 16):
for sample in "${ASSEMBLING_SAMPLES[@]}"; do
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
for sample in "${ASSEMBLING_SAMPLES[@]}"; do
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


# BIN_TARGET is real (task 18; criteria recalibrated by task 23):
for sample in "${ASSEMBLING_SAMPLES[@]}"; do
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
    # C3 must select on its own evidence — the assertion whose absence
    # let plant_pt pass while selecting nothing (task 23 §3).
    n_selected=$(jq -r '.n_target_selected // 0' "$meta")
    if (( n_selected < 1 )); then
        echo "FAIL: $sample n_target_selected=${n_selected} (expected >= 1)"
        FAILED=1
    else
        echo "OK:   $sample n_target_selected=${n_selected}"
    fi

    # No emitted contig may be a sibling organelle (spec §3.7.1).
    n_sibling=$(jq -r '
        [.contigs[]?
         | select(.classification == "sibling_organelle")
         | .contig_id] as $sib
        | [.contigs_selected[]? | select(. as $c | $sib | index($c))]
        | length' "$meta")
    if (( n_sibling > 0 )); then
        echo "FAIL: $sample emitted ${n_sibling} sibling_organelle contig(s)"
        FAILED=1
    else
        echo "OK:   $sample no sibling_organelle contigs emitted"
    fi

    # animal_mt: circularity must come from Flye's circ. column, not
    # from the end-overlap fallback (spec §3.7.4).
    if [[ "$target" == "animal_mt" ]]; then
        circular=$(jq -r '.circular' "$meta")
        circ_method=$(jq -r '.circular_method // "missing"' "$meta")
        if [[ "$circular" != "true" ]]; then
            echo "FAIL: $sample animal_mt target not detected circular"
            FAILED=1
        elif [[ "$circ_method" != "flye_circ" ]]; then
            echo "FAIL: $sample circular_method=${circ_method} (expected flye_circ)"
            FAILED=1
        else
            echo "OK:   $sample circular via ${circ_method}"
        fi
    fi

    # Post-assembly cross-check on the same signal the gate now uses
    # (spec §2.1.5, task 25 §3.1): how much assembled sequence binned as
    # the sibling organelle, and did that trip the report warning.
    # Note: jq's // treats `false` as empty, and sibling_carryover_warning
    # is legitimately false on a clean assembly — test for the key, not
    # for a truthy value.
    sib_frac=$(jq -r '
        if has("sibling_carryover")
           and (.sibling_carryover | has("sibling_organelle_fraction"))
        then .sibling_carryover.sibling_organelle_fraction
        else "missing" end' "$meta")
    sib_warn=$(jq -r '
        if has("sibling_carryover")
           and (.sibling_carryover | has("sibling_carryover_warning"))
        then .sibling_carryover.sibling_carryover_warning
        else "missing" end' "$meta")
    if [[ "$sib_frac" == "missing" || "$sib_warn" == "missing" ]]; then
        echo "FAIL: $sample sibling_carryover not recorded in bin_metadata.json"
        FAILED=1
    else
        echo "OK:   $sample sibling_carryover fraction=${sib_frac} warning=${sib_warn}"
    fi
done


# Plastid canonicalisation (C4) is real (task 20):
for sample in "${ASSEMBLING_SAMPLES[@]}"; do
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
            # Task 24: the substitution is conditional on C3 having
            # selected a plant_pt contig. Assert both together so a
            # confident plastome can never again ship off a 3-edge
            # graph alone (spec §3.6 step 5).
            applied=$(jq -r '
                .plastid_canonicalisation.substitution_applied' "$meta")
            n_selected=$(jq -r '.n_target_selected // 0' "$meta")
            src=$(jq -r '.target_source // "missing"' "$meta")
            if [[ "$applied" == "true" && $n_selected -ge 1 ]]; then
                echo "OK:   $sample substitution_applied with n_target_selected=${n_selected}"
            else
                echo "FAIL: $sample substitution_applied=${applied} n_target_selected=${n_selected}"
                FAILED=1
            fi
            if [[ "$src" == "c4_plastid_path1" ]]; then
                echo "OK:   $sample target_source=${src}"
            else
                echo "FAIL: $sample target_source=${src} (expected c4_plastid_path1)"
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

# BLAST_VALIDATE is real (task 21). Asserted over ASSEMBLING_SAMPLES,
# not SAMPLES: a sample that soft-fails the coverage gate never
# reaches BIN_TARGET and so produces no blast.tsv at all.
declare -A EXPECTED_TITLE_SUBSTR=(
    [INT-ANIMAL-01]="Acyrthosiphon"
    [INT-PLANT-01-pt]="plastid|chloroplast"
)
for sample in "${ASSEMBLING_SAMPLES[@]}"; do
    tsv="$OUTDIR/$sample/blast_validate/${sample}.blast.tsv"
    tgt="$OUTDIR/$sample/bin_target/target.fasta"
    if [[ ! -e "$tsv" ]]; then
        echo "FAIL: $sample blast.tsv missing"
        FAILED=1
        continue
    fi
    # An empty target.fasta (no C3 selection, or a withheld C4
    # substitution) legitimately yields an empty TSV — only demand
    # hits when there was actually a query.
    if [[ ! -s "$tgt" ]]; then
        echo "SKIP: $sample target.fasta empty — no BLAST query expected"
        continue
    fi
    if [[ ! -s "$tsv" ]]; then
        echo "FAIL: $sample blast.tsv empty but target.fasta is non-empty"
        FAILED=1
        continue
    fi
    # Column 8 is stitle. Grep is line-wise; case-insensitive substr match.
    if ! head -n 5 "$tsv" | cut -f8 \
            | grep -qEi "${EXPECTED_TITLE_SUBSTR[$sample]}"; then
        echo "FAIL: $sample top-5 hits do not mention '${EXPECTED_TITLE_SUBSTR[$sample]}'"
        FAILED=1
    else
        echo "OK:   $sample BLAST top hits look on-target"
    fi
done

# MINIPROT_CDS + EXTRACT_BARCODES are real (task 30). Asserted over
# ASSEMBLING_SAMPLES, not SAMPLES — INT-PLANT-01-mt never reaches
# BIN_TARGET (see ASSEMBLING_SAMPLES above), so it has no cds.gff /
# barcodes.fasta at all (task 21 §4's standing consequence).
#
# expected_loci.txt carries a broader gene-completeness list than the
# barcode panel (assets/loci.json) — e.g. it names plant_pt tRNA genes
# that miniprot cannot recover at all. Only the intersection is a
# barcode candidate; even within that intersection, recovery depends
# on real assembly identity against a comprehensive, multi-species
# protein panel and is not guaranteed per-locus (task 29's own go/no-go
# measurement already flagged animal_mt CYTB as borderline under the
# broad panel — task 30 outcomes confirms it misses the identity floor
# on this fixture). So individual misses are WARN, not FAIL; only an
# empty barcodes.fasta (nothing recovered at all) is a hard failure.
for sample in "${ASSEMBLING_SAMPLES[@]}"; do
    cds_gff="$OUTDIR/$sample/cds/${sample}.cds.gff"
    barcodes="$OUTDIR/$sample/barcodes/barcodes.fasta"
    validation_tsv="$OUTDIR/$sample/barcodes/${sample}.validation.tsv"

    if [[ ! -s "$barcodes" ]]; then
        echo "FAIL: $sample barcodes.fasta missing or empty"
        FAILED=1
        continue
    fi

    target="${SAMPLE_TARGET[$sample]}"
    expected_loci="tests/integration/expected/${target}/expected_loci.txt"
    found=0; candidates=0
    while IFS= read -r locus; do
        if ! jq -e --arg g "$locus" \
                "[.${target}[] | ascii_downcase] | index(\$g | ascii_downcase)" \
                assets/loci.json >/dev/null; then
            continue  # not a barcode-panel locus — nothing to recover
        fi
        (( candidates++ )) || true
        if grep -q ">${locus}_" "$barcodes"; then
            (( found++ )) || true
        else
            echo "WARN: $sample missing barcode $locus (not a hard fail — task 30 outcomes)"
        fi
    done < "$expected_loci"
    if (( found > 0 )); then
        echo "OK:   $sample ${found}/${candidates} barcode-panel loci recovered"
    else
        echo "FAIL: $sample recovered 0/${candidates} barcode-panel loci"
        FAILED=1
    fi

    # validation.tsv carries one row per panel locus (found or not).
    if [[ ! -s "$validation_tsv" ]]; then
        echo "FAIL: $sample validation.tsv missing or empty"
        FAILED=1
    else
        n_loci=$(jq -r ".${target} | length" assets/loci.json)
        n_rows=$(( $(wc -l < "$validation_tsv") - 1 ))
        if (( n_rows < n_loci )); then
            echo "FAIL: $sample validation.tsv has ${n_rows} rows, expected >= ${n_loci} panel loci"
            FAILED=1
        else
            echo "OK:   $sample validation.tsv has ${n_rows} rows (>= ${n_loci} panel loci)"
        fi
    fi

    # Coherence invariant (task 30 §3): every barcodes.fasta record
    # traces to a cds.gff feature at identical coordinates — the whole
    # point of the unified miniprot pass, asserted in CI, not just
    # unit tests. Read seqid/start/end from validation.tsv's "pass"
    # rows rather than parsing the FASTA header — both locus and
    # seqid can themselves contain underscores (e.g. "contig_10"),
    # which makes the <locus>_<seqid>_<start>_<end> header ambiguous
    # to split back apart.
    if [[ ! -s "$cds_gff" ]]; then
        echo "FAIL: $sample cds.gff missing or empty"
        FAILED=1
        continue
    fi
    mismatch=0
    while IFS=$'\t' read -r seqid start end; do
        if ! awk -F'\t' -v s="$seqid" -v a="$start" -v b="$end" \
                '$1==s && $3=="CDS" && $4==a && $5==b {found=1}
                 END{exit !found}' "$cds_gff"; then
            echo "FAIL: $sample barcode at ${seqid}:${start}-${end} has no matching cds.gff row"
            mismatch=1
            FAILED=1
        fi
    done < <(awk -F'\t' 'NR>1 && $2=="pass" {print $4"\t"$5"\t"$6}' "$validation_tsv")
    if (( mismatch == 0 )); then
        echo "OK:   $sample all barcodes.fasta records match a cds.gff row"
    fi
done

# ANNOTATE is real (task 31). animal_mt: cds.gff CDS + MITOS2 tRNA/rRNA
# merge, checked against expected/animal_mt/annotation_bounds.json.
# plant_pt: no non-CDS annotator yet (task 32) — CDS-only.
annotation_summary="$OUTDIR/INT-ANIMAL-01/annotation/annotation_summary.json"
annotation_gff="$OUTDIR/INT-ANIMAL-01/annotation/INT-ANIMAL-01.gff"
bounds="tests/integration/expected/animal_mt/annotation_bounds.json"
if [[ ! -s "$annotation_summary" ]]; then
    echo "FAIL: INT-ANIMAL-01 annotation_summary.json missing or empty"
    FAILED=1
else
    status=$(jq -r '.status' "$annotation_summary")
    if [[ "$status" != "ok" ]]; then
        echo "FAIL: INT-ANIMAL-01 annotation status='$status', expected 'ok'"
        FAILED=1
    else
        echo "OK:   INT-ANIMAL-01 annotation status='ok'"
    fi

    cds_source=$(jq -r '.cds_source' "$annotation_summary")
    if [[ "$cds_source" != "miniprot" ]]; then
        echo "FAIL: INT-ANIMAL-01 cds_source='$cds_source', expected 'miniprot'"
        FAILED=1
    fi

    n_cds=$(jq -r '.feature_counts.CDS' "$annotation_summary")
    n_trna=$(jq -r '.feature_counts.tRNA' "$annotation_summary")
    n_rrna=$(jq -r '.feature_counts.rRNA' "$annotation_summary")
    min_cds=$(jq -r '.min_cds' "$bounds")
    min_trna=$(jq -r '.min_trna' "$bounds")
    min_rrna=$(jq -r '.min_rrna' "$bounds")

    for pair in "CDS:$n_cds:$min_cds" "tRNA:$n_trna:$min_trna" "rRNA:$n_rrna:$min_rrna"; do
        IFS=: read -r label n min <<< "$pair"
        if (( n < min )); then
            echo "FAIL: INT-ANIMAL-01 ${label} count ${n} below floor ${min}"
            FAILED=1
        else
            echo "OK:   INT-ANIMAL-01 ${label} count ${n} (>= ${min})"
        fi
    done

    while IFS= read -r gene; do
        if ! grep -qi "${gene}" "$annotation_gff"; then
            echo "FAIL: INT-ANIMAL-01 required gene '$gene' absent from annotation GFF"
            FAILED=1
        else
            echo "OK:   INT-ANIMAL-01 required gene '$gene' present in annotation GFF"
        fi
    done < <(jq -r '.required_genes[]' "$bounds")

    n_conflicts=$(jq -r '.cds_crosscheck.coordinate_conflicts | length' "$annotation_summary")
    max_conflicts=$(jq -r '.max_cds_crosscheck_conflicts' "$bounds")
    if (( n_conflicts > max_conflicts )); then
        echo "FAIL: INT-ANIMAL-01 cds_crosscheck conflicts ${n_conflicts} exceeds max ${max_conflicts}"
        FAILED=1
    else
        echo "OK:   INT-ANIMAL-01 cds_crosscheck conflicts ${n_conflicts} (<= ${max_conflicts})"
    fi

    # Provenance: every CDS row in the merged annotation GFF traces
    # exactly to a cds.gff row (task 31 §1 — no CDS feature ever
    # originates from MITOS2).
    cds_gff="$OUTDIR/INT-ANIMAL-01/cds/INT-ANIMAL-01.cds.gff"
    mismatch=0
    while IFS=$'\t' read -r seqid start end; do
        if ! awk -F'\t' -v s="$seqid" -v a="$start" -v b="$end" \
                '$1==s && $3=="CDS" && $4==a && $5==b {found=1}
                 END{exit !found}' "$cds_gff"; then
            echo "FAIL: INT-ANIMAL-01 annotation CDS at ${seqid}:${start}-${end} has no matching cds.gff row"
            mismatch=1
            FAILED=1
        fi
    done < <(awk -F'\t' '$3=="CDS" {print $1"\t"$4"\t"$5}' "$annotation_gff")
    if (( mismatch == 0 )); then
        echo "OK:   INT-ANIMAL-01 all annotation CDS features trace to cds.gff"
    fi
fi

plant_summary="$OUTDIR/INT-PLANT-01-pt/annotation/annotation_summary.json"
if [[ ! -s "$plant_summary" ]]; then
    echo "FAIL: INT-PLANT-01-pt annotation_summary.json missing or empty"
    FAILED=1
else
    status=$(jq -r '.status' "$plant_summary")
    non_cds_source=$(jq -r '.non_cds_source' "$plant_summary")
    n_cds=$(jq -r '.feature_counts.CDS' "$plant_summary")
    n_trna=$(jq -r '.feature_counts.tRNA' "$plant_summary")
    n_rrna=$(jq -r '.feature_counts.rRNA' "$plant_summary")

    if [[ "$status" != "ok_cds_only" ]]; then
        echo "FAIL: INT-PLANT-01-pt annotation status='$status', expected 'ok_cds_only'"
        FAILED=1
    else
        echo "OK:   INT-PLANT-01-pt annotation status='ok_cds_only'"
    fi
    if [[ "$non_cds_source" != "null" ]]; then
        echo "FAIL: INT-PLANT-01-pt non_cds_source='$non_cds_source', expected null"
        FAILED=1
    fi
    if (( n_cds > 0 )); then
        echo "OK:   INT-PLANT-01-pt annotation has ${n_cds} CDS features"
    else
        echo "FAIL: INT-PLANT-01-pt annotation has 0 CDS features"
        FAILED=1
    fi
    if (( n_trna == 0 && n_rrna == 0 )); then
        echo "OK:   INT-PLANT-01-pt annotation has 0 tRNA/rRNA (no plastid non-CDS annotator yet)"
    else
        echo "FAIL: INT-PLANT-01-pt annotation has tRNA=${n_trna} rRNA=${n_rrna}, expected 0/0"
        FAILED=1
    fi
fi

if [[ "$FAILED" -eq 0 ]]; then
    echo "All assertions passed."
fi
exit "$FAILED"
