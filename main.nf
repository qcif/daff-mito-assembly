#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

import groovy.json.JsonSlurper

include { VALIDATE_SAMPLESHEET     } from './modules/local/validate_samplesheet'
include { PARSE_SAMPLESHEET        } from './modules/local/parse_samplesheet'
include { NANOPLOT_RAW             } from './modules/local/nanoplot_raw'
include { CHOPPER                  } from './modules/local/chopper'
include { FILTLONG                 } from './modules/local/filtlong'
include { NANOPLOT_CLEAN           } from './modules/local/nanoplot_clean'
include { RECRUIT                  } from './modules/local/recruit'
include { COVERAGE_GATE            } from './modules/local/coverage_gate'
include { METAFLYE                 } from './modules/local/metaflye'
include { MEDAKA                   } from './modules/local/medaka'
include { BANDAGE_NG               } from './modules/local/bandage_ng'
include { BIN_TARGET               } from './modules/local/bin_target'
include { BLAST_VALIDATE           } from './modules/local/blast_validate'
include { ANNOTATE                 } from './modules/local/annotate'
include { ANNOTATION_SCORING       } from './subworkflows/local/annotation_scoring'
include { MINIPROT_CDS             } from './modules/local/miniprot_cds'
include { SELECT_GENETIC_CODE      } from './modules/local/select_genetic_code'
include { EXTRACT_BARCODES         } from './modules/local/extract_barcodes'
include { ORGANELLE_MAP            } from './modules/local/organelle_map'
include { COLLATE                  } from './modules/local/collate'
include { RUN_REPORT               } from './modules/local/run_report'

// ---------------------------------------------------------------------------
// Pre-flight validation (null + existence checks before any process runs)
// ---------------------------------------------------------------------------

// Reference-data params, each pointing at a directory the refdata bundle
// is expected to supply (scripts/fetch_refs.sh). A missing one must abort
// here: Nextflow stages a directory input as a symlink, and `ln -s`
// succeeds on a target that does not exist, so the failure would
// otherwise surface deep inside a container as a tool-specific "no such
// directory" — which is how an incomplete bundle silently cost ANNOTATE
// every tRNA and rRNA for nine nightly runs.
def validateParams() {
    def refdataParams = [
        'organelle_refs', 'protein_panel', 'blast_db', 'annotate_refs',
        'refs_manifest',
    ]

    if (!params.samplesheet) { error "ERROR: --samplesheet is required" }
    if (!params.data_dir)    { error "ERROR: --data_dir is required" }
    if (!file(params.data_dir).exists()) {
        error "ERROR: --data_dir '${params.data_dir}' does not exist"
    }

    refdataParams.each { name ->
        def value = params[name]
        if (!value) { error "ERROR: --${name} is required" }
        if (!file(value).exists()) {
            error (
                "ERROR: --${name} '${value}' does not exist — the refdata "
                + "bundle is missing or incomplete (see scripts/fetch_refs.sh)")
        }
    }

    // `min_cov` was retired in task 35 (spec §2.1.1): hard_min_cov and
    // warn_cov each carry part of its old meaning, so a config still
    // setting min_cov is ambiguous and must not be silently absorbed
    // into either floor.
    params.coverage_limits.each { target, limits ->
        if (limits.containsKey('min_cov')) {
            error (
                "ERROR: params.coverage_limits.${target}.min_cov was "
                + "retired in task 35. Set hard_min_cov (terminal floor) "
                + "and warn_cov (advisory floor) instead — see spec "
                + "§2.1.1.")
        }
        if (limits.hard_min_cov > limits.warn_cov) {
            error (
                "ERROR: hard_min_cov must not exceed warn_cov for "
                + "${target} — see spec §2.1.1.")
        }
    }
}

// ---------------------------------------------------------------------------
// Workflow
// ---------------------------------------------------------------------------

workflow {

    validateParams()

    // Workflow-wide value channels for shared reference files
    ch_samplesheet   = Channel.value(file(params.samplesheet))
    ch_data_dir      = Channel.value(file(params.data_dir))
    ch_organelle_refs = Channel.value(file(params.organelle_refs))
    ch_protein_panel  = Channel.value(file(params.protein_panel))
    ch_annotate_refs  = Channel.value(file(params.annotate_refs))
    ch_gene_sets      = Channel.value(file(params.gene_sets))
    ch_sample_metadata_schema =
        Channel.value(file(params.sample_metadata_schema))
    // Bundle-root manifest (spec §4.4), distinct from the four
    // subdirectory params above — C6 needs the reference-bundle
    // version + generated_at for metadata.json's provenance section
    // (CONSTITUTION rule 16, task 42 §5.2). Staged as a proper `path`
    // input (spec §1a pattern 1) rather than a bare ${params.x}
    // string, which can't be staged on a remote executor.
    ch_refs_manifest  = Channel.value(file(params.refs_manifest))

    // Stage 0: validate CSV holistically; emit normalised JSON array.
    // A non-zero exit here aborts the run before any per-sample work starts.
    VALIDATE_SAMPLESHEET(ch_samplesheet, ch_data_dir)

    // Fan-out: one (meta, reads) tuple per validated sample row.
    // Reads paths in the JSON are already resolved absolute paths.
    ch_samples = VALIDATE_SAMPLESHEET.out.json
        | splitJson()
        | map { row ->
            def meta = [
                sample_id          : row.sample_id,
                assembly_target    : row.assembly_target,
                sample_info        : row.sample_info         ?: '',
                sample_type        : row.sample_type         ?: '',
                sample_receipt_date: row.sample_receipt_date ?: '',
                storage_location   : row.storage_location    ?: '',
            ]
            def reads = row.reads.collect { file(it) }
            [ meta, reads ]
        }

    // Stage 0 (per-sample): concatenate multi-file reads into one FASTQ
    PARSE_SAMPLESHEET(ch_samples)

    // Stages 1–4: QC
    NANOPLOT_RAW(PARSE_SAMPLESHEET.out.reads)
    CHOPPER(PARSE_SAMPLESHEET.out.reads)
    FILTLONG(CHOPPER.out.reads)
    NANOPLOT_CLEAN(FILTLONG.out.reads)

    // Stage 5: positive recruitment against target organelle reference
    RECRUIT(FILTLONG.out.reads, ch_organelle_refs)

    // Stage 6: coverage gate — fail | degraded passthrough | passthrough | subsample
    COVERAGE_GATE(RECRUIT.out.reads, ch_organelle_refs)

    // Branch: ok/low_coverage → assembly chain; fail → COLLATE directly.
    // Parsed explicit allowlist, not a substring/prefix test (spec
    // §2.1.3, task 35): the three statuses share no common prefix, and a
    // status added later must not land on whichever side its spelling
    // happens to match (rule 18). Unrecognised statuses fail closed.
    ch_gated = COVERAGE_GATE.out.gated
        .branch { meta, gated_fq, status_json, coverage_json ->
            def status = new JsonSlurper().parse(status_json).status
            ok:     status in ['ok', 'low_coverage']
            failed: true
        }

    // Stages 7–9: assemble + (optional) polish + graph viz
    ch_for_assembly = ch_gated.ok
        .map { meta, gated_fq, status_json, coverage_json -> [ meta, gated_fq ] }

    METAFLYE(ch_for_assembly)

    ch_assembly = params.polish
        ? MEDAKA(METAFLYE.out.assembly, ch_for_assembly.map { it[1] }).assembly
        : METAFLYE.out.assembly

    BANDAGE_NG(ch_assembly)

    // Stage 10: bin contigs
    BIN_TARGET(BANDAGE_NG.out.assembly, ch_organelle_refs)

    // Stage 11: BLAST validation
    BLAST_VALIDATE(BIN_TARGET.out.binned)

    // Stage 12: one miniprot pass per configured genetic-code table over
    // the comprehensive protein panel — feeds both EXTRACT_BARCODES
    // (stage 13) and ANNOTATE (stage 13a, task 31). animal_mt has >1
    // configured table (task 38 §2), so its candidates route through
    // SELECT_GENETIC_CODE (C9); single-table targets (plant_pt,
    // plant_mt) are already resolved by MINIPROT_CDS itself.
    MINIPROT_CDS(BLAST_VALIDATE.out.validated, ch_protein_panel)

    ch_cds_candidates = MINIPROT_CDS.out.candidates
        .branch { meta, target_fasta, gffs ->
            multi:  params.genetic_code_tables[meta.assembly_target].size() > 1
            single: true
        }

    SELECT_GENETIC_CODE(ch_cds_candidates.multi)

    ch_cds_single = MINIPROT_CDS.out.resolved
        .filter { meta, cds_gff, genetic_code_json ->
            params.genetic_code_tables[meta.assembly_target].size() == 1
        }

    ch_cds = BLAST_VALIDATE.out.validated
        .map { meta, target_fasta, blast_tsv -> [ meta, target_fasta ] }
        .join(ch_cds_single.mix(SELECT_GENETIC_CODE.out.selected), by: 0)

    // Stage 13: barcode subset of stage 12's CDS features + ORF validation
    EXTRACT_BARCODES(ch_cds)

    // Stage 13a / 14: annotate (merge cds.gff + MITOS2 non-CDS) + visualise
    ANNOTATE(ch_cds, ch_annotate_refs)

    // Stage 13a cont'd — confidence scoring (task 40): translate every
    // CDS feature uniformly, blastp it against the same per-gene
    // protein panel MINIPROT_CDS aligned to, and join pident/qcovhsp/
    // bitscore back onto the GFF + annotation_summary.json. The
    // subworkflow's output is the annotation stage's final, published
    // artifact.
    ANNOTATION_SCORING(
        ANNOTATE.out.annotation,
        ch_cds.map { meta, target_fasta, cds_gff, genetic_code_json ->
            [ meta, target_fasta ] },
        ch_protein_panel)

    ORGANELLE_MAP(ANNOTATION_SCORING.out.annotation)

    // Stage 15: collate per-sample bundle (C6, task 42).
    // BANDAGE_NG's graph PNG, METAFLYE's assembly_info.txt and the
    // genetic-code selection (MINIPROT_CDS single-table / C9
    // clade-trial) all already flow *through* the assembly chain but
    // were previously dropped before reaching COLLATE (task 42 §2.2) —
    // recovered here rather than re-run.
    ch_graph_png = BANDAGE_NG.out.assembly
        .map { meta, assembly, gfa, info, graph_png -> [ meta, graph_png ] }
    ch_assembly_info = METAFLYE.out.assembly
        .map { meta, fasta, gfa, info -> [ meta, info ] }
    ch_genetic_code = ch_cds
        .map { meta, target_fasta, cds_gff, genetic_code_json ->
            [ meta, genetic_code_json ] }

    ch_ok_gate_meta = ch_gated.ok
        .map { meta, gated_fq, status_json, coverage_json ->
            [ meta, status_json, coverage_json ]
        }

    ch_ok_inputs = ch_ok_gate_meta
        .join(NANOPLOT_RAW.out.reports,      by: 0)
        .join(NANOPLOT_CLEAN.out.reports,    by: 0)
        .join(BIN_TARGET.out.binned,         by: 0)
        .join(BLAST_VALIDATE.out.validated,  by: 0)
        .join(EXTRACT_BARCODES.out.barcodes, by: 0)
        .join(ANNOTATION_SCORING.out.annotation, by: 0)
        .join(ORGANELLE_MAP.out.map,         by: 0)
        .join(ch_graph_png,                  by: 0)
        .join(BIN_TARGET.out.metadata,       by: 0)
        .join(ch_assembly_info,              by: 0)
        .join(ch_genetic_code,               by: 0)
        // BIN_TARGET.out.isoforms only emits on the plant_pt canonical
        // branch (task 24 §3.2) — `remainder: true` keeps every other
        // sample in the join instead of silently dropping it, with the
        // missing isoforms field filled in as null below.
        .join(BIN_TARGET.out.isoforms, by: 0, remainder: true)
        .map { meta, status_json, coverage_json,
               nanoplot_raw, nanoplot_clean,
               target_fasta, secondaries,
               target_fasta2, blast_tsv,
               barcodes_fasta, coords_gff, validation_tsv,
               annotation_gff, annotation_summary,
               organelle_map_svg,
               graph_png, bin_metadata_json, assembly_info,
               genetic_code_json, isoforms ->
            [ meta, status_json, coverage_json,
              nanoplot_raw, nanoplot_clean,
              target_fasta, secondaries, blast_tsv,
              barcodes_fasta, coords_gff, validation_tsv,
              annotation_gff, annotation_summary, organelle_map_svg,
              graph_png, bin_metadata_json, assembly_info,
              genetic_code_json, isoforms ?: [] ]
        }

    ch_failed_inputs = ch_gated.failed
        .join(NANOPLOT_RAW.out.reports,   by: 0)
        .join(NANOPLOT_CLEAN.out.reports, by: 0)
        .map { meta, gated_fq, status_json, coverage_json, nanoplot_raw, nanoplot_clean ->
            [ meta, status_json, coverage_json,
              nanoplot_raw, nanoplot_clean,
              [], [], [], [], [], [], [], [], [],
              [], [], [], [], [] ]
        }

    COLLATE(
        ch_ok_inputs.mix(ch_failed_inputs),
        ch_gene_sets, ch_sample_metadata_schema, ch_refs_manifest)

    // Stage 16: run-level report + manifest
    ch_metadata = COLLATE.out.bundle
        .map { meta, metadata, report -> metadata }
        .collect()

    RUN_REPORT(ch_metadata, ch_samplesheet)
}
