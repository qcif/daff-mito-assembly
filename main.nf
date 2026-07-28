#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VALIDATE_SAMPLESHEET } from './modules/local/validate_samplesheet'
include { PARSE_SAMPLESHEET    } from './modules/local/parse_samplesheet'
include { NANOPLOT_RAW         } from './modules/local/nanoplot_raw'
include { CHOPPER              } from './modules/local/chopper'
include { FILTLONG             } from './modules/local/filtlong'
include { NANOPLOT_CLEAN       } from './modules/local/nanoplot_clean'
include { RECRUIT              } from './modules/local/recruit'
include { COVERAGE_GATE        } from './modules/local/coverage_gate'
include { METAFLYE             } from './modules/local/metaflye'
include { MEDAKA               } from './modules/local/medaka'
include { BANDAGE_NG           } from './modules/local/bandage_ng'
include { BIN_TARGET           } from './modules/local/bin_target'
include { BLAST_VALIDATE       } from './modules/local/blast_validate'
include { ANNOTATE             } from './modules/local/annotate'
include { MINIPROT_EXTRACT     } from './modules/local/miniprot_extract'
include { ORGANELLE_MAP        } from './modules/local/organelle_map'
include { COLLATE              } from './modules/local/collate'
include { RUN_REPORT           } from './modules/local/run_report'

// ---------------------------------------------------------------------------
// Pre-flight validation (null + existence checks before any process runs)
// ---------------------------------------------------------------------------

def validateParams() {
    if (!params.samplesheet) { error "ERROR: --samplesheet is required" }
    if (!params.data_dir)    { error "ERROR: --data_dir is required" }
    if (!file(params.data_dir).exists()) {
        error "ERROR: --data_dir '${params.data_dir}' does not exist"
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

    // Stage 6: coverage gate — soft-fail | passthrough | subsample
    COVERAGE_GATE(RECRUIT.out.reads)

    // Branch: ok → assembly chain; failed → COLLATE directly
    ch_gated = COVERAGE_GATE.out.gated
        .branch { meta, gated_fq, status_json, coverage_json ->
            ok:     status_json.text.contains('"status": "ok"')
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

    // Stages 12–14: annotate + extract barcodes + visualise
    ANNOTATE(BLAST_VALIDATE.out.validated)
    MINIPROT_EXTRACT(BLAST_VALIDATE.out.validated)
    ORGANELLE_MAP(ANNOTATE.out.annotation)

    // Stage 15: collate per-sample bundle
    ch_ok_gate_meta = ch_gated.ok
        .map { meta, gated_fq, status_json, coverage_json ->
            [ meta, status_json, coverage_json ]
        }

    ch_ok_inputs = ch_ok_gate_meta
        .join(NANOPLOT_RAW.out.reports,      by: 0)
        .join(NANOPLOT_CLEAN.out.reports,    by: 0)
        .join(BIN_TARGET.out.binned,         by: 0)
        .join(BLAST_VALIDATE.out.validated,  by: 0)
        .join(MINIPROT_EXTRACT.out.barcodes, by: 0)
        .join(ANNOTATE.out.annotation,       by: 0)
        .join(ORGANELLE_MAP.out.map,         by: 0)
        .map { meta, status_json, coverage_json,
               nanoplot_raw, nanoplot_clean,
               target_fasta, secondaries,
               target_fasta2, blast_tsv,
               barcodes_fasta, coords_gff, validation_tsv,
               annotation_gff, annotation_gbk,
               organelle_map_svg ->
            [ meta, status_json, coverage_json,
              nanoplot_raw, nanoplot_clean,
              target_fasta, blast_tsv,
              barcodes_fasta, coords_gff, validation_tsv,
              annotation_gff, annotation_gbk, organelle_map_svg,
              [] ]
        }

    ch_failed_inputs = ch_gated.failed
        .join(NANOPLOT_RAW.out.reports,   by: 0)
        .join(NANOPLOT_CLEAN.out.reports, by: 0)
        .map { meta, gated_fq, status_json, coverage_json, nanoplot_raw, nanoplot_clean ->
            [ meta, status_json, coverage_json,
              nanoplot_raw, nanoplot_clean,
              [], [], [], [], [], [], [], [], [] ]
        }

    COLLATE(ch_ok_inputs.mix(ch_failed_inputs))

    // Stage 16: run-level report + manifest
    ch_metadata = COLLATE.out.bundle
        .map { bundle -> bundle[4] }
        .collect()

    RUN_REPORT(ch_metadata, ch_samplesheet)
}
