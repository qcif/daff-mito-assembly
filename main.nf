#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

import groovy.json.JsonSlurper

include { VALIDATE_SAMPLESHEET     } from './modules/local/validate_samplesheet'
include { PARSE_SAMPLESHEET        } from './modules/local/parse_samplesheet'
include { NANOPLOT_RAW             } from './modules/local/nanoplot_raw'
include { DEDUPLICATE              } from './modules/local/deduplicate'
include { CHOPPER                  } from './modules/local/chopper'
include { FILTLONG                 } from './modules/local/filtlong'
include { NANOPLOT_CLEAN           } from './modules/local/nanoplot_clean'
include { RECRUIT                  } from './modules/local/recruit'
include { COVERAGE_GATE            } from './modules/local/coverage_gate'
include { METAFLYE                 } from './modules/local/metaflye'
include { MEDAKA                   } from './modules/local/medaka'
include { ALLOCATE_GRAPH_SENTINELS } from './modules/local/allocate_graph_sentinels'
include { BANDAGE_NG               } from './modules/local/bandage_ng'
include { ANNOTATE_GRAPH_SVG       } from './modules/local/annotate_graph_svg'
include { BIN_TARGET               } from './modules/local/bin_target'
include { BLAST_VALIDATE           } from './modules/local/blast_validate'
include { BARRNAP                  } from './modules/local/barrnap'
include { ANNOTATE                 } from './modules/local/annotate'
include { ANNOTATION_SCORING       } from './subworkflows/local/annotation_scoring'
include { MINIPROT_CDS             } from './modules/local/miniprot_cds'
include { ROTATE_ANCHORS            } from './modules/local/rotate_anchors'
include { ROTATE_TARGET             } from './modules/local/rotate_target'
include { SELECT_GENETIC_CODE      } from './modules/local/select_genetic_code'
include { EXTRACT_BARCODES         } from './modules/local/extract_barcodes'
include { ORGANELLE_MAP            } from './modules/local/organelle_map'
include { COLLATE                  } from './modules/local/collate'
include { RUN_REPORT               } from './modules/local/run_report'

// ---------------------------------------------------------------------------
// Pre-flight validation (null + existence checks before any process runs)
// ---------------------------------------------------------------------------

// Turn a samplesheet reads entry into a Path: a URI as-is, anything else
// relative to --data_dir (C1 guarantees --data_dir is set in that case).
def resolveRead(String entry) {
    return entry.contains('://')
        ? file(entry)
        : file(params.data_dir).resolve(entry)
}

// Existence check that names missing credentials rather than reporting
// every remote read as absent.
def existsOrFail(Path path) {
    try {
        return path.exists()
    }
    catch (Exception e) {
        def reason = "${e.class.simpleName}: ${e.message}"
        error "ERROR: cannot check '${path}' (${reason}) - for az:// reads, check the Azure credentials/profile"
    }
}

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
    if (params.data_dir && !file(params.data_dir).exists()) {
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

    def sizes = params.expected_target_size
    params.coverage_limits.keySet().each { target ->
        def entry = sizes?.get(target)
        if (entry == null) {
            error (
                "ERROR: params.expected_target_size.${target} is missing "
                + "— every assembly target needs min_bp/max_bp.")
        }
        if (!(entry.min_bp instanceof Integer) || entry.min_bp <= 0) {
            error (
                "ERROR: params.expected_target_size.${target}.min_bp must "
                + "be a positive integer.")
        }
        if (entry.max_bp != null && (
                !(entry.max_bp instanceof Integer)
                || entry.max_bp < entry.min_bp)) {
            error (
                "ERROR: params.expected_target_size.${target}.max_bp must "
                + "be null or an integer >= min_bp.")
        }
    }
    params.coverage_limits.keySet().each { target ->
        if (!params.tandem_collapse.containsKey(target)) {
            error (
                "ERROR: params.tandem_collapse.${target} is missing — "
                + "every assembly target needs an entry (null to disable).")
        }
        def cfg = params.tandem_collapse[target]
        if (cfg == null) { return }
        def prefix = "params.tandem_collapse.${target}"
        if (!(cfg.kmer instanceof Integer) || cfg.kmer < 15 || cfg.kmer > 63) {
            error "ERROR: ${prefix}.kmer must be an integer in [15, 63]."
        }
        if (!(cfg.min_support instanceof Number) || cfg.min_support <= 0
                || cfg.min_support > 1) {
            error "ERROR: ${prefix}.min_support must be a number in (0, 1]."
        }
        ['tolerance_bp', 'min_period_bp'].each { key ->
            if (!(cfg[key] instanceof Integer) || cfg[key] <= 0) {
                error "ERROR: ${prefix}.${key} must be a positive integer."
            }
        }
    }
    def fraction = params.partial_recovery_min_fraction
    if (!(fraction instanceof Number) || fraction <= 0 || fraction > 1) {
        error (
            "ERROR: params.partial_recovery_min_fraction must be a number "
            + "in (0, 1].")
    }
}

// `workflow.commitId` is null whenever the pipeline runs from a working
// directory rather than a cloned/tagged checkout (task 43a §2026-09-02
// outcome, tasks/todo.md "Run-level provenance"), which is COLLATE's
// (and, until now, would-be RUN_REPORT's) everyday case in local dev.
// Falls back to a live `git rev-parse HEAD` against the pipeline's own
// projectDir before giving up — resolves a real hash for that common
// case instead of shipping "unknown" in run_manifest.json's provenance
// (CONSTITUTION rule 16) whenever a real one is available.
def resolvePipelineCommit() {
    if (workflow.commitId) { return workflow.commitId }
    if (workflow.revision) { return workflow.revision }
    try {
        def proc = [
            'git', '-C', workflow.projectDir.toString(),
            'rev-parse', 'HEAD',
        ].execute()
        proc.waitForOrKill(5000)
        if (proc.exitValue() == 0) {
            return proc.text.trim()
        }
    } catch (Exception e) {
        // git unavailable, or projectDir isn't a git checkout — fall
        // through to the honest "unknown" rather than raising.
    }
    return 'unknown'
}

// ---------------------------------------------------------------------------
// Workflow
// ---------------------------------------------------------------------------

workflow {

    validateParams()

    // Workflow-wide value channels for shared reference files
    ch_samplesheet   = Channel.value(file(params.samplesheet))
    ch_organelle_refs = Channel.value(file(params.organelle_refs))
    ch_protein_panel  = Channel.value(file(params.protein_panel))
    ch_annotate_refs  = Channel.value(file(params.annotate_refs))
    ch_gene_sets      = Channel.value(file(params.gene_sets))
    ch_blast_db       = Channel.value(file(params.blast_db))
    ch_locus_panel    = Channel.value(file(params.locus_panel))
    ch_sample_metadata_schema =
        Channel.value(file(params.sample_metadata_schema))
    ch_run_manifest_schema =
        Channel.value(file(params.run_manifest_schema))
    // Bundle-root manifest (spec §4.4), distinct from the four
    // subdirectory params above — C6 needs the reference-bundle
    // version + generated_at for metadata.json's provenance section
    // (CONSTITUTION rule 16, task 42 §5.2), and C7 inlines the manifest
    // in full into run_manifest.json (task 45). Staged as a proper
    // `path` input (spec §1a pattern 1) rather than a bare ${params.x}
    // string, which can't be staged on a remote executor.
    ch_refs_manifest  = Channel.value(file(params.refs_manifest))
    // Report renderer assets (task 43a §5.1) — kept out of bin/ since
    // Nextflow stages the whole bin/ dir onto every process, and
    // static/js/ alone is 1.6 MB of vendored front-end libraries only
    // COLLATE and RUN_REPORT need.
    ch_report_templates =
        Channel.value(file("${projectDir}/scripts/report/templates"))
    ch_report_static =
        Channel.value(file("${projectDir}/scripts/report/static"))

    // Stage 0: validate CSV holistically; emit normalised JSON array.
    // A non-zero exit here aborts the run before any per-sample work starts.
    VALIDATE_SAMPLESHEET(ch_samplesheet, params.data_dir as boolean)

    // Resolve every reads entry on the launcher into a Path: relative
    // entries against --data_dir, URIs as-is. Missing files fail the whole
    // run once, listing every entry, before any per-sample work starts.
    ch_samples = VALIDATE_SAMPLESHEET.out.json
        | splitJson()
        | toList()
        | flatMap { rows ->
            def missing = []
            def resolved = rows.collect { row ->
                def reads = row.reads.collect { entry ->
                    def path = resolveRead(entry)
                    if (!existsOrFail(path)) {
                        missing << "${row.sample_id}: ${entry}"
                    }
                    path
                }
                [ row, reads ]
            }
            if (missing) {
                error "ERROR: reads not found:\n  " + missing.join("\n  ")
            }
            resolved
        }
        | map { row, reads ->
            def meta = [
                sample_id          : row.sample_id,
                assembly_target    : row.assembly_target,
                sample_info        : row.sample_info         ?: '',
                sample_type        : row.sample_type         ?: '',
                sample_receipt_date: row.sample_receipt_date ?: '',
                storage_location   : row.storage_location    ?: '',
            ]
            [ meta, reads ]
        }

    // Stage 0 (per-sample): check read integrity and concatenate
    // multi-file reads into one FASTQ (C1b). PARSE_SAMPLESHEET must remain
    // the single consumer of raw reads: on a remote executor each process
    // input is staged to its node, so any other process taking raw reads
    // would copy them again.
    PARSE_SAMPLESHEET(ch_samples)

    // Branch on intake.json: explicit allowlist, anything unrecognised
    // fails closed to the invalid side (spec §2.1.3, rule 18), where C6
    // turns it into an `invalid_input` sample. `dedup` is valid reads
    // that repeat read IDs identically (task 71).
    ch_intake = PARSE_SAMPLESHEET.out.intake
        .branch { meta, intake_json ->
            def intake = new JsonSlurper().parse(intake_json)
            dedup:   intake.status == 'valid' && intake.dedup_required == true
            valid:   intake.status == 'valid'
            invalid: true
        }

    ch_plain_reads = ch_intake.valid
        .join(PARSE_SAMPLESHEET.out.reads, by: 0)
        .map { meta, intake_json, reads -> [ meta, reads ] }

    // Stage 0b: drop repeated read IDs. NANOPLOT_RAW keeps the submitted
    // reads; everything else sees the deduplicated ones. The expected
    // removal count (intake's identical repeats) is cross-checked inside
    // the process.
    ch_dup_reads = ch_intake.dedup
        .join(PARSE_SAMPLESHEET.out.reads, by: 0)
        .map { meta, intake_json, reads ->
            [ meta, reads,
              new JsonSlurper().parse(intake_json).duplicates.identical ]
        }
    DEDUPLICATE(ch_dup_reads)

    ch_submitted_reads = ch_plain_reads.mix(
        ch_dup_reads.map { meta, reads, expected -> [ meta, reads ] })

    // Stages 1–4: QC
    NANOPLOT_RAW(ch_submitted_reads)
    CHOPPER(ch_plain_reads.mix(DEDUPLICATE.out.reads))
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

    // Stages 7–8: assemble + (optional) polish
    ch_for_assembly = ch_gated.ok
        .map { meta, gated_fq, status_json, coverage_json -> [ meta, gated_fq ] }

    METAFLYE(ch_for_assembly)

    ch_assembly = params.polish
        ? MEDAKA(METAFLYE.out.assembly, ch_for_assembly.map { it[1] }).assembly
        : METAFLYE.out.assembly

    // Stage 9: bin contigs — reads the assembly directly, not routed
    // through the graph chain (task 52: that plumbing only existed to
    // hand BIN_TARGET a 5th tuple element it never read).
    BIN_TARGET(ch_assembly, ch_organelle_refs)

    // Stage 10: node-labelled, hover-interactive assembly graph SVG
    // (task 47) — sentinel-colour round-trip across three processes
    // since BandageNG's biocontainer has no Python either side of the
    // render (see modules/local/bandage_ng.nf). ALLOCATE_GRAPH_
    // SENTINELS and BANDAGE_NG need only the assembly, so they run
    // concurrently with BIN_TARGET; only the final annotate step waits
    // on BIN_TARGET's bin_metadata.json to colour nodes by binning
    // bucket (task 52 §3.4) — this diagnostic chain is off the
    // critical path to BLAST_VALIDATE and beyond.
    ALLOCATE_GRAPH_SENTINELS(ch_assembly)
    BANDAGE_NG(ALLOCATE_GRAPH_SENTINELS.out.sentinels)
    ANNOTATE_GRAPH_SVG(
        BANDAGE_NG.out.rendered.join(BIN_TARGET.out.metadata, by: 0))

    // Stage 9a: rotate circular animal_mt contigs to start at COX1 (or
    // the next anchor gene) so annotation sees whole genes and every run
    // shares one coordinate frame. Must precede BLAST_VALIDATE and
    // MINIPROT_CDS (task 76).
    ROTATE_ANCHORS(BIN_TARGET.out.binned, ch_protein_panel)
    ROTATE_TARGET(
        BIN_TARGET.out.binned
            .join(BIN_TARGET.out.metadata, by: 0)
            .join(ROTATE_ANCHORS.out.anchors, by: 0)
            .map { meta, target_fasta, secondaries, bin_metadata, gffs ->
                [ meta, target_fasta, secondaries, bin_metadata, gffs,
                  params.rotation_anchors[meta.assembly_target],
                  params.tandem_collapse[meta.assembly_target] ]
            })

    // Stage 11: BLAST validation
    BLAST_VALIDATE(ROTATE_TARGET.out.binned, ch_blast_db)

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
    EXTRACT_BARCODES(ch_cds, ch_locus_panel)

    // Stage 13a / 14: annotate (merge cds.gff + MITOS2 non-CDS) + visualise
    ch_barrnap_in = ch_cds
        .filter { meta, target_fasta, cds_gff, genetic_code_json ->
            params.annotate[meta.assembly_target].non_cds_tool == 'barrnap'
        }
        .map { meta, target_fasta, cds_gff, genetic_code_json ->
            [ meta, target_fasta ] }
    BARRNAP(ch_barrnap_in)

    // Samples with no barrnap step get the shared placeholder file so
    // the ANNOTATE tuple has one shape (rule 13).
    ch_no_file = file("${projectDir}/assets/NO_FILE")
    ch_annotate_in = ch_cds
        .join(BARRNAP.out.result, by: 0, remainder: true)
        .map { row ->
            def (meta, target_fasta, cds_gff, genetic_code_json,
                 barrnap_gff, barrnap_exit) = row
            [ meta, target_fasta, cds_gff, genetic_code_json,
              barrnap_gff ?: ch_no_file, barrnap_exit ?: ch_no_file ]
        }
    ANNOTATE(ch_annotate_in, ch_annotate_refs, ch_gene_sets)

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

    // ORGANELLE_MAP is the sole stage aware of the plastid quadripartite
    // structure (spec §3.6 step 5) — it needs bin_metadata.json for the
    // true genome length and, on the plant_pt canonical branch, the
    // LSC/IR/SSC segment lengths to remap features onto path2.
    ORGANELLE_MAP(
        ANNOTATION_SCORING.out.annotation
            .join(BIN_TARGET.out.metadata, by: 0),
        ch_locus_panel)

    // Stage 15: collate per-sample bundle (C6, task 42).
    // ANNOTATE_GRAPH_SVG's node-labelled graph SVG, METAFLYE's
    // assembly_info.txt and the genetic-code selection (MINIPROT_CDS
    // single-table / C9 clade-trial) all already flow *through* the
    // assembly chain but were previously dropped before reaching
    // COLLATE (task 42 §2.2) — recovered here rather than re-run.
    ch_graph_svg = ANNOTATE_GRAPH_SVG.out.graph
    ch_assembly_info = METAFLYE.out.assembly
        .map { meta, fasta, gfa, info -> [ meta, info ] }
    ch_genetic_code = ch_cds
        .map { meta, target_fasta, cds_gff, genetic_code_json ->
            [ meta, genetic_code_json ] }

    // QC-only stages (NanoPlot) join with `remainder: true` so a failed
    // report leaves an empty slot instead of dropping the sample. A
    // remainder join also emits right-hand items with no left-hand
    // match (a null in position 1, the gate/intake slot); those belong
    // to samples on the other gate branch and are filtered out
    // straight after each join.
    ch_ok_inputs = ch_gated.ok
        .map { meta, gated_fq, status_json, coverage_json ->
            [ meta, status_json, coverage_json ]
        }
        .join(PARSE_SAMPLESHEET.out.intake,  by: 0)
        .map { meta, status_json, coverage_json, intake_json ->
            [ meta, intake_json, status_json, coverage_json ]
        }
        .join(NANOPLOT_RAW.out.reports,      by: 0, remainder: true)
        .filter { it[1] != null }
        .join(NANOPLOT_CLEAN.out.reports,    by: 0, remainder: true)
        .filter { it[1] != null }
        .join(DEDUPLICATE.out.stats,         by: 0, remainder: true)
        .filter { it[1] != null }
        .join(RECRUIT.out.stats,             by: 0)
        .join(ROTATE_TARGET.out.binned,      by: 0)
        .join(BLAST_VALIDATE.out.validated,  by: 0)
        .join(EXTRACT_BARCODES.out.barcodes, by: 0)
        .join(ANNOTATION_SCORING.out.annotation, by: 0)
        .join(ORGANELLE_MAP.out.map,         by: 0)
        .join(ch_graph_svg,                  by: 0)
        .join(BIN_TARGET.out.metadata,       by: 0)
        .join(ch_assembly_info,              by: 0)
        .join(ch_genetic_code,               by: 0)
        .join(ROTATE_TARGET.out.rotation,    by: 0)
        // BIN_TARGET.out.isoforms only emits on the plant_pt canonical
        // branch (task 24 §3.2) — `remainder: true` keeps every other
        // sample in the join instead of silently dropping it, with the
        // missing isoforms field filled in as null below. Isoforms exist
        // only for samples on this branch, so nothing leaks right-only.
        .join(BIN_TARGET.out.isoforms, by: 0, remainder: true)
        .map { meta, intake_json, status_json, coverage_json,
               nanoplot_raw, nanoplot_clean, dedup_stats, recruit_stats,
               target_fasta, secondaries,
               target_fasta2, blast_tsv,
               barcodes_fasta, coords_gff, validation_tsv,
               annotation_gff, annotation_summary,
               organelle_map_svg,
               graph_svg, bin_metadata_json, assembly_info,
               genetic_code_json, rotation_json, isoforms ->
            [ meta, intake_json, status_json, coverage_json,
              nanoplot_raw ?: [], nanoplot_clean ?: [], dedup_stats ?: [],
              recruit_stats,
              target_fasta, secondaries, blast_tsv,
              barcodes_fasta, coords_gff, validation_tsv,
              annotation_gff, annotation_summary, organelle_map_svg,
              graph_svg, bin_metadata_json, assembly_info,
              genetic_code_json, rotation_json, isoforms ?: [] ]
        }

    ch_failed_inputs = ch_gated.failed
        .join(PARSE_SAMPLESHEET.out.intake,  by: 0)
        .map { meta, gated_fq, status_json, coverage_json, intake_json ->
            [ meta, intake_json, status_json, coverage_json ]
        }
        .join(NANOPLOT_RAW.out.reports,      by: 0, remainder: true)
        .filter { it[1] != null }
        .join(NANOPLOT_CLEAN.out.reports,    by: 0, remainder: true)
        .filter { it[1] != null }
        .join(DEDUPLICATE.out.stats,         by: 0, remainder: true)
        .filter { it[1] != null }
        .join(RECRUIT.out.stats,             by: 0)
        .map {
            meta, intake_json, status_json, coverage_json,
            nanoplot_raw, nanoplot_clean, dedup_stats, recruit_stats ->
            [ meta, intake_json, status_json, coverage_json,
              nanoplot_raw ?: [], nanoplot_clean ?: [], dedup_stats ?: [],
              recruit_stats,
              [], [], [], [], [], [], [], [], [],
              [], [], [], [], [], [] ]
        }

    // Samples whose reads failed the intake check go straight to COLLATE
    // with intake.json alone: no QC, recruitment or gate ran.
    ch_invalid_inputs = ch_intake.invalid
        .map { meta, intake_json ->
            [ meta, intake_json ] + (1..21).collect { [] }
        }

    COLLATE(
        ch_ok_inputs.mix(ch_failed_inputs, ch_invalid_inputs),
        ch_gene_sets, ch_locus_panel, ch_sample_metadata_schema,
        ch_refs_manifest, ch_report_templates, ch_report_static)

    // Stage 16: run-level report + manifest
    // `ifEmpty([])` keeps RUN_REPORT running when no sample reached
    // COLLATE: every samplesheet row then reconciles as `error`.
    ch_metadata = COLLATE.out.bundle
        .map { meta, metadata, report -> metadata }
        .collect()
        .ifEmpty([])

    RUN_REPORT(
        ch_metadata, ch_samplesheet, ch_refs_manifest,
        ch_run_manifest_schema, ch_report_templates, ch_report_static,
        resolvePipelineCommit())
}

// ---------------------------------------------------------------------------
// End-of-run signal for tasks that failed but were ignored (per-sample crash
// isolation, conf/base.config). Exit status is unchanged.
// ---------------------------------------------------------------------------

workflow.onComplete {
    def ignored = workflow.stats.ignoredCount
    if (ignored > 0) {
        log.warn(
            "\n" + '=' * 72 + "\n" +
            "WF5 WARNING: ${ignored} task(s) failed and were ignored.\n" +
            "The affected samples are listed in run-report.html and " +
            "run_manifest.json\n(status 'error' or 'invalid_input'). " +
            "Logs for the failed tasks are in ${params.outdir}/errors/.\n" +
            "To list the failed tasks with their sample tags:\n" +
            "  nextflow log ${workflow.runName} -f name,status,workdir " +
            "-F \"status == 'FAILED'\"\n" +
            '=' * 72)
    }
}
