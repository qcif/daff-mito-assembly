// Stage 15 — per-sample bundle aggregation + metadata.json (C6).
// C6 custom logic — see spec §2 stage 15, §6a.5, task 42.
//
// Derives the five-value `sample_status` (CONSTITUTION.md principle 7)
// from sample_status.json (gate) / bin_metadata.json (assembly) /
// validation.tsv (barcodes) and dispatches to the full or minimal
// per-sample bundle (task 42 §3). `report.html` is rendered from
// metadata.json by bin/report/ (task 43a) — templates/static are
// staged in from `report_templates`/`report_static` rather than
// bin/ (rule 13: bin/ is staged onto every process, and the vendored
// front-end assets are 1.6 MB only this process needs).
//
// Three separate output declarations, not one: metadata.json +
// report.html always exist; organelle_assembly.fasta /
// organelle_annotation.gff / barcodes.fasta exist together or not at
// all (full vs. minimal bundle); diagnostics/ exists whenever any
// upstream stage ran, independently of bundle kind. `optional: true`
// sits at the TUPLE level for each of the latter two — per-path
// optional inside a tuple is a no-op in Nextflow 25.10.2 (confirmed
// with a standalone repro during task 38); see
// modules/local/miniprot_cds.nf's `resolved` output for the same
// pattern. Folding all five paths into one tuple with tuple-level
// optional would make metadata.json/report.html vanish along with the
// bundle on the minimal path, which would silently break RUN_REPORT
// (task 45) — every sample, including `fail`/`no_assembly`, must
// still emit a `metadata.json` for C7 to read.

process COLLATE {
    tag          "${meta.sample_id}"
    label        'process_low'
    publishDir   "${params.outdir}/${meta.sample_id}",
                 mode: 'copy'

    input:
    tuple val(meta),
          path(status_json),
          path(coverage_json),
          path(nanoplot_raw),
          path(nanoplot_clean),
          path(recruit_stats),
          path(target_fasta),
          path(secondaries_tsv),
          path(blast_tsv),
          path(barcodes_fasta),
          path(coords_gff),
          path(validation_tsv),
          path(annotation_gff),
          path(annotation_summary),
          path(organelle_map_svg),
          path(graph_png),
          path(bin_metadata_json),
          path(assembly_info),
          path(genetic_code_json),
          path(plastid_isoforms)
    path gene_sets
    path sample_metadata_schema
    path refs_manifest
    path report_templates
    path report_static

    output:
    tuple val(meta), path("metadata.json"), path("report.html"), emit: bundle
    tuple val(meta),
          path("organelle_assembly.fasta"),
          path("organelle_annotation.gff"),
          path("barcodes.fasta"),
          optional: true, emit: full_bundle
    tuple val(meta), path("diagnostics"), optional: true, emit: diagnostics

    script:
    def optArg = { flag, val -> val ? "--${flag} ${val}" : '' }
    // Written to a file rather than interpolated onto the command line:
    // sample_info/sample_type/storage_location are unconstrained
    // submitter free text (spec §0), and shell-quoting arbitrary text
    // directly into the script block would be a command-injection risk.
    def metaJson = groovy.json.JsonOutput.toJson([
        sample_id          : meta.sample_id,
        assembly_target    : meta.assembly_target,
        sample_info        : meta.sample_info,
        sample_type        : meta.sample_type,
        sample_receipt_date: meta.sample_receipt_date,
        storage_location   : meta.storage_location,
    ])
    // Resolved workflow params, for the report's "view all parameters"
    // modal (task 43a §4.4) — same command-injection reasoning as
    // metaJson above: params values are arbitrary and must not be
    // shell-interpolated.
    def paramsJson = groovy.json.JsonOutput.toJson(params)
    """
    cat <<'META_EOF' > meta.json
${metaJson}
META_EOF

    cat <<'PARAMS_EOF' > params.json
${paramsJson}
PARAMS_EOF

    collate.py \\
        --meta-json meta.json \\
        --pipeline-commit '${workflow.commitId ?: workflow.revision ?: "unknown"}' \\
        --workflow-start '${workflow.start.format("yyyy-MM-dd'T'HH:mm:ss")}' \\
        --status-json ${status_json} \\
        --coverage-json ${coverage_json} \\
        --recruit-stats ${recruit_stats} \\
        --nanoplot-raw ${nanoplot_raw} \\
        --nanoplot-clean ${nanoplot_clean} \\
        --gene-sets ${gene_sets} \\
        --schema ${sample_metadata_schema} \\
        --params-json params.json \\
        --report-templates ${report_templates} \\
        --report-static ${report_static} \\
        ${optArg('bin-metadata-json', bin_metadata_json)} \\
        ${optArg('assembly-info', assembly_info)} \\
        ${optArg('target-fasta', target_fasta)} \\
        ${optArg('secondaries-tsv', secondaries_tsv)} \\
        ${optArg('blast-tsv', blast_tsv)} \\
        ${optArg('barcodes-fasta', barcodes_fasta)} \\
        ${optArg('validation-tsv', validation_tsv)} \\
        ${optArg('annotation-gff', annotation_gff)} \\
        ${optArg('annotation-summary', annotation_summary)} \\
        ${optArg('genetic-code-json', genetic_code_json)} \\
        ${optArg('organelle-map-svg', organelle_map_svg)} \\
        ${optArg('graph-png', graph_png)} \\
        ${optArg('plastid-isoforms', plastid_isoforms)} \\
        ${optArg('refs-manifest', refs_manifest)}
    """

    stub:
    """
    touch metadata.json report.html
    """
}
