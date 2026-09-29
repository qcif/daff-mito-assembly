// Stage 16 — cross-sample run summary + manifest (C7).
// C7 custom logic — see plan.md §2.2, §6a.4, task 45.
//
// Join point across all samples: does not care whether each sample
// succeeded or soft-failed, only that a metadata.json exists (spec
// §2.1.3). COLLATE always emits one, even on the `fail` / `no_assembly`
// branches (modules/local/collate.nf), so a samplesheet sample_id with
// no matching bundle here crashed somewhere upstream of COLLATE —
// run_report.py reconciles against the samplesheet to catch it and
// reports it as `error`, the one sample_status value C7 derives itself
// (CONSTITUTION principle 7).
//
// `pipeline_commit` is resolved in main.nf (resolvePipelineCommit())
// rather than inline here — it falls back to a live `git rev-parse`
// against the pipeline's own checkout when `workflow.commitId` is null
// (the everyday local-dev case), so run_manifest.json's provenance
// carries a real commit hash wherever one is available (tasks/todo.md
// "Run-level provenance").

process RUN_REPORT {
    label        'process_low'
    publishDir   "${params.outdir}", mode: 'copy'

    input:
    // stageAs '?/metadata.json' gives each file a unique subdir (0/, 1/, ...) so
    // Nextflow doesn't complain about filename collisions across samples.
    path bundles, stageAs: '?/metadata.json'
    path samplesheet
    path refs_manifest
    path run_manifest_schema
    path report_templates
    path report_static
    val  pipeline_commit

    output:
    path "run_manifest.json", emit: manifest
    path "run-report.html",   emit: report

    script:
    """
    run_report.py \\
        --metadata-json */metadata.json \\
        --samplesheet ${samplesheet} \\
        --refs-manifest ${refs_manifest} \\
        --pipeline-commit '${pipeline_commit}' \\
        --workflow-start '${workflow.start.format("yyyy-MM-dd'T'HH:mm:ss")}' \\
        --schema ${run_manifest_schema} \\
        --report-templates ${report_templates} \\
        --report-static ${report_static}
    """

    stub:
    """
    touch run_manifest.json run-report.html
    """
}
