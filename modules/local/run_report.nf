// Stage 16 — cross-sample run summary + manifest.
// C7 custom logic — see plan.md §2.2, §6a.4.
// Join point: collects all per-sample bundles, emits run-level outputs.

process RUN_REPORT {
    label        'process_low'
    publishDir   "${params.outdir}", mode: 'copy'

    input:
    // stageAs '?/metadata.json' gives each file a unique subdir (0/, 1/, ...) so
    // Nextflow doesn't complain about filename collisions across samples.
    path bundles, stageAs: '?/metadata.json'
    path samplesheet

    output:
    path "run_manifest.json", emit: manifest
    path "run-report.html",   emit: report

    script:
    """
    # STUB — real implementation in P4 (bin/run_report.py)
    # Reads all metadata.json files, classifies sample statuses,
    # emits run_manifest.json + run-report.html.
    touch run_manifest.json run-report.html
    """

    stub:
    """
    touch run_manifest.json run-report.html
    """
}
