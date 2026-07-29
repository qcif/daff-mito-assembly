// Stage 6 — coverage estimation + subsample/soft-fail gate.
// C2 custom logic — see plan.md §2.2 and §2.1.
// Always exits 0; gate decision is written to sample_status.json (data, not error).
// errorStrategy 'ignore' guards against unexpected seqkit/seqtk crashes only.

process COVERAGE_GATE {
    tag            "${meta.sample_id}"
    label          'process_low'
    errorStrategy  'ignore'
    // sample_status.json always published — needed by integration assertions
    // and later re-bundled by COLLATE.
    publishDir     "${params.outdir}/${meta.sample_id}/coverage_gate",
                   mode: 'copy', pattern: 'sample_status.json'
    publishDir     "${params.outdir}/${meta.sample_id}/coverage_gate",
                   mode: 'copy', enabled: params.publish_intermediates,
                   pattern: '*.fastq.gz,coverage.json'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.sample_id}.gated.fastq.gz"),
          path("sample_status.json"), path("coverage.json"), emit: gated

    stub:
    """
    touch ${meta.sample_id}.gated.fastq.gz
    echo '{"status": "ok", "estimated_cov": 150.0}' > sample_status.json
    echo '{"pre_subsample_cov": 150.0, "post_subsample_cov": 150.0}' > coverage.json
    """

    script:
    def limits = params.coverage_limits[meta.assembly_target]
    """
    coverage_gate.py \\
        --reads ${reads} \\
        --sample-id ${meta.sample_id} \\
        --nominal-size ${limits.nominal_size} \\
        --min-cov ${limits.min_cov} \\
        --max-cov ${limits.max_cov} \\
        --seed ${params.seqtk_seed} \\
        --out-fastq ${meta.sample_id}.gated.fastq.gz
    """
}
