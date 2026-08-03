// Stage 6 — coverage estimation + subsample/soft-fail gate.
// C2 custom logic — see plan.md §2.2 and §2.1.
// Where the declared target has a sibling organelle panel, recruited
// reads are split by panel and only target-assigned bases feed the
// estimate (spec §2.1.5, task 25) — hence the organelle_refs input.
// Always exits 0; gate decision is written to sample_status.json (data, not error).
// errorStrategy 'ignore' guards against unexpected seqkit/seqtk/minimap2 crashes only.

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
    path organelle_refs

    output:
    tuple val(meta), path("${meta.sample_id}.gated.fastq.gz"),
          path("sample_status.json"), path("coverage.json"), emit: gated

    script:
    def limits = params.coverage_limits[meta.assembly_target]
    """
    coverage_gate.py \\
        --reads ${reads} \\
        --sample-id ${meta.sample_id} \\
        --assembly-target ${meta.assembly_target} \\
        --ref-dir ${organelle_refs} \\
        --nominal-size ${limits.nominal_size} \\
        --min-cov ${limits.min_cov} \\
        --max-cov ${limits.max_cov} \\
        --seed ${params.seqtk_seed} \\
        --threads ${task.cpus} \\
        --out-fastq ${meta.sample_id}.gated.fastq.gz
    """

    stub:
    """
    touch ${meta.sample_id}.gated.fastq.gz
    echo '{"status": "ok", "estimated_cov": 150.0}' > sample_status.json
    echo '{"pre_subsample_cov": 150.0, "post_subsample_cov": 150.0}' > coverage.json
    """
}
