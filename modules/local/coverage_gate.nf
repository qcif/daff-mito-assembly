// Stage 6 — coverage estimation + subsample/soft-fail gate.
// C2 custom logic — see plan.md §2.2 and §2.1.
// Always exits 0; gate decision is written to sample_status.json (data, not error).
// errorStrategy 'ignore' guards against unexpected seqkit/seqtk crashes only.

process COVERAGE_GATE {
    tag            "${meta.sample_id}"
    label          'process_low'
    errorStrategy  'ignore'
    publishDir     "${params.outdir}/${meta.sample_id}/coverage_gate",
                   mode: 'copy', enabled: params.publish_intermediates

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
    """
    # STUB — real implementation in P1 (bin/coverage_gate.py)
    # Writes sample_status.json with status: ok | low_coverage
    # Subsamples with seqtk if estimated_cov > MAX
    touch ${meta.sample_id}.gated.fastq.gz
    echo '{"status": "ok", "estimated_cov": 150.0, "min_required": 30, "max_allowed": 500}' > sample_status.json
    echo '{"pre_subsample_cov": 150.0, "post_subsample_cov": 150.0, "seed": ${params.seqtk_seed}}' > coverage.json
    """
}
