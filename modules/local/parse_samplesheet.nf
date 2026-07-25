// C1 custom logic — see plan.md §2.2
// Full implementation: bin/parse_samplesheet.py (P1)
// Validates the CSV, resolves pipe-delimited reads against --data_dir,
// concatenates multi-file samples, emits one tuple(meta, reads) per row.

process PARSE_SAMPLESHEET {
    tag   "${meta.sample_id}"
    label 'process_low'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.sample_id}.reads.fastq.gz"), emit: reads

    stub:
    """
    touch ${meta.sample_id}.reads.fastq.gz
    """

    script:
    """
    # STUB — real implementation in P1 (bin/parse_samplesheet.py)
    touch ${meta.sample_id}.reads.fastq.gz
    """
}
