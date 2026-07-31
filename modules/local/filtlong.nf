// Stage 3 — identity-weighted top-quality read selection.
// Tool: Filtlong. See spec §2 stage 3.

process FILTLONG {
    tag          "${meta.sample_id}"
    label        'process_low'
    publishDir   "${params.outdir}/${meta.sample_id}/qc/filtlong",
                 mode: 'copy', enabled: params.publish_intermediates

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.sample_id}.filtlong.fastq.gz"), emit: reads

    script:
    """
    filtlong \\
        --min_length ${params.min_read_length} \\
        --keep_percent ${params.filtlong_keep_percent} \\
        ${reads} \\
        | gzip > ${meta.sample_id}.filtlong.fastq.gz
    """

    stub:
    """
    touch ${meta.sample_id}.filtlong.fastq.gz
    """
}
