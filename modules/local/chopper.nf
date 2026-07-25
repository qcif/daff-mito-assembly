// Stage 2 — length/quality filter.
// Tool: chopper. See plan.md §2 stage 2.

process CHOPPER {
    tag          "${meta.sample_id}"
    label        'process_low'
    publishDir   "${params.outdir}/${meta.sample_id}/qc/chopper",
                 mode: 'copy', enabled: params.publish_intermediates

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.sample_id}.chopper.fastq.gz"), emit: reads

    stub:
    """
    touch ${meta.sample_id}.chopper.fastq.gz
    """

    script:
    """
    # STUB — real implementation in P1
    # chopper -q ${params.min_mean_q} -l ${params.min_read_length} \
    #     < <(zcat ${reads}) | gzip > ${meta.sample_id}.chopper.fastq.gz
    touch ${meta.sample_id}.chopper.fastq.gz
    """
}
