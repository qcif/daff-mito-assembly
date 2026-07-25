// Stage 5 — positive recruitment against kingdom organelle reference panel.
// Tools: minimap2 + samtools + seqtk. See plan.md §2 stage 5.
// kingdom_refs staged via Channel.value(file(params.kingdom_refs)) in main.nf.

process RECRUIT {
    tag          "${meta.sample_id}"
    label        'process_medium'
    publishDir   "${params.outdir}/${meta.sample_id}/recruit",
                 mode: 'copy', enabled: params.publish_intermediates

    input:
    tuple val(meta), path(reads)
    path kingdom_refs

    output:
    tuple val(meta), path("${meta.sample_id}.recruited.fastq.gz"), emit: reads

    stub:
    """
    touch ${meta.sample_id}.recruited.fastq.gz
    """

    script:
    """
    # STUB — real implementation in P1
    # minimap2 -ax map-ont -t ${task.cpus} ${kingdom_refs} ${reads} \\
    #     | samtools view -b -F 4 -q 1 -@ ${task.cpus} \\
    #     | samtools fastq -@ ${task.cpus} \\
    #     | gzip > ${meta.sample_id}.recruited.fastq.gz
    touch ${meta.sample_id}.recruited.fastq.gz
    """
}
