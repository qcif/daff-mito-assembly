// Stage 5 — positive recruitment against target organelle reference.
// Tools: minimap2 + samtools + seqtk (mulled biocontainer).
// See spec §2 stage 5 and brief §3.2–3.3.

process RECRUIT {
    tag          "${meta.sample_id}"
    label        'process_medium'
    publishDir   "${params.outdir}/${meta.sample_id}/recruit",
                 mode: 'copy', enabled: params.publish_intermediates

    input:
    tuple val(meta), path(reads)
    path organelle_refs

    output:
    tuple val(meta), path("${meta.sample_id}.recruited.fastq.gz"), emit: reads

    script:
    """
    minimap2 -ax map-ont -t ${task.cpus} \\
             ${organelle_refs}/${meta.assembly_target}.mmi ${reads} \\
        | samtools view -F 4 -q 1 -@ ${task.cpus} \\
        | cut -f1 \\
        | sort -u > ids.txt

    seqtk subseq ${reads} ids.txt \\
        | gzip > ${meta.sample_id}.recruited.fastq.gz
    """

    stub:
    """
    touch ${meta.sample_id}.recruited.fastq.gz
    """
}
