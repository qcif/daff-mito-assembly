// C1 per-sample stage — concatenates pipe-delimited reads into one FASTQ.
// Single-file rows are re-emitted under the canonical filename for
// downstream tag consistency. stageAs 'input/*' avoids basename collisions
// when two reads files share a name from different subdirs.

process PARSE_SAMPLESHEET {
    tag        "${meta.sample_id}"
    label      'process_low'
    publishDir "${params.outdir}/${meta.sample_id}/qc", mode: 'copy',
               enabled: params.publish_intermediates

    // stageAs 'input*/*' puts each file in a numbered subdir (input1/, input2/…)
    // so identical basenames from different sources don't collide on staging.
    input:
    tuple val(meta), path(reads, stageAs: 'input*/*')

    output:
    tuple val(meta), path("${meta.sample_id}.reads.fastq.gz"), emit: reads

    script:
    // Byte-concatenation of gzip streams is a valid multi-member gzip file;
    // zcat and standard readers decompress them transparently.
    """
    cat input*/* > ${meta.sample_id}.reads.fastq.gz
    """

    stub:
    """
    touch ${meta.sample_id}.reads.fastq.gz
    """
}
