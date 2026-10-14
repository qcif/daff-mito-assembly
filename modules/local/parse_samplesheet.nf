// C1b per-sample stage — validates raw reads and concatenates them into
// one FASTQ. intake.json is always emitted; the reads output only exists
// when every input passed (bin/intake_reads.py). stageAs 'input*/*'
// avoids basename collisions when two reads files share a name from
// different subdirs.

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
    tuple val(meta), path("${meta.sample_id}.reads.fastq.gz"),
          optional: true, emit: reads
    tuple val(meta), path("intake.json"), emit: intake

    script:
    """
    intake_reads.py \\
        --reads ${reads} \\
        --out-reads ${meta.sample_id}.reads.fastq.gz \\
        --out-status intake.json
    """

    stub:
    """
    touch ${meta.sample_id}.reads.fastq.gz
    echo '{"\$schema":"wf5/intake/v1","status":"valid","files":[],"total_records":0}' \\
        > intake.json
    """
}
