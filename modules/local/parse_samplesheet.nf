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
    // STUB-02 exercises the DEDUPLICATE branch.
    def dedup = meta.sample_id == 'STUB-02'
    def duplicates = dedup
        ? ',"dedup_required":true,"duplicates":{"identical":2,"conflicting":0,"distinct_ids":2,"examples":[]}'
        : ''
    """
    touch ${meta.sample_id}.reads.fastq.gz
    echo '{"\$schema":"wf5/intake/v2","status":"valid","files":[],"total_records":0${duplicates}}' \\
        > intake.json
    """
}
