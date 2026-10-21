// Stage 0b — remove repeated read IDs (conditional on intake.json's
// dedup_required). Only identical repeats reach here: conflicting ones
// make the sample invalid at intake. `seqkit rmdup` without -n keys on the
// ID (first header token), the rule intake_reads.py and Filtlong use.
// The removed count must equal intake's identical count, else detection
// and repair disagree and the task fails.

process DEDUPLICATE {
    tag        "${meta.sample_id}"
    label      'process_low'
    publishDir "${params.outdir}/${meta.sample_id}/qc/dedup",
               mode: 'copy', enabled: params.publish_intermediates

    input:
    tuple val(meta), path(reads), val(expected_removed)

    output:
    tuple val(meta), path("${meta.sample_id}.dedup.fastq.gz"), emit: reads
    tuple val(meta), path("dedup_stats.json"), emit: stats

    script:
    """
    count() { seqkit stats -T -j ${task.cpus} "\$1" | awk 'NR==2{print \$4}'; }
    n_in=\$(count ${reads})
    seqkit rmdup -j ${task.cpus} -o ${meta.sample_id}.dedup.fastq.gz ${reads}
    n_out=\$(count ${meta.sample_id}.dedup.fastq.gz)
    removed=\$((n_in - n_out))
    printf '{"records_in": %s, "records_out": %s, "removed": %s, "expected_removed": %s, "seqkit_version": "%s"}\\n' \\
        "\$n_in" "\$n_out" "\$removed" ${expected_removed} "\$(seqkit version)" > dedup_stats.json
    [ "\$removed" -eq ${expected_removed} ] \\
        || { echo "removed \$removed records, intake counted ${expected_removed}" >&2; exit 1; }
    """

    stub:
    """
    touch ${meta.sample_id}.dedup.fastq.gz
    echo '{"records_in": 0, "records_out": 0, "removed": 0, "expected_removed": ${expected_removed}, "seqkit_version": "seqkit stub"}' \\
        > dedup_stats.json
    """
}
