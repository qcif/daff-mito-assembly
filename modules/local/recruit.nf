// Stage 5 — positive recruitment against target organelle reference.
// Tools: minimap2 + seqtk + python3 (mulled biocontainer, shared with C2).
// See spec §2 stage 5, §9 item 1, and brief §3.2–3.3.
//
// Selection is on merged aligned extent, not mapping quality. The old
// `samtools view -F 4 -q 1` filter inverts on the broad, redundant panels
// task 28 introduced — MAPQ measures which *reference genome* a read came
// from, a question recruitment does not ask, and on `INT-PLANT-01-mt` it
// preferentially discarded mitochondrial reads (task 28 §9.5). Thresholds
// are per-target and live in params.recruit_thresholds, not in the script
// (CONSTITUTION rule 18); samtools is no longer needed since the filter
// reads PAF directly.

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
    tuple val(meta), path("recruit_stats.json"), emit: stats

    script:
    def th = params.recruit_thresholds[meta.assembly_target]
    """
    minimap2 -x map-ont -t ${task.cpus} \\
             ${organelle_refs}/${meta.assembly_target}.mmi ${reads} \\
        > hits.paf

    recruit_filter.py \\
        --paf hits.paf \\
        --min-aligned-frac ${th.min_aligned_frac} \\
        --min-aligned-bp ${th.min_aligned_bp} \\
        --stats recruit_stats.json \\
        --out ids.txt

    seqtk subseq ${reads} ids.txt \\
        | gzip > ${meta.sample_id}.recruited.fastq.gz
    """

    stub:
    """
    touch ${meta.sample_id}.recruited.fastq.gz recruit_stats.json
    """
}
