// Stage 14 — annotated organelle map (diagnostic visualisation).
// Tool choice deferred — see plan.md §8.
// Input: annotated GenBank from ANNOTATE. Output: inline-ready SVG.

process ORGANELLE_MAP {
    tag          "${meta.sample_id}"
    label        'process_low'
    publishDir   "${params.outdir}/${meta.sample_id}/annotation",
                 mode: 'copy', enabled: params.publish_intermediates

    input:
    tuple val(meta), path(gff), path(gbk)

    output:
    tuple val(meta), path("${meta.sample_id}.map.svg"), emit: map

    stub:
    """
    touch ${meta.sample_id}.map.svg
    """

    script:
    """
    # STUB — real implementation in P4, once tool is selected (plan.md §8)
    touch ${meta.sample_id}.map.svg
    """
}
