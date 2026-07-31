// Stage 9 — assembly graph visualisation (diagnostic).
// Tool: BandageNG. See plan.md §2 stage 9.

process BANDAGE_NG {
    tag          "${meta.sample_id}"
    label        'process_low'
    publishDir   "${params.outdir}/${meta.sample_id}/assembly",
                 mode: 'copy', enabled: params.publish_intermediates

    input:
    tuple val(meta), path(assembly), path(gfa), path(info)

    output:
    tuple val(meta), path(assembly), path(gfa), path(info),
          path("${meta.sample_id}.graph.png"), emit: assembly

    script:
    """
    BandageNG image ${gfa} ${meta.sample_id}.graph.png \\
        --height 600 \\
        --width 800
    """

    stub:
    """
    touch ${meta.sample_id}.graph.png
    """
}
