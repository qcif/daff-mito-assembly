// Stage 9 — C12 (part 2/2): node identity for BandageNG's SVG.
// Post-processes BANDAGE_NG's sentinel-coloured raw SVG into the
// published, hover-interactive `<sample_id>.graph.svg` (task 47 §3) —
// see bin/annotate_graph_svg.py's module docstring for the round-trip
// design. Own process for the same reason as
// ALLOCATE_GRAPH_SENTINELS: the BandageNG biocontainer has no Python.
// Shares EXTRACT_BARCODES's wf5/scripts image.

process ANNOTATE_GRAPH_SVG {
    tag          "${meta.sample_id}"
    label        'process_low'
    publishDir   "${params.outdir}/${meta.sample_id}/assembly",
                 mode: 'copy', enabled: params.publish_intermediates

    input:
    tuple val(meta), path(assembly), path(gfa), path(info),
          path(sentinels_csv), path(raw_svg)

    output:
    tuple val(meta), path(assembly), path(gfa), path(info),
          path("${meta.sample_id}.graph.svg"), emit: assembly

    script:
    """
    annotate_graph_svg.py annotate \\
        --svg ${raw_svg} \\
        --sentinels ${sentinels_csv} \\
        --assembly-info ${info} \\
        --sample-id ${meta.sample_id} \\
        --out ${meta.sample_id}.graph.svg
    """

    stub:
    """
    touch ${meta.sample_id}.graph.svg
    """
}
