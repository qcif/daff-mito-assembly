// Stage 10 — C12 (part 2/2): node identity + binning-bucket colour for
// BandageNG's SVG. Post-processes BANDAGE_NG's sentinel-coloured raw
// SVG into the published, hover-interactive `<sample_id>.graph.svg`
// (task 47 §3; bucket colouring added task 52 §3.4) — see
// bin/annotate_graph_svg.py's module docstring for the round-trip
// design. Own process for the same reason as
// ALLOCATE_GRAPH_SENTINELS: the BandageNG biocontainer has no Python.
// Shares EXTRACT_BARCODES's wf5/scripts image.
//
// Joined to BIN_TARGET.out.metadata by `meta` in main.nf (task 52
// §3.2) so colouring can read bin_metadata.json — a sample whose
// BIN_TARGET process errors outright loses its graph along with
// everything else downstream, same as before this task.

process ANNOTATE_GRAPH_SVG {
    tag          "${meta.sample_id}"
    label        'process_low'
    publishDir   "${params.outdir}/${meta.sample_id}/assembly",
                 mode: 'copy', enabled: params.publish_intermediates

    input:
    tuple val(meta), path(assembly), path(gfa), path(info),
          path(sentinels_csv), path(raw_svg), path(bin_metadata_json)

    output:
    tuple val(meta), path("${meta.sample_id}.graph.svg"), emit: graph

    script:
    """
    annotate_graph_svg.py annotate \\
        --svg ${raw_svg} \\
        --sentinels ${sentinels_csv} \\
        --assembly-info ${info} \\
        --sample-id ${meta.sample_id} \\
        --bin-metadata ${bin_metadata_json} \\
        --out ${meta.sample_id}.graph.svg
    """

    stub:
    """
    touch ${meta.sample_id}.graph.svg
    """
}
