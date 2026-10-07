// Stage 10 — assembly graph visualisation (diagnostic). Tool: BandageNG.
// See spec/02-stages.md stage 10.
//
// Renders SVG, not PNG, so the report can bind per-node tooltips
// (task 47). BandageNG's own SVG carries no node identity (task 47
// §2), so this render is sandwiched between ALLOCATE_GRAPH_SENTINELS
// (writes the `name,color` CSV this process's `--color` consumes) and
// ANNOTATE_GRAPH_SVG (rewrites the sentinel fills back to a neutral
// colour and attaches the real node identity) — the round-trip is
// spread across three processes rather than one because the
// BandageNG biocontainer has no Python either side of the render
// (rule 14, task 47 §4.2 route (b)).

process BANDAGE_NG {
    tag          "${meta.sample_id}"
    label        'process_low'

    input:
    tuple val(meta), path(assembly), path(gfa), path(info), path(sentinels_csv)

    output:
    tuple val(meta), path(assembly), path(gfa), path(info),
          path(sentinels_csv), path("${meta.sample_id}.graph.raw.svg"),
          emit: rendered

    script:
    """
    BandageNG image ${gfa} ${meta.sample_id}.graph.raw.svg \\
        --height 600 \\
        --width 800 \\
        --color ${sentinels_csv}
    """

    stub:
    """
    touch ${meta.sample_id}.graph.raw.svg
    """
}
