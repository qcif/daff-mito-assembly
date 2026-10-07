// Stage 10 — C12 (part 1/2): sentinel colour allocation ahead of
// BandageNG's render. Own process, not folded into BANDAGE_NG, because
// the BandageNG biocontainer has no Python (task 47 §2/§3) — the same
// constraint that forced C9 (SELECT_GENETIC_CODE) into its own step
// (rule 14). Shares EXTRACT_BARCODES's wf5/scripts image.

process ALLOCATE_GRAPH_SENTINELS {
    tag          "${meta.sample_id}"
    label        'process_low'

    input:
    tuple val(meta), path(assembly), path(gfa), path(info)

    output:
    tuple val(meta), path(assembly), path(gfa), path(info),
          path("${meta.sample_id}.sentinels.csv"), emit: sentinels

    script:
    """
    annotate_graph_svg.py allocate \\
        --gfa ${gfa} \\
        --out ${meta.sample_id}.sentinels.csv
    """

    stub:
    """
    touch ${meta.sample_id}.sentinels.csv
    """
}
