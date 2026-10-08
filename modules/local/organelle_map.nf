// Stage 14 — annotated organelle map (diagnostic visualisation).
// In-house Jinja+SVG renderer — see plan.md §8 item 4 and task
// 44_organelle_map.md §1 for the tool-choice decision.
// Input: GFF3 + annotation_summary.json from ANNOTATION_SCORING (task 40),
// plus bin_metadata.json (C3/C4) for genome length and, on the plant_pt
// canonical branch, the LSC/IR/SSC segment lengths needed to remap
// features onto the plastid's second isoform (spec §3.6 step 5). Also
// reads params.locus_panel to label barcode-panel loci by name, the
// same match EXTRACT_BARCODES uses (task 54 §5.4).
// Output: inline-ready SVG.

process ORGANELLE_MAP {
    tag          "${meta.sample_id}"
    label        'process_low'
    publishDir   "${params.outdir}/${meta.sample_id}/annotation",
                 mode: 'copy', enabled: params.publish_intermediates

    input:
    tuple val(meta), path(gff), path(annotation_summary), path(bin_metadata)

    output:
    tuple val(meta), path("${meta.sample_id}.map.svg"), emit: map

    script:
    """
    organelle_map.py \\
        --gff ${gff} \\
        --annotation-summary ${annotation_summary} \\
        --bin-metadata ${bin_metadata} \\
        --locus-panel ${file(params.locus_panel)} \\
        --assembly-target ${meta.assembly_target} \\
        --sample-id ${meta.sample_id} \\
        --out ${meta.sample_id}.map.svg
    """

    stub:
    """
    touch ${meta.sample_id}.map.svg
    """
}
