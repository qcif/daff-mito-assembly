// Stage 12 — full organelle annotation.
// Target-dispatched: animal_mt → MITOS2; plant_pt / plant_mt → TBD
// (see spec §8). See spec §2 stage 12, §3.2.

process ANNOTATE {
    tag          "${meta.sample_id}"
    label        'process_medium'
    publishDir   "${params.outdir}/${meta.sample_id}/annotation",
                 mode: 'copy', enabled: params.publish_intermediates

    input:
    tuple val(meta), path(target_fasta), path(blast_tsv)

    output:
    tuple val(meta), path("${meta.sample_id}.gff"), path("${meta.sample_id}.gbk"), emit: annotation

    stub:
    """
    touch ${meta.sample_id}.gff ${meta.sample_id}.gbk
    """

    script:
    // Tool dispatched by meta.assembly_target in P4
    """
    # STUB — real implementation in P4
    # animal_mt → MITOS2; plant_pt / plant_mt → TBD (see spec §8)
    touch ${meta.sample_id}.gff ${meta.sample_id}.gbk
    """
}
