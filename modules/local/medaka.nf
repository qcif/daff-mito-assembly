// Stage 8 — optional homopolymer polish.
// Tool: medaka. Opt-in via --polish flag (default false). See plan.md §2 stage 8.

process MEDAKA {
    tag          "${meta.sample_id}"
    label        'process_high'
    publishDir   "${params.outdir}/${meta.sample_id}/assembly/medaka",
                 mode: 'copy', enabled: params.publish_intermediates

    input:
    tuple val(meta), path(assembly), path(gfa), path(info)
    path reads

    output:
    tuple val(meta),
          path("${meta.sample_id}.polished.fasta"),
          path(gfa),
          path(info), emit: assembly

    script:
    """
    # STUB — real implementation in P2
    # medaka_consensus -i ${reads} -d ${assembly} -o medaka_out -t ${task.cpus}
    # cp medaka_out/consensus.fasta ${meta.sample_id}.polished.fasta
    touch ${meta.sample_id}.polished.fasta
    """

    stub:
    """
    touch ${meta.sample_id}.polished.fasta
    """
}
