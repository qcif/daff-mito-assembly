// Stage 7 — de novo organelle assembly.
// Tool: Flye (--meta --nano-hq). See plan.md §2 stage 7, §3.4, §3.5.

process METAFLYE {
    tag          "${meta.sample_id}"
    label        'process_high'
    publishDir   "${params.outdir}/${meta.sample_id}/assembly",
                 mode: 'copy', enabled: params.publish_intermediates

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta),
          path("assembly.fasta"),
          path("assembly_graph.gfa"),
          path("assembly_info.txt"), emit: assembly

    stub:
    """
    touch assembly.fasta assembly_graph.gfa assembly_info.txt
    """

    script:
    // genome_size hint and --asm-coverage set per meta.assembly_target in P2
    """
    # STUB — real implementation in P2
    touch assembly.fasta assembly_graph.gfa assembly_info.txt
    """
}
