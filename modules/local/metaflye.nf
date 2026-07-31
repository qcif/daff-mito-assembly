// Stage 7 — de novo organelle assembly.
// Tool: Flye (--meta --nano-hq). See spec §2 stage 7, §3.4, §3.5.

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

    script:
    // --asm-coverage is incompatible with --meta in Flye 2.9.6 (spec §3.5 deferred).
    def genome_size = params.assembly_size_hints[meta.assembly_target]
    """
    flye --meta --nano-hq ${reads} \\
         --genome-size ${genome_size} \\
         --threads ${task.cpus} \\
         --out-dir flye_out
    mv flye_out/assembly.fasta       assembly.fasta
    mv flye_out/assembly_graph.gfa   assembly_graph.gfa
    mv flye_out/assembly_info.txt    assembly_info.txt
    """

    stub:
    """
    touch assembly.fasta assembly_graph.gfa assembly_info.txt
    """
}
