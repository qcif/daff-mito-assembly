// Stage 11 — BLAST identity sanity check against target organelle RefSeq DB.
// Tool: blastn. See plan.md §2 stage 11, §4.2.
// blast_db staged via ${file(params.blast_db)} in script block (workflow-wide ref).

process BLAST_VALIDATE {
    tag          "${meta.sample_id}"
    label        'process_medium'
    publishDir   "${params.outdir}/${meta.sample_id}/blast_validate",
                 mode: 'copy', enabled: params.publish_intermediates

    input:
    tuple val(meta), path(target_fasta), path(secondaries)

    output:
    tuple val(meta), path(target_fasta), path("${meta.sample_id}.blast.tsv"), emit: validated

    stub:
    """
    touch ${meta.sample_id}.blast.tsv
    """

    script:
    """
    # STUB — real implementation in P3
    # blastn -query ${target_fasta} -db ${file(params.blast_db)} \\
    #     -outfmt '6 qaccver saccver pident length qcovs evalue bitscore stitle' \\
    #     -max_target_seqs 5 -num_threads ${task.cpus} \\
    #     -out ${meta.sample_id}.blast.tsv
    touch ${meta.sample_id}.blast.tsv
    """
}
