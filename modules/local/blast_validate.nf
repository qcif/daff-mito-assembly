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

    script:
    def db_name = [
        animal_mt: 'refseq_mt_metazoa',
        plant_pt:  'refseq_pt',
        plant_mt:  'refseq_mt_viridiplantae',
    ][meta.assembly_target]
    """
    if [[ -s ${target_fasta} ]]; then
        blastn \\
            -query ${target_fasta} \\
            -db ${file(params.blast_db)}/${db_name} \\
            -outfmt '6 qaccver saccver pident length qcovs evalue bitscore stitle' \\
            -max_target_seqs 5 \\
            -num_threads ${task.cpus} \\
            -out ${meta.sample_id}.blast.tsv
    else
        # Empty target.fasta from BIN_TARGET → empty BLAST output; not an error.
        : > ${meta.sample_id}.blast.tsv
    fi
    """

    stub:
    """
    touch ${meta.sample_id}.blast.tsv
    """
}
