// Stage 13a cont'd — confidence scoring, step 2/3: blastp each gene's
// translated queries against that same gene's protein-panel file —
// the panel miniprot itself aligned to (task 40 §2/§4). One blastp
// call per gene present in this sample's queries/, not one per
// feature: cheap (a handful of genes, ~10 panel sequences each), and
// keeps every hit scoped to its own gene family rather than the whole
// panel, so a scored hit is never contaminated by a homologous but
// wrong gene. Reuses BLAST_VALIDATE's already-pinned BLAST 2.17.0
// image — no new container (task 40 §4, option 1).

process SCORE_ANNOTATION_BLASTP {
    tag          "${meta.sample_id}"
    label        'process_low'
    publishDir   "${params.outdir}/${meta.sample_id}/annotation/scoring",
                 mode: 'copy', enabled: params.publish_intermediates

    input:
    tuple val(meta), path(queries_dir)
    path protein_panel

    output:
    tuple val(meta), path("blast.tsv"), emit: scored

    script:
    """
    : > blast.tsv
    shopt -s nullglob
    for query in ${queries_dir}/*.faa; do
        gene=\$(basename "\${query}" .faa)
        panel="${protein_panel}/${meta.assembly_target}/\${gene}.faa"
        if [[ -s "\${query}" && -s "\${panel}" ]]; then
            blastp -query "\${query}" -subject "\${panel}" \\
                -outfmt '6 qseqid sseqid pident qcovhsp evalue bitscore' \\
                -max_target_seqs 5 -num_threads ${task.cpus} \\
                >> blast.tsv
        fi
    done
    """

    stub:
    """
    touch blast.tsv
    """
}
