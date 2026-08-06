// Stage 12 — one broad miniprot pass over the comprehensive organellar
// protein panel. Produces every protein-coding feature found on the
// binned target contig(s); feeds both EXTRACT_BARCODES (stage 13) and
// ANNOTATE (stage 13a, task 31). See spec §2 stage 12, §4.3.

process MINIPROT_CDS {
    tag          "${meta.sample_id}"
    label        'process_medium'
    publishDir   "${params.outdir}/${meta.sample_id}/cds",
                 mode: 'copy', enabled: params.publish_intermediates

    input:
    tuple val(meta), path(target_fasta), path(blast_tsv)
    path protein_panel

    output:
    tuple val(meta), path(target_fasta),
          path("${meta.sample_id}.cds.gff"), emit: cds

    script:
    """
    if [[ -s ${target_fasta} ]]; then
        cat ${protein_panel}/${meta.assembly_target}/*.faa > panel.faa
        miniprot --gff -t ${task.cpus} ${target_fasta} panel.faa \\
            > ${meta.sample_id}.cds.gff
    else
        # Empty target.fasta from BIN_TARGET / withheld C4 substitution
        # (task 24 §3.2) → empty cds.gff; not an error.
        : > ${meta.sample_id}.cds.gff
    fi
    """

    stub:
    """
    touch ${meta.sample_id}.cds.gff
    """
}
