// Stage 13 — barcode locus extraction + ORF validation.
// Tool: miniprot. C5 custom validation logic. See plan.md §2 stage 13, §3.3.
// locus_panel staged via ${file(params.locus_panel)} in script block.

process MINIPROT_EXTRACT {
    tag          "${meta.sample_id}"
    label        'process_medium'
    publishDir   "${params.outdir}/${meta.sample_id}/barcodes",
                 mode: 'copy', enabled: params.publish_intermediates

    input:
    tuple val(meta), path(target_fasta), path(blast_tsv)

    output:
    tuple val(meta),
          path("barcodes.fasta"),
          path("${meta.sample_id}.coords.gff"),
          path("${meta.sample_id}.validation.tsv"), emit: barcodes

    stub:
    """
    touch barcodes.fasta ${meta.sample_id}.coords.gff ${meta.sample_id}.validation.tsv
    """

    script:
    """
    # STUB — real implementation in P3 (bin/validate_barcodes.py)
    # miniprot -t ${task.cpus} ${target_fasta} ${file(params.locus_panel)} > aln.paf
    # bin/validate_barcodes.py --kingdom ${meta.kingdom} --paf aln.paf ...
    touch barcodes.fasta ${meta.sample_id}.coords.gff ${meta.sample_id}.validation.tsv
    """
}
