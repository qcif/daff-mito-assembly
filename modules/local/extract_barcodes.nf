// Stage 13 — barcode locus extraction + ORF validation (C5).
// Subsets MINIPROT_CDS's cds.gff to the loci named in assets/loci.json
// (params.locus_panel) and validates each: length, protein-identity
// floor, ORF/internal-stop check under the target-appropriate NCBI
// genetic code (clade trial on animal_mt). Never re-aligns or
// re-derives coordinates — see spec §2 stage 13, §2.2 C5.

process EXTRACT_BARCODES {
    tag          "${meta.sample_id}"
    label        'process_low'
    publishDir   "${params.outdir}/${meta.sample_id}/barcodes",
                 mode: 'copy', enabled: params.publish_intermediates

    input:
    tuple val(meta), path(target_fasta), path(cds_gff), path(genetic_code_json)

    output:
    tuple val(meta),
          path("barcodes.fasta"),
          path("${meta.sample_id}.coords.gff"),
          path("${meta.sample_id}.validation.tsv"), emit: barcodes

    script:
    def codes = params.genetic_code_tables[meta.assembly_target].join(',')
    """
    validate_barcodes.py \\
        --cds-gff ${cds_gff} \\
        --target-fasta ${target_fasta} \\
        --assembly-target ${meta.assembly_target} \\
        --locus-panel ${file(params.locus_panel)} \\
        --genetic-codes ${codes} \\
        --min-identity 60 \\
        --out-fasta barcodes.fasta \\
        --out-coords ${meta.sample_id}.coords.gff \\
        --out-tsv ${meta.sample_id}.validation.tsv
    """

    stub:
    """
    touch barcodes.fasta ${meta.sample_id}.coords.gff ${meta.sample_id}.validation.tsv
    """
}
