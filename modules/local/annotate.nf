// Stage 13a — CDS + non-CDS annotation merge (C8).
// CDS features come from cds.gff (MINIPROT_CDS), always. Non-CDS
// (tRNA/rRNA) features come from a target-appropriate specialist
// annotator where one exists — today MITOS2 for animal_mt. See spec
// §2 stage 13a, §2.2 C8, task 31.

process ANNOTATE {
    tag          "${meta.sample_id}"
    label        'process_medium'
    publishDir   "${params.outdir}/${meta.sample_id}/annotation",
                 mode: 'copy', enabled: params.publish_intermediates

    input:
    tuple val(meta), path(target_fasta), path(cds_gff)
    path annotate_refs

    output:
    tuple val(meta),
          path("${meta.sample_id}.gff"),
          path("annotation_summary.json"), emit: annotation

    script:
    def cfg = params.annotate[meta.assembly_target]
    // Config proxy for C5's per-run clade-trial choice — ANNOTATE has
    // no dependency on EXTRACT_BARCODES, so it cannot see what table
    // C5 actually picked (task 31 §3). Last entry of the configured
    // trial order (params.genetic_code_tables) is the representative
    // value recorded for comparison.
    def barcodeCode = params.genetic_code_tables[meta.assembly_target][-1]
    """
    if [[ "${cfg.non_cds_tool}" == "mitos2" && -s ${cds_gff} ]]; then
        mkdir -p mitos_out
        runmitos --input ${target_fasta} --outdir mitos_out \\
            --code ${cfg.genetic_code} \\
            --refdir ${annotate_refs} --refseqver mitos/${cfg.refseqver} \\
            2> mitos.stderr || true
        NONCDS="--mitos-dir mitos_out --genetic-code-annotate ${cfg.genetic_code} --reference-data ${cfg.refseqver}"
    else
        NONCDS=""
    fi

    annotate_summary.py \\
        --cds-gff ${cds_gff} \\
        --sample-id ${meta.sample_id} \\
        --assembly-target ${meta.assembly_target} \\
        --gene-sets ${file(params.gene_sets)} \\
        --genetic-code-barcodes ${barcodeCode} \\
        --out-gff ${meta.sample_id}.gff \\
        --out-summary annotation_summary.json \\
        \$NONCDS
    """

    stub:
    """
    touch ${meta.sample_id}.gff annotation_summary.json
    """
}
