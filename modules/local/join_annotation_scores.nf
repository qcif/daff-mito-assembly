// Stage 13a cont'd — confidence scoring, step 3/3: join blastp scores
// back onto the merged GFF + annotation_summary.json (task 40). This
// is the annotation stage's final, published output — every CDS
// feature now carries pident/qcovhsp/bitscore regardless of
// provenance, or an explicit null where nothing scored (task 40 §5:
// scores are data, never a gate — no feature is dropped here).

process JOIN_ANNOTATION_SCORES {
    tag          "${meta.sample_id}"
    label        'process_low'
    publishDir   "${params.outdir}/${meta.sample_id}/annotation",
                 mode: 'copy', enabled: params.publish_intermediates

    input:
    tuple val(meta), path(annotation_gff), path(annotation_summary),
          path(blast_tsv)

    output:
    tuple val(meta),
          path("${meta.sample_id}.gff"),
          path("annotation_summary.json"), emit: annotation

    script:
    """
    annotation_scores.py \\
        --merged-gff ${annotation_gff} \\
        --annotation-summary ${annotation_summary} \\
        --blast-tsv ${blast_tsv} \\
        --out-gff ${meta.sample_id}.gff \\
        --out-summary annotation_summary.json
    """

    stub:
    """
    touch ${meta.sample_id}.gff annotation_summary.json
    """
}
