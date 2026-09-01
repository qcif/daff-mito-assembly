// Stage 13a cont'd — confidence scoring, step 1/3: translate every CDS
// feature in the merged annotation GFF (task 40 §2). Runs in the
// wf5-scripts image (python3.12 + biopython) — the pinned BLAST
// biocontainer the next step uses carries no Python at all, so
// translation has to happen in a separate process (task 40 §4).

process TRANSLATE_ANNOTATION_CDS {
    tag          "${meta.sample_id}"
    label        'process_low'
    publishDir   "${params.outdir}/${meta.sample_id}/annotation/scoring",
                 mode: 'copy', enabled: params.publish_intermediates

    input:
    tuple val(meta), path(annotation_gff), path(annotation_summary),
          path(target_fasta)

    output:
    tuple val(meta), path("queries"), emit: queries

    script:
    """
    mkdir -p queries
    translate_annotation_cds.py \\
        --merged-gff ${annotation_gff} \\
        --target-fasta ${target_fasta} \\
        --annotation-summary ${annotation_summary} \\
        --out-dir queries
    """

    stub:
    """
    mkdir -p queries
    """
}
