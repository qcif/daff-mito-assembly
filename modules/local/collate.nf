// Stage 15 — per-sample bundle aggregation + report render.
// C6 custom logic — see plan.md §2.2, §6a.5.
// Handles both ok and soft-failed (low_coverage) samples.
// Invokes wf-report-boilerplate/report.py to produce report.html.

process COLLATE {
    tag          "${meta.sample_id}"
    label        'process_low'
    publishDir   "${params.outdir}/${meta.sample_id}",
                 mode: 'copy'

    input:
    tuple val(meta),
          path(status_json),
          path(coverage_json),
          path(nanoplot_raw),
          path(nanoplot_clean),
          path(target_fasta),
          path(blast_tsv),
          path(barcodes_fasta),
          path(coords_gff),
          path(validation_tsv),
          path(annotation_gff),
          path(annotation_gbk),
          path(organelle_map_svg),
          path(graph_png)

    output:
    // Paths marked optional: true are absent for soft-failed (low-coverage) samples;
    // real logic (P4 bin/collate.py) creates them conditionally.
    tuple val(meta),
          path("organelle_assembly.fasta",   optional: true),
          path("organelle_annotation.gff",   optional: true),
          path("barcodes.fasta",             optional: true),
          path("metadata.json"),
          path("report.html"), emit: bundle

    stub:
    """
    touch organelle_assembly.fasta organelle_annotation.gff barcodes.fasta
    touch metadata.json report.html
    """

    script:
    """
    # STUB — real implementation in P4 (bin/collate.py)
    # Dispatches on status_json to emit full or minimal (low-coverage) bundle.
    touch organelle_assembly.fasta organelle_annotation.gff barcodes.fasta
    touch metadata.json report.html
    """
}
