// Stage 15 — per-sample bundle aggregation + report render.
// C6 custom logic — see plan.md §2.2, §6a.5.
// Handles all three gate statuses (spec §2.1.3, task 35): `fail` gets
// the minimal bundle (no assembly exists); `low_coverage` and `ok` both
// get the full bundle, with `low_coverage` carrying the warning into
// metadata.json. See tasks/todo.md's COLLATE carry-forward for the
// real (P4) implementation of this dispatch.
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
          path(annotation_summary),
          path(organelle_map_svg),
          path(graph_png)

    output:
    // Paths marked optional: true are absent for gate-failed (`fail` status)
    // samples; real logic (P4 bin/collate.py) creates them conditionally.
    tuple val(meta),
          path("organelle_assembly.fasta",   optional: true),
          path("organelle_annotation.gff",   optional: true),
          path("barcodes.fasta",             optional: true),
          path("metadata.json"),
          path("report.html"), emit: bundle

    script:
    """
    # STUB — real implementation in P4 (bin/collate.py)
    # Dispatches on status_json to emit full (ok/low_coverage) or
    # minimal (fail) bundle.
    touch organelle_assembly.fasta organelle_annotation.gff barcodes.fasta
    touch metadata.json report.html
    """

    stub:
    """
    touch organelle_assembly.fasta organelle_annotation.gff barcodes.fasta
    touch metadata.json report.html
    """
}
