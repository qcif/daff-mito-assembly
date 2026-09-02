// Stage 10 — per-contig binning: merged homology to the declared panel,
// discriminated against the sibling organelle panel(s), ranked on coverage.
// C3 custom logic — see spec §2.2, §3.3, §3.7. Thresholds come from
// params.bin_target_thresholds (spec §3.7.6), not from the Python script.
// Plant-pt canonicalisation (C4) runs in-process from bin_target.py on the
// plant_pt branch; see spec/plastid-canonicalisation.md and task 20.
// Circularity is read from Flye's circ. column, with the end-overlap
// self-alignment as a fallback (spec §3.7.4).

process BIN_TARGET {
    tag          "${meta.sample_id}"
    label        'process_medium'
    publishDir   "${params.outdir}/${meta.sample_id}/bin_target",
                 mode: 'copy', enabled: params.publish_intermediates

    input:
    tuple val(meta), path(assembly), path(gfa), path(info), path(graph_png)
    path organelle_refs

    output:
    tuple val(meta), path("target.fasta"), path("secondaries.tsv"), emit: binned
    tuple val(meta), path("bin_metadata.json"), emit: metadata
    // `meta` on both of these (task 42 §2.2 defects 1/2) so COLLATE
    // (and task 44_organelle_map.md for isoforms) can `join` them
    // per-sample — a bare `path(...)` output can't be joined at all.
    tuple val(meta), path("plastid_isoforms/"), optional: true, emit: isoforms

    script:
    def codes = params.genetic_code_tables[meta.assembly_target].join(',')
    def th    = params.bin_target_thresholds[meta.assembly_target]
    """
    bin_target.py \\
        --assembly ${assembly} \\
        --assembly-info ${info} \\
        --gfa ${gfa} \\
        --ref-dir ${organelle_refs} \\
        --sample-id ${meta.sample_id} \\
        --assembly-target ${meta.assembly_target} \\
        --genetic-codes ${codes} \\
        --min-identity ${th.min_identity} \\
        --min-aligned-frac ${th.min_aligned_frac} \\
        --emit ${th.emit} \\
        --max-contigs ${th.max_contigs} \\
        --low-coverage-fraction ${th.low_coverage_fraction} \\
        --sibling-warn-fraction ${th.sibling_warn_fraction} \\
        --out-target target.fasta \\
        --out-secondaries secondaries.tsv \\
        --out-metadata bin_metadata.json
    """

    stub:
    """
    touch target.fasta secondaries.tsv bin_metadata.json
    """
}
