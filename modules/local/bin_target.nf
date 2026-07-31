// Stage 10 — per-contig binning: coverage spike ∩ ref identity ∩ ORF integrity.
// C3 custom logic — see spec §2.2, §3.3.
// Plant-cp canonicalisation (C4) and its path1/path2 isoform outputs
// deferred to a follow-up task; the commented emit block below stays
// commented until then.
// Animal-mt: end-overlap circularity check recorded in bin_metadata.json.

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
    path("bin_metadata.json"), emit: metadata
    // TODO(task-20 C4): uncomment when plastid_canonicalise.py lands.
    // path "plastid_isoforms/", optional: true, emit: isoforms

    script:
    def codes = params.genetic_code_tables[meta.assembly_target].join(',')
    """
    bin_target.py \\
        --assembly ${assembly} \\
        --assembly-info ${info} \\
        --organelle-ref ${organelle_refs}/${meta.assembly_target}.mmi \\
        --sample-id ${meta.sample_id} \\
        --assembly-target ${meta.assembly_target} \\
        --genetic-codes ${codes} \\
        --out-target target.fasta \\
        --out-secondaries secondaries.tsv \\
        --out-metadata bin_metadata.json
    """

    stub:
    """
    touch target.fasta secondaries.tsv bin_metadata.json
    """
}
