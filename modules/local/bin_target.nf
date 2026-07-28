// Stage 10 — per-contig binning: coverage spike ∩ ref identity ∩ ORF integrity.
// C3 + C4 custom logic — see plan.md §2.2, §3.3, §3.6.
// Plant-cp branch: quadripartite canonicalisation → path1.fasta + path2.fasta.
// Animal branch: end-overlap circularity check.

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
    // Plant-cp isoforms (P3): emitted as a separate named output once the real
    // canonicalisation logic (bin/plastid_canonicalise.py) is wired.
    // path "plastid_isoforms/", optional: true, emit: isoforms

    stub:
    """
    touch target.fasta secondaries.tsv
    """

    script:
    """
    # STUB — real implementation in P3 (bin/bin_target.py + bin/plastid_canonicalise.py)
    touch target.fasta secondaries.tsv
    """
}
