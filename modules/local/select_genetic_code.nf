// Stage 12, animal_mt only — C9: picks between MINIPROT_CDS's per-table
// candidate cds.gff files (bin/select_genetic_code.py). Runs in its own
// container (rule 14): the miniprot biocontainer has no Python, so this
// custom-logic component shares EXTRACT_BARCODES's stdlib-capable
// wf5/scripts image instead. See spec §2 stage 12, §2.2 C9, task 38 §2/§3.

process SELECT_GENETIC_CODE {
    tag          "${meta.sample_id}"
    label        'process_low'
    publishDir   "${params.outdir}/${meta.sample_id}/cds",
                 mode: 'copy', enabled: params.publish_intermediates

    input:
    tuple val(meta), path(target_fasta), path(candidate_gffs)

    output:
    tuple val(meta), path("${meta.sample_id}.cds.gff"),
          path("genetic_code.json"), emit: selected

    script:
    """
    CANDIDATES=""
    for f in ${candidate_gffs}; do
        T=\$(echo \$f | sed -E 's/cds\\.([0-9]+)\\.gff/\\1/')
        CANDIDATES="\$CANDIDATES --candidate \${T}:\${f}"
    done
    select_genetic_code.py \${CANDIDATES} \\
        --out-gff ${meta.sample_id}.cds.gff \\
        --out-json genetic_code.json
    """

    stub:
    """
    touch ${meta.sample_id}.cds.gff genetic_code.json
    """
}
