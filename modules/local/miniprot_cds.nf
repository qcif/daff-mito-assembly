// Stage 12 — one broad miniprot pass over the comprehensive organellar
// protein panel. Produces every protein-coding feature found on the
// binned target contig(s); feeds both EXTRACT_BARCODES (stage 13) and
// ANNOTATE (stage 13a, task 31). See spec §2 stage 12, §4.3.
//
// miniprot's -T selects the NCBI translation table it conceptually
// translates the target under; it defaults to table 1 and every
// organelle but plant_mt needs a different one (task 38 §1). This
// process only runs miniprot — once per configured table — and never
// decides between candidates itself: the miniprot biocontainer carries
// no Python, so the animal_mt clade-trial selection is a separate
// process, SELECT_GENETIC_CODE (C9), in main.nf (task 38 §2/§3).
//
// plant_pt/plant_mt are single-valued: `resolved.cds.gff` /
// `resolved.genetic_code.json` are written directly here (plain
// shell, no selector needed) and SELECT_GENETIC_CODE is skipped for
// those samples entirely — no new failure surface for those branches.

process MINIPROT_CDS {
    tag          "${meta.sample_id}"
    label        'process_medium'
    publishDir   "${params.outdir}/${meta.sample_id}/cds",
                 mode: 'copy', enabled: params.publish_intermediates

    input:
    tuple val(meta), path(target_fasta), path(blast_tsv)
    path protein_panel

    output:
    tuple val(meta), path(target_fasta),
          path("candidates/cds.*.gff"), emit: candidates
    tuple val(meta), path("${meta.sample_id}.cds.gff"),
          path("genetic_code.json"), optional: true, emit: resolved

    script:
    def codes = params.genetic_code_tables[meta.assembly_target].join(',')
    """
    mkdir -p candidates
    for T in \$(echo ${codes} | tr ',' ' '); do
        if [[ -s ${target_fasta} ]]; then
            cat ${protein_panel}/${meta.assembly_target}/*.faa > panel.faa
            miniprot --gff -T \${T} -t ${task.cpus} ${target_fasta} panel.faa \\
                > candidates/cds.\${T}.gff
        else
            # Empty target.fasta from BIN_TARGET / withheld C4 substitution
            # (task 24 §3.2) → empty candidate gff(s); not an error.
            : > candidates/cds.\${T}.gff
        fi
    done

    if [[ "${codes}" != *,* ]]; then
        cp candidates/cds.${codes}.gff ${meta.sample_id}.cds.gff
        echo '{"\$schema": "wf5/genetic-code-selection/v1", "selected_table": ${codes}, "criterion": "single configured table for ${meta.assembly_target} — no trial", "candidates": [{"table": ${codes}, "gff": "candidates/cds.${codes}.gff"}]}' > genetic_code.json
    fi
    # animal_mt clade trial: leave ${meta.sample_id}.cds.gff /
    # genetic_code.json unwritten — SELECT_GENETIC_CODE (C9) publishes
    # them into this same directory once it has picked a winner, and
    # `optional: true` above keeps this process from publishing empty
    # placeholders that would race with C9's real output.
    """

    stub:
    """
    mkdir -p candidates
    touch candidates/cds.1.gff ${meta.sample_id}.cds.gff genetic_code.json
    """
}
