// Stage 9a — miniprot pass over the anchor-gene panels only, once per
// configured genetic-code table, so ROTATE_TARGET (C13) can find where
// to start position 1 on circular animal_mt contigs. Runs before
// MINIPROT_CDS because rotation must precede annotation. See spec §2
// stage 9a, §3.8.

process ROTATE_ANCHORS {
    tag          "${meta.sample_id}"
    label        'process_low'
    publishDir   "${params.outdir}/${meta.sample_id}/rotate",
                 mode: 'copy', enabled: params.publish_intermediates

    input:
    tuple val(meta), path(target_fasta), path(secondaries)
    path protein_panel

    output:
    tuple val(meta), path("anchors/anchors.*.gff"), emit: anchors

    script:
    def codes   = params.genetic_code_tables[meta.assembly_target].join(' ')
    def anchors = (params.rotation_anchors[meta.assembly_target] ?: []).join(' ')
    """
    mkdir -p anchors
    : > panel.faa
    for G in ${anchors}; do
        PANEL=${protein_panel}/${meta.assembly_target}/\${G}.faa
        if [[ ! -s \${PANEL} ]]; then
            echo "ERROR: rotation anchor \${G} has no panel file \${PANEL}" >&2
            exit 1
        fi
        cat \${PANEL} >> panel.faa
    done

    for T in ${codes}; do
        if [[ -s ${target_fasta} && -s panel.faa ]]; then
            miniprot --gff -T \${T} -t ${task.cpus} ${target_fasta} panel.faa \\
                > anchors/anchors.\${T}.gff
        else
            : > anchors/anchors.\${T}.gff
        fi
    done
    """

    stub:
    """
    mkdir -p anchors
    touch anchors/anchors.1.gff
    """
}
