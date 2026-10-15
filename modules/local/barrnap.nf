// Stage 13a — tRNA/rRNA calling for plant_pt with barrnap (task 65).
// Runs only for targets whose params.annotate entry names
// non_cds_tool 'barrnap'; ANNOTATE (C8) merges the GFF.
//
// A barrnap failure is data, not a pipeline abort (principle 8): the
// task always exits 0 and hands its real exit code to ANNOTATE in
// barrnap_exit.txt. Peak memory is ~1.6 GB measured on the fixture
// (Rfam CM load); the 6 GB limit lives in conf/base.config. An empty target_fasta (BIN_TARGET withheld a
// substitution) gives an empty GFF and exit code 0.

process BARRNAP {
    tag          "${meta.sample_id}"
    label        'process_low'
    publishDir   "${params.outdir}/${meta.sample_id}/annotation",
                 mode: 'copy', enabled: params.publish_intermediates

    input:
    tuple val(meta), path(target_fasta)

    output:
    tuple val(meta), path("barrnap.gff"), path("barrnap_exit.txt"), emit: result
    path "barrnap.log", emit: log

    script:
    def cfg = params.annotate[meta.assembly_target]
    """
    if [[ -s ${target_fasta} ]]; then
        set +e
        barrnap --kingdom ${cfg.kingdom} --rrna --trna --threads ${task.cpus} \\
            ${target_fasta} > barrnap.gff 2> barrnap.log
        BARRNAP_EXIT=\$?
        set -e
    else
        : > barrnap.gff
        echo "empty target fasta: barrnap not run" > barrnap.log
        BARRNAP_EXIT=0
    fi
    echo "--- barrnap exit code: \${BARRNAP_EXIT}" >> barrnap.log
    echo \${BARRNAP_EXIT} > barrnap_exit.txt
    """

    stub:
    """
    touch barrnap.gff
    echo "stub" > barrnap.log
    echo 0 > barrnap_exit.txt
    """
}
