// Stage 4 — post-filter read quality metrics for pre/post comparison.
// Tool: NanoPlot. See spec §2 stage 4.

process NANOPLOT_CLEAN {
    tag          "${meta.sample_id}"
    label        'process_low'
    publishDir   "${params.outdir}/${meta.sample_id}/qc/nanoplot_clean",
                 mode: 'copy', enabled: params.publish_intermediates

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("nanoplot_clean/"), emit: reports

    stub:
    """
    mkdir -p nanoplot_clean
    touch nanoplot_clean/NanoPlot-report.html nanoplot_clean/NanoStats.txt
    """

    script:
    """
    NanoPlot \\
        --threads ${task.cpus} \\
        --fastq ${reads} \\
        --outdir nanoplot_clean \\
        --tsv_stats \\
        --no_static \\
        --format png \\
        --title "${meta.sample_id} — clean"
    """
}
