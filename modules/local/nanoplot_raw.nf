// Stage 1 — baseline read quality metrics before any filtering.
// Tool: NanoPlot. See plan.md §2 stage 1.

process NANOPLOT_RAW {
    tag          "${meta.sample_id}"
    label        'process_low'
    publishDir   "${params.outdir}/${meta.sample_id}/qc/nanoplot_raw",
                 mode: 'copy', enabled: params.publish_intermediates

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("nanoplot_raw/"), emit: reports

    stub:
    """
    mkdir -p nanoplot_raw
    touch nanoplot_raw/NanoPlot-report.html
    touch nanoplot_raw/NanoStats.txt
    """

    script:
    """
    # STUB — real implementation in P1
    mkdir -p nanoplot_raw
    touch nanoplot_raw/NanoPlot-report.html
    touch nanoplot_raw/NanoStats.txt
    """
}
