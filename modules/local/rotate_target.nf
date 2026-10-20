// Stage 9a — C13 (bin/rotate_target.py): rotates each circular animal_mt
// target contig to start at its anchor gene (COX1, then fallbacks) on the
// forward strand, so every run lands in the same coordinate frame.
// Targets with no configured anchors pass through unchanged. Runs in the
// wf5/scripts image because the miniprot biocontainer has no Python
// (rules 11, 12). See spec §2 stage 9a, §2.2 C13, §3.8.

process ROTATE_TARGET {
    tag          "${meta.sample_id}"
    label        'process_low'
    publishDir   "${params.outdir}/${meta.sample_id}/rotate",
                 mode: 'copy', enabled: params.publish_intermediates

    input:
    tuple val(meta), path(target_fasta), path(secondaries),
          path(bin_metadata), path(anchor_gffs), val(anchors)

    output:
    tuple val(meta), path("rotated/target.fasta"), path(secondaries), emit: binned
    tuple val(meta), path("rotation.json"), emit: rotation

    script:
    def gffArgs = anchor_gffs.collect { f ->
        "--gff ${f.name.replaceAll(/^anchors\.(\d+)\.gff$/, '$1')}:${f}"
    }.join(' ')
    """
    mkdir -p rotated
    rotate_target.py \\
        --target-fasta ${target_fasta} \\
        --bin-metadata ${bin_metadata} \\
        --anchors '${anchors.join(',')}' \\
        ${gffArgs} \\
        --out-fasta rotated/target.fasta \\
        --out-json rotation.json
    """

    stub:
    """
    mkdir -p rotated
    touch rotated/target.fasta rotation.json
    """
}
