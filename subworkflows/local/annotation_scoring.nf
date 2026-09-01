// Confidence scoring subworkflow (task 40) — wraps the three-process
// translate -> blastp -> join chain behind one (meta, gff, summary) in,
// (meta, gff, summary) out interface, so main.nf reads at the stage
// level and the plumbing needed only because the pinned BLAST
// biocontainer has no Python (task 40 §4) stays internal to this file.

include { TRANSLATE_ANNOTATION_CDS } from '../../modules/local/translate_annotation_cds'
include { SCORE_ANNOTATION_BLASTP  } from '../../modules/local/score_annotation_blastp'
include { JOIN_ANNOTATION_SCORES   } from '../../modules/local/join_annotation_scores'

workflow ANNOTATION_SCORING {
    take:
    ch_annotation      // (meta, gff, annotation_summary) — ANNOTATE.out.annotation
    ch_target_fasta    // (meta, target_fasta)
    ch_protein_panel   // value channel, protein panel root dir

    main:
    ch_translate_in = ch_annotation.join(ch_target_fasta, by: 0)
    TRANSLATE_ANNOTATION_CDS(ch_translate_in)

    SCORE_ANNOTATION_BLASTP(
        TRANSLATE_ANNOTATION_CDS.out.queries, ch_protein_panel)

    ch_join_in = ch_annotation
        .join(SCORE_ANNOTATION_BLASTP.out.scored, by: 0)
    JOIN_ANNOTATION_SCORES(ch_join_in)

    emit:
    annotation = JOIN_ANNOTATION_SCORES.out.annotation
}
