# Task 33 — Stage 13a: `plant_mt` non-CDS annotation

**Phase:** P4 (from [spec §6](../spec/06-phases.md)), but **blocked** —
see §1. Do not start before its preconditions are met.

## 0. Overview

[Task 30](30_unified_locus_pass.md) gives every target a protein-coding
annotation from the unified miniprot pass, and
[task 31](31_annotate_merge.md) publishes it. A `plant_mt` sample
therefore already receives a real annotation of its ~40 mitochondrial
protein-coding genes, with `status: "ok_cds_only"`. What it lacks is
tRNA and rRNA features.

This is the **lowest-value** of the three targets' remaining gaps:

- [spec §3.1](../spec/03-organelles.md#31-what-differs-between-assembly-targets)
  records that plant mitochondrial loci are "rarely used for barcoding" —
  `plant_mt` exists in the pipeline for completeness.
- Stage 13a is supplementary to the contractual barcode path
  ([spec §2 stage 13](../spec/02-stages.md#2-stage-detail)).
- A plant mitogenome is 200 kb to several Mb of which ~10 % is coding,
  so a feature map is a sparse object regardless.

And it is the **hardest** target to get right, for a reason specific to
this organelle: plant mitogenomes carry plastid-derived insertions
(NUPTs) that will match plastid gene models. A careless tRNA/rRNA search
will annotate plastid-derived tRNAs sitting inside a NUPT as
mitochondrial — the gene-model form of the carry-over problem the
pipeline already fights at
[recruitment](../spec/02-stages.md#215-sibling-organelle-carry-over-in-the-estimate)
and [binning](../spec/03-organelles.md#371-homology-is-measured-against-sibling-panels-too).

This task exists to hold the decision open explicitly rather than let it
be forgotten, and to record what must be true before it is worth making.

**Goal:** add tRNA/rRNA features to the `plant_mt` annotation, once §1's
preconditions are met — or record a decision not to.

---

## 1. Blocked on two things

**Precondition A — a `plant_mt` integration fixture that reaches
`BIN_TARGET`.** `INT-PLANT-01-mt` soft-fails the coverage gate at 28.48×
against the 30× MIN and never gets past stage 6
([task 25](completed/25_coverage_gate_carryover.md),
[task 28](completed/28_plastid_masked_mt_panel.md)). It produces no
`target.fasta`, so a `plant_mt` non-CDS annotator could be built but not
validated on anything except synthetic input. Given the NUPT risk above,
shipping an un-integration-tested annotator for this target is not a
trade worth taking.

The replacement fixture is already tracked in
[`tasks/todo.md`](todo.md) under *Integration fixtures*. That item is
this task's gate; do not duplicate it, and do not add a fourth fixture
or relax the coverage gate for this task's benefit — the soft-fail is
the *correct* verdict on that sample.

**Precondition B — [task 32](32_research_plastid_noncds.md) has
resolved.** If task 32 accepts its §2 option (tRNAscan-SE + barrnap in a
`mulled-build` bundle), this task is **nearly free**: the same two tools
serve both plant targets, and the work reduces to flipping
`params.annotate.plant_mt.non_cds_tool` and adding the NUPT check in §2.
If task 32 selects something plastid-specific instead, this task needs
its own tool decision. Either way, do not decide here first.

## 2. If it proceeds

- **Tools:** whatever task 32 selected, most likely `tRNAscan-SE -O`
  (organellar mode) plus `barrnap`. Plant mitochondrial rRNAs are
  bacterial-like, as plastid ones are, so the bacterial models are
  probably right — verify, and record what was used.
- **The NUPT check is this task's real content.** Cross-reference called
  tRNA/rRNA features against the plastid panel: a feature that matches a
  plastid model better than a mitochondrial one, inside a region with
  plastid homology, is a NUPT artifact. Do **not** silently drop it —
  emit it flagged (`putative_nupt: true` in the GFF attributes and a
  count in `annotation_summary.json`), consistent with how
  [spec §3.7.3](../spec/03-organelles.md#373-coverage-is-a-ranking-signal-not-a-gate)
  handles the same uncertainty at binning: flag rather than filter, so
  the operator sees it ([principle 7](../CONSTITUTION.md)).
- **Do not de-duplicate repeated features.** Plant mitogenomes carry
  genuine large repeats and multiple isoforms; a gene appearing twice
  may be real. Same trap flagged for the plastid inverted repeat in
  [task 31 §4](31_annotate_merge.md).
- **Framework is unchanged.** The output contract, C8, the
  `params.annotate` table and the process wiring all come from
  [task 31](31_annotate_merge.md); this fills in one config entry and
  one merge branch. Any new component gets its own `scripts/tests/`
  module at **100 % branch coverage** ([rule 14](../CONSTITUTION.md)).
- **Gene set:** `assets/organelle_gene_sets.json`
  ([task 29 §5](completed/29_comprehensive_protein_panel.md)) gains `plant_mt`
  tRNA/rRNA expectations so completeness reports as it does for the
  other targets.

## 3. Doing nothing is a legitimate outcome

`ok_cds_only` is a complete, correct, honestly-labelled annotation of
every protein-coding gene found. For a target whose loci are rarely used
for barcoding, on a stage that is supplementary by design,
[rule 19](../CONSTITUTION.md) may well favour leaving it there. If that
is the conclusion, record it in
[spec §3.1/§3.2](../spec/03-organelles.md#31-what-differs-between-assembly-targets)
with the reasoning and close
[spec §8 item 3](../spec/07-open-questions.md#8-remaining-open-questions)
— a deliberate, recorded "no" closes an open question as well as a yes
does.

## 4. Outcomes

_(fill on completion, or on a recorded decision not to implement.)_
