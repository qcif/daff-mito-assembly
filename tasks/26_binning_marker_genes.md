# Task 26 — Marker-gene presence as a `BIN_TARGET` criterion

**Phase:** deferred — P-tune-C ([spec §9 grouping](../spec/07-open-questions.md#9-fine-tuning-post-prototype-benchmarking)),
not P3. This task is **conditional**: §2 is a go/no-go gate that may
close it without implementation.
**Goal:** Decide whether C3 should carry a third, gene-content-based
selection criterion, and if so implement it. Carved out of
[task 23 §5.1](completed/23_bin_target_recalibration.md) so that task can land
without blocking on a reference-bundle addition and a design-principle
boundary call.

**Prerequisite:** [task 23](completed/23_bin_target_recalibration.md) landed and
`-profile integration` green. This task cannot be scoped before then —
§2 depends on observing how much of the problem task 23 already solves.

**Spec basis:**
[§3.7.5](../spec/03-organelles.md#375-orf-integrity-as-implemented-is-vacuous)
(the criterion this would replace),
[§3.7.6](../spec/03-organelles.md#376-revised-per-target-criteria)
(the criteria table it would extend),
[§4.3](../spec/04-reference-data.md#43-protein-panel-for-miniprot_extract)
(the existing protein-panel build machinery),
[§9 item 10a](../spec/07-open-questions.md#9-fine-tuning-post-prototype-benchmarking)
(NUMT/NUPT discrimination, which this is the proposed answer to).

## 1. The gap

Task 23 removes `orf_ok` from C3's selection intersection because
[§3.7.5](../spec/03-organelles.md#375-orf-integrity-as-implemented-is-vacuous)
showed it to be vacuous — every contig in every fixture cleared the
100 aa floor (125–773 aa), and in the `plant_mt` fixture the longest ORF
sat on a *plastid* contig. It is removed rather than recalibrated
because the measure is wrong, not badly tuned.

That leaves C3's remaining criteria on two axes:

| Axis | Signals | Source |
|---|---|---|
| Nucleotide homology to whole organelle genomes | merged aligned fraction vs declared panel; same vs each sibling panel | minimap2, task 23 §2.1–2.2 |
| Structure / abundance | coverage dominance; circularity | task 23 §2.3–2.4, ranking + metadata only |
| **Gene content** | **none** | **this task** |

Removing a criterion without replacing it weakens C3's biological
evidence. The specific exposure is **NUMTs/NUPTs** — nuclear insertions
of organellar DNA, which are homologous by descent and therefore
indistinguishable from the real organelle on the first axis. Task 23
§5.2 accepts this and mitigates by *flagging*
(`low_coverage_candidate`), per
[principle 7](../CONSTITUTION.md); the coverage gate that incidentally
defended against it is gone for the reasons in
[§3.7.3](../spec/03-organelles.md#373-coverage-is-a-ranking-signal-not-a-gate).

Marker-gene presence is orthogonal to both existing axes: it asks
whether the contig *carries an intact set of target-appropriate
organellar protein-coding genes*, which a degraded pseudogene insertion
fails and a genuine organelle passes. It would also give `plant_mt` a
positive discriminator instead of the weak-homology one it currently
relies on (0.165–0.353 merged fraction against a non-conspecific
reference).

## 2. Go/no-go gate — do this first

**Do not build anything before this resolves.** The evidence in
[§3.7.6](../spec/03-organelles.md#376-revised-per-target-criteria)
suggests the sibling-panel test alone may be sufficient: the plastid
contigs in `INT-PLANT-01-mt` score 1.000 against `plant_pt.mmi` versus
0.012–0.064 against `plant_mt.mmi`. That is not a marginal call needing
a tiebreaker — it is a chasm. If it holds on real data, the bundle
addition below is not worth its maintenance cost
([rule 19](../CONSTITUTION.md)).

Assess, once task 23 has run against the P3 datasets rather than the
three integration fixtures:

1. **Is `plant_mt` selection stable?** Does the set of selected contigs
   change under small perturbations of `min_aligned_frac` around 0.15?
   A stable selection means the homology signal is doing real work; an
   unstable one means the threshold is arbitrary and needs a second
   opinion.
2. **Does the declared-vs-sibling margin stay wide?** Plot merged
   aligned fraction against the declared panel vs every sibling panel
   across all P3 samples (this is
   [§9 item 10](../spec/07-open-questions.md#9-fine-tuning-post-prototype-benchmarking)'s
   experiment — share the harness). If the two populations remain
   separable with no overlap, criterion (c) is already sufficient.
3. **Is there measurable NUMT/NUPT contamination?** Requires a sample
   with known NUMT content
   ([§9 item 10a](../spec/07-open-questions.md#9-fine-tuning-post-prototype-benchmarking)).
   If no P3 dataset can demonstrate the failure mode, this task cannot
   demonstrate its own fix either — say so and close it rather than
   building against a hypothetical.

**Outcome:** either *close as unnecessary* with the evidence recorded in
§9 below, or proceed to §3 with a named failure mode the change is
required to fix.

## 3. Design-principle boundary — needs explicit sign-off

Using a protein panel inside binning couples stage 10 to reference data
that stage 13 owns, which sits against
[principle 6](../CONSTITUTION.md) ("assemble first, then extract").

**Argument that it is acceptable:** C3 would use the panel as *evidence
for binning*, not to emit barcodes. The contractual extraction stays at
stage 13, `MINIPROT_EXTRACT` is unchanged, and the binning gene set is
a distinct bundle artefact from `assets/loci.json` — so
[principle 9](../CONSTITUTION.md) is unaffected and the panel can still
evolve with Taxodactyl without re-running assembly.

**Argument against:** it gives a deliberately-scoped stage a second
reference dependency and makes the binning decision sensitive to panel
curation, which is a different maintenance cadence
([§4.3](../spec/04-reference-data.md#43-protein-panel-for-miniprot_extract)
says annually) than the recruit panel.

Judgement is that this is acceptable, but it is a design-principle
boundary and wants explicit sign-off per the constitution's amendment
procedure **before** implementation, not after. Record the decision and
its rationale in [CONSTITUTION.md](../CONSTITUTION.md) principle 6 as a
clarifying note if approved.

## 4. Why not reuse `BLAST_VALIDATE`

The obvious cheap route — C3 uses the RefSeq organelle nucleotide DBs
already built for stage 11
([§4.2](../spec/04-reference-data.md#42-refseq-organelle-dbs-for-blast_validate))
— is rejected, for two reasons. Record them here so the question is not
reopened without new argument.

1. **Same axis, not a new one.** `blastn` against
   `refseq_mt_viridiplantae` measures whole-genome nucleotide homology,
   which is what minimap2 already measures. It is a *more sensitive*
   instrument on that axis (word size 11, local gapped alignment, no
   minimizer seeding) and would likely raise `plant_mt`'s merged
   fraction — but it cannot separate a NUMT from a mitogenome, because
   a NUMT *is* mt sequence by descent and hits RefSeq mt cleanly. It
   does not fill the §1 gap.
2. **It destroys stage 11's independence.** Task 21 positions
   `BLAST_VALIDATE` as an independent check that catches mis-binning
   ([§4.2](../spec/04-reference-data.md#42-refseq-organelle-dbs-for-blast_validate):
   "validation must be against a clean, curated reference set to catch
   mis-binning"). If C3 bins using those same DBs, stage 11 validates a
   decision made with its own evidence against its own database, agrees
   by construction, and the report gains a confirmation that carries no
   information. That is a direct erosion of
   [principle 7](../CONSTITUTION.md) — the operator sees two green
   signals where there is one.

**What BLAST should do instead** (cheap, and worth filing separately if
this task closes at §2): run stage 11 over C3's *rejected* contigs as
well as `target.fasta`, so `secondaries.tsv` carries a human-readable
`stitle` per rejection — *"contig_1, 69 kb, 168×, top hit: chloroplast,
complete genome"*. That is an auditor seeing at a glance why the
plastid was excluded ([rule 16](../CONSTITUTION.md) provenance), and
because it feeds no decision, independence is preserved. It also
supplies [task 25 §3.1](25_coverage_gate_carryover.md)'s
sibling-organelle-fraction reporting for free.

## 5. Reference data — the actual cost

The expensive part is the panel, not the tool. §6's tool choice is
downstream of this and does not change it.

The `plant_mt` barcode panel is far too sparse to bin on — 2 loci
(cox1, nad1) per
[§3.2](../spec/03-organelles.md#32-per-stage-parameter-selection), or 5
(cox1, cob, nad1, atp1, matR) per the
[§4.3](../spec/04-reference-data.md#43-protein-panel-for-miniprot_extract)
schema. **These two sections disagree; reconcile them as part of this
task regardless of its outcome** — it is a one-line spec fix and the
inconsistency will mislead someone.

A binning gene set needs roughly the full organellar proteome:

| Target | Indicative set | Count |
|---|---|---|
| `plant_mt` | nad1–9, cox1–3, atp1/4/6/8/9, ccmB/C/Fc/Fn, cob, matR, mttB, rps/rpl | ~40 |
| `plant_pt` | psa/psb, pet, atp, rbcL, rpo, ndh, rps/rpl | ~80 |
| `animal_mt` | the 13 canonical PCGs | 13 |

Build implications, per [§4](../spec/04-reference-data.md) and
[rule 10](../CONSTITUTION.md) versioning:

- This is the **same build machinery** as
  [§4.3](../spec/04-reference-data.md#43-protein-panel-for-miniprot_extract)
  pointed at a longer symbol list — RefSeq protein query by gene symbol
  + kingdom taxid restriction. That materially lowers the cost estimate
  carried in task 23 §5.1.
- Emit to a **distinct bundle path** (suggest
  `refs/v<ver>/binning_proteins/<target>.faa`), not into
  `proteins/<origin>/`. Keeping the artefacts separate is what makes
  the [principle 9](../CONSTITUTION.md) argument in §3 hold, and lets
  the two panels version independently.
- Fewer representatives per gene than the extraction panel needs —
  presence/absence tolerates 3–5 per gene where extraction wants 5–10
  ([§9 item 11](../spec/07-open-questions.md#9-fine-tuning-post-prototype-benchmarking)).
  Concatenate per target for a single alignment invocation.
- Bundle version bumps; `manifest.json` records per-gene accessions,
  fetch date, and exact query string, as
  [§4.3](../spec/04-reference-data.md#43-protein-panel-for-miniprot_extract)
  already does for the locus panel.
- `scripts/build_refs.sh` extended; the new tree must reach the Azure
  fixture bucket before integration can use it (see
  [task 12](completed/12_azure_integration_fixtures.md) for the upload
  workflow).

## 6. Tool choice — `tblastn` vs `miniprot`

Both need the §5 panel. Decide with the panel cost already sunk.

| | `tblastn` | `miniprot` |
|---|---|---|
| Container | BLAST+ already SHA-pinned for stage 11 — but still a mulled rebuild to add it to C3 alongside `mappy`/`biopython` ([rule 12](../CONSTITUTION.md)) | New tool, mulled rebuild either way |
| Output needed here | presence/absence + bitscore — adequate | precise CDS coords, splice sites, frameshifts — more than binning needs |
| Divergent homologs | weaker | stronger (63.1 % base sensitivity vs Spaln2's 55.7 %, [spec/miniprot.md](../spec/miniprot.md)) — matters against non-conspecific references |
| Genetic code | `-db_gencode` | `--trans`, and already exercised for the animal table 2/5 trial (§3.3) |
| Maintenance | two protein-alignment idioms in the codebase | one tool, one idiom, shared with stage 13 |

**Lean: miniprot**, despite being more than binning strictly needs. The
container rebuild is required either way, so `tblastn`'s "already
pinned" advantage largely evaporates; sensitivity on distant homologs
is the property this criterion depends on; and a single
protein-alignment idiom is the lower-maintenance outcome
([rule 19](../CONSTITUTION.md)). Confirm during implementation — if
miniprot's runtime on a 15-contig raw assembly against ~80 plastid
proteins proves material, revisit.

**Note the direction of alignment.** Stage 13 aligns proteins to
`target.fasta` (one binned contig set). Here the subject is **every
METAFLYE contig**, pre-selection — a larger, messier input. Measure
runtime on the `plant_mt` fixture before assuming it is trivial.

## 7. Change to `bin/bin_target.py`

Criterion (c) in
[§3.7.6](../spec/03-organelles.md#376-revised-per-target-criteria)
becomes a fourth condition, not a replacement for the sibling test:

```
# pseudocode — a contig is a target candidate where all of:
#   1. merged aligned frac vs declared panel >= min_aligned_frac
#   2. merged identity vs declared panel      >= min_identity
#   3. merged aligned frac vs declared panel  >  every sibling panel
#   4. marker_gene_count                      >= min_marker_genes   <-- new

hits    = align_proteins(contig_fasta, panel_faa, trans=code)
genes   = {h.query_gene for h in hits if h.identity >= MIN and h.cov >= MIN}
n_genes = len(genes)
```

Design constraints:

- **`min_marker_genes` is per-target and lives in `nextflow.config`**,
  alongside the `params.bin_target_thresholds` map task 23 introduces —
  not as a Python constant ([rule 18](../CONSTITUTION.md)). A plant
  mitogenome fragmented into sub-genomic contigs may carry 1–2 genes
  per contig, so the floor is low and target-specific; `animal_mt`,
  where the whole genome is one contig, can demand many more.
- **Do not apply the floor per-contig on `plant_mt` without thought.**
  `plant_mt` emits *all* candidates, and genuine sub-genomic contigs may
  be gene-poor. Consider a set-level rule (the emitted set collectively
  carries ≥ N genes) rather than a per-contig one, and record which was
  used. Getting this wrong reintroduces exactly the false negative that
  task 23 exists to remove.
- **Record, never silently filter.** Per contig in `secondaries.tsv` and
  `bin_metadata.json`: `marker_genes_found` (the symbol list, not just a
  count), `marker_gene_count`, and the panel version used
  ([rule 16](../CONSTITUTION.md)). An auditor must be able to see *which*
  genes drove the call.
- **Guard a missing panel** exactly as task 23 §2.2 guards a missing
  sibling `.mmi`: warn to stderr, record `marker_panel: null`, fall back
  to the three-criterion decision, exit 0. An older bundle must not fail
  the sample ([principle 8](../CONSTITUTION.md)).
- Keep `orf_aa_len` and `selected_genetic_code` as reported diagnostics —
  task 23 already demoted them and this does not restore them.

## 8. Tests

**Unit** — `scripts/tests/test_bin_target.py`, extending task 23's
cases, to [rule 14](../CONSTITUTION.md) 100 % branch coverage of the
selection path. Synthetic inline sequences only, no checked-in fixture
data ([rule 19](../CONSTITUTION.md)):

| # | Scenario | Expected |
|---|---|---|
| 1 | Contig passes criteria 1–3, carries ≥ `min_marker_genes` | `target_candidate` |
| 2 | Contig passes 1–3, carries 0 marker genes | rejected; new classification value (suggest `no_marker_genes`), distinct from `off_target` and `sibling_organelle` |
| 3 | Contig carries marker genes but fails the sibling test | `sibling_organelle` — criterion 3 still dominates |
| 4 | Marker panel absent from bundle | warns, `marker_panel: null`, selection falls back to 1–3, exit 0 |
| 5 | `plant_mt` set-level rule: 3 gene-poor contigs collectively above the floor | all emitted (guards the §7 sub-genomic false negative) |
| 6 | Per-target `min_marker_genes` from CLI args | same contig admitted under `plant_mt`, rejected under `animal_mt` |
| 7 | Marker hit below identity/coverage floor | not counted toward `marker_gene_count` |

Case 5 is the regression test for the failure mode most likely to be
introduced by this change.

Run via `scripts/pytest.sh`; flake8 per the project venv at 79 cols.

**Integration** — `tests/integration/assertions.sh`, BIN_TARGET block:

- every sample still selects ≥ 1 contig (task 23's assertion must not
  regress — this is the primary risk of adding a criterion);
- `INT-PLANT-01-mt` still emits contig_3 + contig_5 and not contig_1 or
  contig_6;
- `bin_metadata.json` records a non-empty `marker_genes_found` for every
  emitted contig on `INT-ANIMAL-01` and `INT-PLANT-01-pt`.

`expected/*/bin_bounds.json` should need no change. **If adding this
criterion moves any fixture outside its bounds, that is a finding, not
a bound to retune** — report it before touching the file.

**Fixture prerequisite:** the `binning_proteins/` tree must be present
in the integration reference bundle before these assertions can run. If
it is not, the integration block is gated the same way task 21 gates on
the `validate/` tree.

## 9. Exit criteria

Either:

- **Closed at §2** — evidence recorded that the sibling-panel test alone
  is sufficient, with the stability and margin measurements from §2.1–2.3
  written into an Outcomes section here and
  [§9 item 10a](../spec/07-open-questions.md#9-fine-tuning-post-prototype-benchmarking)
  updated to reflect that NUMT discrimination remains accepted-and-flagged;
  the §4 BLAST-over-rejected-contigs improvement filed as its own task;
  the §5 `plant_mt` locus-count discrepancy fixed.

Or:

- §3 sign-off recorded against [principle 6](../CONSTITUTION.md).
- `binning_proteins/` tree built, versioned, manifested, and uploaded;
  `scripts/build_refs.sh` extended.
- C3 container rebuilt and SHA-pinned ([rule 12](../CONSTITUTION.md)).
- `min_marker_genes` in `nextflow.config`, not in `bin/bin_target.py`.
- `bin_metadata.json` records `marker_genes_found`, `marker_gene_count`,
  and panel version per contig.
- All three fixtures still select ≥ 1 contig; no fixture bound retuned
  to accommodate the change.
- 100 % branch coverage on the selection path; `scripts/pytest.sh`
  green; flake8 clean at 79 cols.
- `-stub-run` still green — stub block unchanged.
- [§3.7.5](../spec/03-organelles.md#375-orf-integrity-as-implemented-is-vacuous)
  and [§3.7.6](../spec/03-organelles.md#376-revised-per-target-criteria)
  updated with the criterion as implemented and the observed values that
  set `min_marker_genes`;
  [§9 item 10a](../spec/07-open-questions.md#9-fine-tuning-post-prototype-benchmarking)
  closed.

## 10. Not in scope

- **Any change to `MINIPROT_EXTRACT` or `assets/loci.json`.** The
  extraction contract at stage 13 is untouched; the binning panel is a
  separate artefact ([principle 9](../CONSTITUTION.md)).
- **Annotation-based binning** (MITOS2/GeSeq gene content). Stage 12
  runs downstream of binning; inverting that ordering is a pipeline
  restructure, not a criterion change.
- **Retuning task 23's homology thresholds.** That is
  [§9 item 10](../spec/07-open-questions.md#9-fine-tuning-post-prototype-benchmarking),
  benchmarked on the same P3 datasets but independently.
- **Emitting marker-gene hits as barcodes.** Presence is evidence for a
  binning decision; nothing from this stage reaches the report as a
  recovered locus.

## 11. Outcomes

- §2 gate result (proceed / close), with the measurements: —
- §3 principle-6 sign-off (who, when, rationale): —
- Tool selected and why: —
- Observed `min_marker_genes` per target: —
- Panel build: gene counts, RefSeq release, bundle version: —
