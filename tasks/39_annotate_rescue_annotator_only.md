# Task 39 — `ANNOTATE`: rescue annotator-only CDS calls into the merged GFF

**Phase:** P4a (a revision to the stage 13a merge built in task 31).
**Prerequisite:** task 38_miniprot_genetic_code.md — the cross-check
buckets this task consumes are computed from miniprot alignments, and
task 38 changes which genes appear in them.

**Goal:** Where the specialist annotator calls a protein-coding gene that
miniprot missed entirely, emit it as a real feature with its source
recorded, instead of noting it in a JSON field and dropping it.

---

## 1. Overview

Stage 13a merges two annotations. CDS features come from miniprot
(`cds.gff`); tRNA and rRNA features come from a specialist annotator —
MITOS2 on `animal_mt`. MITOS2's own CDS calls are deliberately **not**
used as features; they are retained only as an independent cross-check,
bucketed in `annotation_summary.json` as `agreed`, `miniprot_only`,
`annotator_only` and `coordinate_conflicts`.

That design is sound where the two methods overlap. It has one blind
spot, and a client query found it.

On a *Bactrocera dorsalis* intercept (`CLIENT-BC05`) the annotation
reported `ATP8` missing and protein-coding completeness 12/13. The gene
is not missing. It is present and intact at `contig_1:11225-11377` — a
53-codon ORF with no internal stops under table 5, the canonical
`MPQMAP…` / `…NWKW` termini, and the classic 7 bp `ATGATAA` ATP8/ATP6
overlap of insect mitogenomes. Read depth across the locus is a flat
~190×, matching its flanks. MITOS2 called it correctly.

miniprot did not, and cannot. `ATP8` is the fastest-evolving metazoan
mitochondrial protein: at **53 aa** with the nearest panel reference
only **47 % identical**, there is not enough signal to seed an
alignment. Re-running miniprot against the `ATP8` panel alone, under
table 5 and with deliberately permissive seeding
(`-S -n 1 -m 5 -p 0.1 --outs 0.1 -k 4 -c 200000`), still returns zero
alignments. Supplying a near-identical query does hit immediately, so
this is a divergence-and-length limit, not a bug in the invocation and
not something a parameter sweep will recover.

So the merge as specified will systematically drop short, fast-evolving
protein-coding genes on any under-referenced taxon — and the
under-referenced intercept is the case this workflow exists to serve.
The information was never lost: it sat in `cds_crosscheck.annotator_only`
while the same JSON reported the gene missing. Under
[principle 18](../CONSTITUTION.md) a result that ships a present,
intact gene as absent is the failure mode to close.

**The change is narrow.** Only the `annotator_only` bucket is rescued —
calls with no miniprot counterpart at all. `coordinate_conflicts` stays
untouched and unactioned: where both methods called a gene and
disagreed on where it is, miniprot still wins and the disagreement is
still just reported. That preserves the single-CDS-source coherence
task 30 bought for the barcodes, because a rescued call is by definition
a locus miniprot produced nothing for, and therefore cannot conflict
with a barcode.

**Exit criteria:**

- `INT-ANIMAL-01` and `CLIENT-BC05`-shaped inputs emit annotator-rescued
  CDS features in the merged GFF, attributed to their source.
- Completeness counts rescued genes, so `ATP8` no longer reports missing
  when it was found.
- Every rescued feature is distinguishable from a miniprot feature in
  both the GFF and the summary.

---

## 2. Four guards — the obstacles are all in MITOS2's output shape

None of these are reasons not to do the rescue. Each is a way to do it
wrong, and each is visible in current fixture output.

### 2.1 Fragmented gene names

When MITOS2 splits a gene across fragments it suffixes the copies. From
`INT-ANIMAL-01`'s MITOS2 output:

```
gene  1719 1919  ID=gene_atp8_1;Name=atp8_1
gene  1897 2031  ID=gene_atp8_0;Name=atp8_0
gene  2881 4197  ID=gene_cox1_0;Name=cox1_0
gene  4190 4381  ID=gene_cox1_1;Name=cox1_1
```

`normalise_gene_symbol()` in `bin/annotate_summary.py` strips
underscores, so `atp8_0` folds to `ATP80` — matching no canonical name
in `assets/organelle_gene_sets.json`. Rescued blind, that emits **two**
ATP8 features *and* still reports ATP8 missing from completeness, which
is worse than the current behaviour.

Needs a fragment-suffix rule (`_\d+$`) applied in the normalisation
path. Apply it carefully: `nad4l` is a real gene name whose tail must
survive, and the existing anticodon-stripping behaviour for tRNAs must
not regress.

### 2.2 Multi-exon calls carry phase on the exon rows

MITOS2 puts the gene span on the `gene` row and the reading frame on
`exon` rows:

```
exon  1197 1481  ...  0  Parent=transcript_atp6;Name=atp6-b
gene  1197 1615  ...  .  ID=gene_atp6;Name=atp6
exon  1496 1615  ...  0  Parent=transcript_atp6;Name=atp6-a
```

`parse_mitos_gff()` currently captures only `gene` rows, which carry
`.` in the phase column. A CDS feature emitted from a `gene` row alone
is not a valid CDS feature. Rescue must assemble the `exon` rows grouped
by `Parent`, preserving phase and ordering, and emit an mRNA/CDS block
matching the shape `write_gff()` already produces for miniprot.

### 2.3 Off-panel genes

MITOS2 calls features outside the canonical set — `INT-ANIMAL-01` has
`gene_lagli`, a LAGLIDADG homing endonuclease. A blanket rescue injects
those into the organelle GFF.

Decide explicitly and state it in the spec. The conservative default,
and the one consistent with [rule 19](../CONSTITUTION.md): rescue only
calls that canonicalise into the target's `protein_coding` set, and
leave everything else in `annotator_only` as today. A rescued feature
should be one the canonical gene set vouches for.

### 2.4 `annotator_only` is a name test, not a position test

The bucket is currently decided purely by absence of a matching **name**
in the miniprot winners. A call named `cox1_0` sitting directly on top
of miniprot's `COX1` is "annotator only" by name while being a duplicate
by position — and §2.1's suffix fix does not fully close this, because
any name that fails to canonicalise lands in the same trap.

Add a coordinate guard: reject any rescue candidate that overlaps an
existing miniprot CDS on the same strand. `overlaps()` already exists in
`bin/annotate_summary.py` and can be reused.

## 3. Provenance

Use the GFF `source` column — it is the format-native place for this and
costs nothing. miniprot features keep `miniprot`; rescued features carry
the annotator (`mitos`). A reader can then separate the two without
consulting the summary JSON.

In `annotation_summary.json`:

- Add a `cds_rescued` list (gene names promoted from `annotator_only` to
  features), so the rescue is auditable and testable.
- `annotator_only` keeps only the calls **not** rescued, with the reason
  they were held back (§2.3 off-panel, §2.4 overlap).
- `cds_source` becomes insufficient as a single string once two sources
  contribute. Either widen it to a list or add a per-source feature
  count; do not leave it reading `"miniprot"` for a mixed GFF.
- `protein_coding_completeness` and `protein_coding_genes_missing` must
  count rescued genes. This is the field the client actually read.

## 4. The integration invariant this breaks

`tests/integration/assertions.sh` asserts, from task 31 §1:

> every CDS row in the merged annotation GFF traces exactly to a
> `cds.gff` row — no CDS feature ever originates from MITOS2

That invariant is the thing this task deliberately relaxes. It must be
**narrowed, not deleted** — its regression value is real. The
replacement: every CDS row with `source == miniprot` traces exactly to a
`cds.gff` row, and every CDS row that does not must appear in
`cds_rescued`. That still catches the accident the original was written
to catch (MITOS2 silently supplanting miniprot on a gene both called)
while permitting the rescue.

## 5. Spec updates

[spec §2 stage 13a](../spec/02-stages.md) currently states MITOS2's CDS
calls are "discarded as features and retained as an independent
cross-check", and the C8 row in
[spec §2.2](../spec/02-stages.md) says the cross-check records
single-source calls "without acting on them". Both are now false and
must be amended: the cross-check still records everything, and acts on
exactly one bucket, under the §2 guards.

Keep the stated rationale for *why* CDS is single-sourced — it is still
correct for every gene both methods find, and it is what guarantees the
barcode is always a feature on the map.

## 6. Unit tests

`scripts/tests/test_annotate_summary.py` extends to cover, at
**100 % branch coverage** ([rule 14](../CONSTITUTION.md)):

- Rescue of a clean single-exon `annotator_only` call (the `ATP8` case).
- Multi-exon rescue with phase preserved (§2.2).
- Fragment-suffixed names canonicalising correctly, and a fragmented
  gene not producing duplicate features (§2.1).
- An off-panel call (`lagli`) held back, with a reason (§2.3).
- An overlapping call held back (§2.4), including the same-name-different
  -strand case.
- Completeness before/after rescue, asserting the gene leaves
  `protein_coding_genes_missing`.
- `cds_source` / per-source counts on a mixed GFF.
- The no-annotator path (`plant_pt`, `plant_mt`) unchanged —
  `status: "ok_cds_only"`, nothing rescued, no new fields populated.

Fixtures stay hand-written minimal GFF fragments
([rule 19](../CONSTITUTION.md)).

## 7. Integration reconciliation

- `tests/integration/expected/animal_mt/annotation_bounds.json` —
  `min_cds` may rise; re-baseline after task 38 has landed, not before,
  or the two changes will be measured on top of each other.
- Add an assertion that `cds_rescued` and the GFF `source` column agree.
- Consider promoting the `CLIENT-BC05` `ATP8` case to a fixture
  assertion. It is the clearest available example of the failure and
  currently exists only as a manual run under `tests/manual/output/`.
  Weigh against the fixture-drift cost ([rule 19](../CONSTITUTION.md))
  and the fixture-provenance item in `tasks/todo.md`.

## 8. Acceptance criteria

- An `annotator_only` protein-coding call with no overlapping miniprot
  CDS is emitted as a CDS feature attributed to its annotator.
- Rescued genes count toward `protein_coding_completeness`.
- All four §2 guards are implemented and individually tested.
- The §4 invariant is narrowed and still green.
- Spec §2 stage 13a and §2.2 C8 reconciled.
- `scripts/pytest.sh` green at 100 % branch coverage for C8; `flake8`
  clean; `-profile stub -stub-run` green; `-profile integration` green.

## 9. Out of scope

- **`coordinate_conflicts`.** Still recorded, still unactioned. Where
  both methods found a gene, miniprot remains authoritative.
- **`miniprot_only`.** Unchanged — miniprot features are already
  emitted.
- **Any decision logic keyed on annotation output.** A rescued or
  missing gene fails nothing and re-bins nothing; counts remain data for
  the report and classification stays centralised in `COLLATE`
  (task 31 §0's standing constraint).
- **Confidence scores for the rescued features** —
  task 40_annotation_confidence_scores.md.
- **Non-CDS rescue.** tRNA/rRNA already come from the annotator; there
  is no equivalent gap.
- **Plant branches.** Neither has a specialist annotator yet
  (32_research_plastid_noncds.md, 37_annotate_plant_mt_noncds.md), so
  there is nothing to rescue. The code path must simply stay inert.
