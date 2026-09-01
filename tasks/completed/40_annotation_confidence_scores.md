# Task 40 — Annotation confidence: one comparable score per CDS feature

**Phase:** P4a (an addition to stage 13a).
**Prerequisite:** task 39_annotate_rescue_annotator_only.md — the merged
GFF only carries mixed provenance once rescue lands, which is what makes
a comparable score necessary rather than merely nice.

**Goal:** Attach a confidence figure to every protein-coding feature in
the merged annotation, on a scale that is comparable between miniprot
and the specialist annotator, so a reviewer can tell a confident call
from a marginal one without knowing which tool produced it.

---

## 1. Overview

After task 39 the annotation GFF carries CDS features from two
independent methods. A reviewer reading it — or the per-sample
`report.html` built on it — has no way to judge any individual call.
Nothing in the current output distinguishes a near-certain `COX1` from a
speculative rescue.

The obvious material is already there but is not usable as-is, because
neither tool emits a portable score:

- **miniprot** column 6 is its own dynamic-programming alignment score.
  It is not comparable to anything outside miniprot, and it is not even
  stable across settings — correcting the translation table in task 38
  moved ATP6 from 776 to 1012 on the same locus with identical
  coordinates.
- **MITOS2**'s exon score column is unusable as a confidence value. On
  `CLIENT-BC05` it reports `atp8` at `1172369.5` and `atp6` at
  `194762004.0`; on `INT-ANIMAL-01`, `atp8_1` at `2511.5` against
  `cox1_0-b` at `749869353.8`. The scale is undefined and length-driven.

What a reviewer actually wants is a single number combining **how
similar** the call is to a known protein and **how much of that protein
it covers**. A gene can look poor on either axis for opposite reasons: a
truncated fragment of a well-referenced gene (high identity, low
coverage) is a very different problem from a complete gene in an
under-referenced clade (low identity, high coverage). The second is
common in this workflow and must not read as failure.

That statistic exists and is standard: the **bitscore**. It is a
log-odds sum over aligned positions, so it grows with both identity and
alignment length, and — unlike an E-value — it is independent of
database size. E-values from the 13-gene protein panel would not be
comparable to E-values from `refseq89m`, so **bitscore is the field to
carry; E-value is not.**

The `CLIENT-BC05` `ATP8` case is the worked example. Under a bitscore +
identity + coverage triplet it reads as ~100 % query coverage at ~47 %
identity — correctly, "confidently the right gene, distant from our
references" — rather than as the absence it currently reports.

**Exit criteria:**

- Every CDS feature in the merged annotation carries identity, query
  coverage and bitscore against the same reference panel, regardless of
  which tool called it.
- The figures appear in `annotation_summary.json` and are available to
  the report.

---

## 2. Mechanism: one `blastp` pass over the emitted CDS

BLAST 2.17.0 is already SHA-pinned in
[conf/containers.config](../conf/containers.config) for
`BLAST_VALIDATE`, so no new dependency is introduced
([rule 12](../CONSTITUTION.md)).

```
# pseudocode
translate every CDS feature in the merged GFF, under the sample's
    selected genetic-code table (task 38)
blastp -query <translated CDS> -subject <the same per-gene panel .faa> \
       -outfmt "6 qseqid sseqid pident qcovhsp evalue bitscore" \
       -max_target_seqs 1
join the best hit per feature back onto the GFF / summary
```

Cost is negligible: on `animal_mt`, thirteen 50–570 aa queries against
the 130-sequence panel.

Two implementation notes:

- **Translate from coordinates, under the selected table.** miniprot's
  `--trans` output is available for its own features, but rescued
  features have no equivalent, and using two different translation paths
  would defeat the point of a comparable score. Translate uniformly from
  the assembly using the table task 38 selects. `bin/validate_barcodes.py`
  already does coordinate-driven, CIGAR-aware translation for C5 —
  reuse rather than reimplement ([rule 19](../CONSTITUTION.md)).
- **Blast against the panel, not `refseq89m`.** The panel is what
  miniprot aligned to; scoring against a different database would
  measure the database, not the call.

## 3. Which number to surface

Report **`pident`, `qcovhsp` and `bitscore` as three fields.** For a
13-row table that is more informative than any single composite, and it
lets a reviewer see the two axes independently — which is exactly the
distinction §1 says matters.

If a single sortable number is wanted (for a report column, or a future
threshold), use **normalised bitscore**: the feature's bitscore divided
by the query's self-bitscore. That lands in 0–1 and is comparable across
genes of very different lengths, where a raw bitscore is not. Add it
only if there is a consumer; do not introduce a derived field with no
reader ([rule 19](../CONSTITUTION.md)).

**Do not carry E-value** as a confidence field, for the database-size
reason in §1. If it is emitted at all, document what database it refers
to.

## 4. Where the process runs — a container decision

[spec §2.2](../spec/02-stages.md) records C8 as sharing the `ANNOTATE`
MITOS2 biocontainer and being stdlib-only. Adding `blastp` to stage 13a
breaks that, and there are three options. Pick one and record it:

1. **A separate process** between `ANNOTATE` and its consumers, using
   the existing pinned BLAST container. Simplest, no new image, one more
   channel hop.
2. **A mulled MITOS2 + BLAST image.** Keeps stage 13a as one process;
   [rule 14](../CONSTITUTION.md) explicitly anticipates mulled images
   where one process needs multiple biocontainer packages. Costs an
   image to build and maintain.
3. **Fold into `BLAST_VALIDATE` (stage 11).** Rejected on ordering —
   stage 11 runs before the annotation exists.

Option 1 is the recommended default under
[rule 19](../CONSTITUTION.md): no new artefact to maintain, and it keeps
the scoring independently testable and independently skippable.

Whichever is chosen, C8 stays stdlib-only and
[rule 11](../CONSTITUTION.md) applies — pinned tag, no `latest`.

## 5. Scope guard: scores are data, not a gate

No feature is dropped, flagged or re-called on the basis of its score,
and no sample status depends on it. This is the same constraint task 31
placed on annotation counts: annotation is supplementary, it never gates
a sample, and the negative-clarity classification stays centralised in
`COLLATE`. Introducing a confidence threshold here would create a second
decision surface — explicitly out of scope, and if one is ever wanted it
needs the benchmark set ([spec §9](../spec/07-open-questions.md)), not
this task.

The one permitted consequence is presentational: the report may sort or
style by score.

## 6. Unit tests

The join/parse logic is custom logic under
[rule 14](../CONSTITUTION.md) and needs a `scripts/tests/` module at
**100 % branch coverage**:

- Best-hit selection when a query has multiple HSPs and multiple
  subjects.
- A query with **no** blastp hit at all — must yield an explicit null,
  not a zero, and must not drop the feature from the GFF.
- Features from both sources scored identically by the same code path.
- Malformed/truncated blastp output.
- Normalised bitscore arithmetic, if §3's optional field is implemented,
  including the self-hit denominator.
- The no-annotator targets (`plant_pt`, `plant_mt`) scoring their
  miniprot-only features normally.

Hand-written tabular fixtures, not captured BLAST output
([rule 19](../CONSTITUTION.md)).

## 7. Integration reconciliation

- Assert that every CDS feature in `INT-ANIMAL-01`'s annotation carries
  a score triplet, and that no feature was dropped by the scoring step
  (count in == count out).
- Keep bounds loose. Identity and coverage against a fixed panel are
  stable, but pinning exact bitscores would make the suite fragile to
  any BLAST version bump — assert presence, type and plausible range,
  not values.

## 8. Spec updates

- [spec §2 stage 13a](../spec/02-stages.md) — the scoring step and where
  it runs (§4).
- [spec §2.2](../spec/02-stages.md) — the new component row, or the
  amended C8 row.
- [spec §6a](../spec/06a-reports.md) — the report surfaces the score
  alongside the source, so a reviewer sees provenance and confidence
  together. This complements the standing `RUN_REPORT` item in
  `tasks/todo.md` about surfacing `cds_crosscheck` disagreements.

## 9. Acceptance criteria

- Every CDS feature in the merged GFF carries `pident`, `qcovhsp` and
  `bitscore` from a single scoring path, whatever tool called it.
- Unscoreable features carry an explicit null and survive into the GFF.
- No feature is filtered and no sample status changes as a result of
  this task.
- Container choice made, recorded in the spec, pinned, and CI's
  container lint green.
- `scripts/pytest.sh` green at 100 % branch coverage; `flake8` clean;
  `-profile stub -stub-run` green; `-profile integration` green.

## 10. Out of scope

- **Any threshold on the score** — §5.
- **Scoring non-CDS features.** tRNA/rRNA confidence is a different
  problem with different tooling (MITOS2's `mitfi` rows do carry real
  E-values); revisit separately if a reviewer asks for it.
- **Replacing the barcode `--min-identity` floor.** The barcode side has
  its own identity gate and its own open question in `tasks/todo.md`
  about whether coverage should replace identity there. This task may
  inform that decision; it does not make it.
- **Broadening the protein panel** — the standing panel-breadth item in
  `tasks/todo.md`.

## 11. Outcomes

Implemented largely as briefed, with one deviation from §4 discovered
during scoping and resolved with the user before writing any code.

### 11.1 Deviation: §4's option 1 needed splitting into three processes

**As briefed:** §4 recommended "a separate process ... using the
existing pinned BLAST container. Simplest, no new image, one more
channel hop" as the default under rule 19.

Checked directly before implementing:

```
$ docker run --rm quay.io/biocontainers/blast@sha256:78e8ebf2... \
    sh -c 'which python3 || echo "no python3"'
no python3
```

The pinned `blast:2.17.0` biocontainer has **no Python interpreter at
all** — not even a stripped one. Task 40's mechanism needs Python for
two things a single BLAST-only process cannot do: coordinate-driven
translation (needs `biopython`) and the GFF/`annotation_summary.json`
join (needs `json`). A literal single-process option 1 was not
possible without either adding Python to that image (which is not
"the existing pinned container" any more) or moving translation/join
logic into bash/awk running inside it (a much larger and more fragile
piece of custom logic than rule 14 anticipates for a mechanical step).

Presented to the user as a choice between (a) three processes reusing
two already-pinned images with zero new build/maintenance cost, or (b)
extending the wf5-scripts Dockerfile with the `blast` package to keep
it to one process. **User chose (a).** Implemented as:

1. `TRANSLATE_ANNOTATION_CDS` (`bin/translate_annotation_cds.py`,
   shares `EXTRACT_BARCODES`'s `neoformit/daff-wf5-scripts` image) —
   parses the merged GFF, splices + translates every CDS feature
   uniformly (reusing `validate_barcodes.extract_cds_seq`, not
   reimplementing it — rule 19), writes one query FASTA per gene.
2. `SCORE_ANNOTATION_BLASTP` (plain shell, reuses `BLAST_VALIDATE`'s
   pinned BLAST 2.17.0 image) — one `blastp -subject` call per gene
   against that gene's own panel file, so a hit can never cross into a
   different gene family.
3. `JOIN_ANNOTATION_SCORES` (`bin/annotation_scores.py`, same
   scripts image as step 1) — best-hit selection (bitscore, ties on
   identity), writes the scored GFF + `annotation_summary.json`'s new
   `cds_scores` list. This step's output is the annotation stage's
   final, published artifact — `ORGANELLE_MAP` and `COLLATE` were
   rewired to consume it instead of `ANNOTATE.out.annotation` directly.

A fourth module, `bin/annotation_gff.py`, factors out the merged-GFF
parser shared by steps 1 and 3 (sibling import, same mechanism
`intervals.py` already establishes) so both agree on exactly one
definition of "one CDS feature" regardless of provenance.

Net: **no new container image** — the outcome rule 19 was steering
toward — at the cost of two extra channel hops beyond what §4
anticipated.

### 11.2 Other notes

- Per §3, only the `pident`/`qcovhsp`/`bitscore` triplet was added —
  normalised bitscore was not implemented, since no consumer asked for
  it yet (rule 19).
- `annotation_summary.json`'s `$schema` bumped `v2` → `v3` in the join
  step, for the same reason task 39/31's additions were versioned.
- Measured on `INT-ANIMAL-01` (`-profile integration`): 13/13 CDS
  features carry a score triplet (`cds_scores` count == GFF mRNA
  count). The rescued `ATP8` scores 42.9 % identity / 64 % query
  coverage / bitscore 31.6 against its own panel — the same shape as
  this task's worked example (confidently the right gene, distant
  from the reference), not a `not_found`.
- `bin/translate_annotation_cds.py` and `bin/annotation_scores.py`
  needed `chmod +x` — missed on first `-profile integration` run
  (exit 126, "Permission denied"), caught and fixed before the run
  that's recorded above.
- New `scripts/tests/` modules (`test_annotation_gff.py`,
  `test_translate_annotation_cds.py`, `test_annotation_scores.py`) hit
  100 % branch coverage; full suite (300 tests) green.
- `tests/integration/assertions.sh` §7 additions: `cds_scores` count
  matches CDS feature count, every score is within a plausible range,
  every mRNA row in the GFF carries a `pident=` attribute. All pass.

### 11.3 Follow-up: wrapped in a subworkflow (2026-09-02)

Post-completion cleanup, prompted by a question about `main.nf`
readability rather than by any defect. The three-process chain was the
one block in `main.nf` with internal plumbing that carried no
stage-level meaning — a `.join()` to re-attach `target_fasta` (needed
because `ANNOTATE` doesn't emit it) and a second `.join()` to
re-attach the blastp TSV before the final join step. Neither join is
something a reader of the top-level flow needs to see.

Added `subworkflows/local/annotation_scoring.nf` — `ANNOTATION_SCORING`
takes `(ch_annotation, ch_target_fasta, ch_protein_panel)` and emits
`annotation` in the same `(meta, gff, summary)` shape `ANNOTATE.out
.annotation` already has, so it's a drop-in decorator on that channel.
`main.nf`'s stage 13a block goes from 3 includes + ~18 wiring lines to
1 include + a single call; `ORGANELLE_MAP` and the `ch_ok_inputs` join
both now name `ANNOTATION_SCORING.out.annotation` instead of the last
process in the chain, so a future 4th scoring step (§3's optional
normalised bitscore, or §10's deferred non-CDS scoring) changes only
the subworkflow file.

The three processes themselves are unchanged and stayed in
`modules/local/` (not moved into the subworkflow directory) —
CI's container-coverage lint globs `modules/local/*.nf` and expects a
`containers.config` entry per stem; moving them would have broken that
check for no benefit, since they still need per-process container
pins regardless of which `.nf` file includes them.

**Verified:** `-stub-run` green end to end. `-profile integration
-resume` re-ran only the three scoring processes and everything
downstream of them (process identity is namespaced by the
subworkflow, so Nextflow's cache key changes) — everything upstream of
`ANNOTATE` stayed cached, confirming the refactor didn't touch
anything outside stage 13a cont'd. `tests/integration/assertions.sh`
still green, including all of §7's new assertions. Container-coverage
and bare-params lint checks re-run locally, still green.
