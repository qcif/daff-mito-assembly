# Task 28 — Plastid-masked `plant_mt` recruitment panel

> **Reconcile with [task 26](26_binning_marker_genes.md) before starting
> either.** Both tasks answer the same underlying problem — *whole-genome
> nucleotide homology cannot cleanly separate the two plant organelles* —
> by different means. Task 26 adds a **gene-content** criterion to C3;
> this task fixes the **reference panel** the existing homology criterion
> measures against. They are not independent:
>
> - Task 26's §2 go/no-go gate asks whether the declared-vs-sibling
>   margin stays wide on real data. This task moves that margin
>   (measured: `plant_mt` contigs from 0.28 → 0.98 merged aligned
>   fraction, with plastid margins held at −0.89 to −0.93). Running
>   task 26's gate against the *current* single-genome panel would
>   answer the wrong question.
> - If this task lands first and the margin is wide and stable, task 26
>   may close unimplemented — which is its own §2 stated preference,
>   and the cheaper outcome under [rule 19](../CONSTITUTION.md).
> - If task 26 lands first, its marker-gene criterion would be tuned
>   against a panel this task then replaces, and would need
>   recalibrating.
>
> **Recommended sequence: 28 before 26.** Task 26's §2 should consume
> this task's post-change measurements. Whoever picks up either task
> first must read both and record the sequencing decision.

**Phase:** P-tune-C ([spec §9 grouping](../spec/07-open-questions.md#9-fine-tuning-post-prototype-benchmarking)).
Reference-data change; no workflow logic moves.
**Goal:** Replace the single-sequence `plant_mt` recruitment panel with a
broad, plastid-masked panel built from RefSeq, without re-opening the
sibling-organelle defects that
[task 23](completed/23_bin_target_recalibration.md) and
[task 25](completed/25_coverage_gate_carryover.md) closed.

**Prerequisite:** [task 25](completed/25_coverage_gate_carryover.md)
landed (C2 gates on target-assigned bases; C3 records
`sibling_carryover`). This task's success metric is expressed in fields
those two tasks introduced.

**Spec basis:**
[§4.1](../spec/04-reference-data.md#41-setup-task-source-reference-panel-from-getorganelle)
(the panel source, and its own "re-evaluate after P1" clause),
[§4.2](../spec/04-reference-data.md#42-refseq-organelle-dbs-for-blast_validate)
(the RefSeq Viridiplantae mitochondrion set this reuses),
[§4.4](../spec/04-reference-data.md#44-consolidated-build-script)
(bundle layout + `manifest.json`),
[§2.1.5](../spec/02-stages.md#215-sibling-organelle-carry-over-in-the-estimate)
(the coverage split that consumes the panel),
[§3.7.1](../spec/03-organelles.md#371-homology-is-measured-against-sibling-panels-too)
(the sibling discrimination that consumes it),
[rule 10](../CONSTITUTION.md) (versioned, immutable, provenance-tracked
bundles), [rule 19](../CONSTITUTION.md) (maintenance cost).

## Overview

Every plant sample that enters this pipeline is filtered against a
reference "panel" — a FASTA of known organelle genomes that we map reads
to in order to decide *this read looks like the organelle we were asked
to assemble*. There is one panel per `assembly_target`.

The `plant_pt` (chloroplast) panel holds 101 genomes. The `animal_mt`
panel holds 1853. The `plant_mt` (plant mitochondrion) panel holds
**one**: a single *Vigna radiata* mitogenome, 398 kb. It came from
GetOrganelleDB's seed database, which ships exactly one seed for plant
mitochondria. GetOrganelle can afford that because its assembler *grows*
outward from the seed iteratively — the seed only has to get it started.
We use the panel differently: as a one-shot filter, and (since task 25)
also to decide which recruited reads count toward the coverage estimate
and which assembled contigs are mitochondrial rather than chloroplast.
For those three jobs, one distant genome is far too thin.

The symptom is measurable. On the `INT-PLANT-01-mt` fixture, genuine
mitochondrial contigs match the *Vigna* panel over only 28–32 % of their
length, because *Datura* and *Vigna* mitogenomes have diverged a long
way outside the conserved genes. Weak homology to the thing you are
looking for means weak recruitment, and it is a large part of why that
sample's recruited read pool came back 70 % chloroplast.

The obvious fix — swap in the 679-genome RefSeq Viridiplantae
mitochondrion set we *already ship* for `BLAST_VALIDATE` — turns out to
be actively dangerous, and understanding why is the whole point of this
task. Plant mitochondrial genomes really do contain chunks of
chloroplast DNA: over evolutionary time, plastid sequence gets copied
into the mitochondrial genome. These insertions are called **NUPTs**
(nuclear/organellar plastid DNA). Across the RefSeq collection they add
up to **26 % of the total sequence**, present in 624 of the 679 genomes.

So a broad, unmodified plant-mitochondrion panel is, in part, a
chloroplast panel. Feed it to the pipeline and chloroplast contigs start
matching the *mitochondrial* reference just as well as they match the
chloroplast one — the signal we use to tell the two organelles apart
collapses to a coin flip.

The fix is to keep the breadth and delete the confusion: take all 679
genomes, find the regions that align to the chloroplast panel, mask them
out (replace with `N`), and index what remains. That leaves 191 Mb of
sequence that is mitochondrial *and only* mitochondrial — still 480×
more reference than we have today.

## 1. Evidence

Measured 2026-08-03 while investigating
[task 25](completed/25_coverage_gate_carryover.md), on the
`INT-PLANT-01-mt` fixture assembly. Merged aligned fraction against each
candidate panel; the bracketed figure is the **margin over the
`plant_pt` panel**, which is the quantity C3 actually decides on
([§3.7.1](../spec/03-organelles.md#371-homology-is-measured-against-sibling-panels-too)).

| Contig | Truth | Vigna-only (current) | RefSeq raw | RefSeq pt-masked |
|---|---|---|---|---|
| contig_2 | mito | 0.278 (+0.226) | 0.999 (+0.946) | 0.996 (+0.944) |
| contig_5 | mito | 0.317 (+0.100) | 1.000 (+0.782) | 0.976 (+0.758) |
| contig_3 | plastid | 0.008 (−0.991) | 1.000 (**−0.000**) | 0.069 (−0.931) |
| contig_1 | plastid | 0.073 (−0.927) | 0.999 (**+0.000**) | 0.110 (−0.889) |

Read-level effect on the [§2.1.5](../spec/02-stages.md#215-sibling-organelle-carry-over-in-the-estimate)
coverage split for the same sample (400 kb nominal, 30× MIN):

| Panel | `estimated_cov` | Gate verdict |
|---|---|---|
| Vigna-only (current) | 27.29× | `low_coverage` |
| RefSeq raw | 52.88× | `ok` ← **false pass** |
| RefSeq pt-masked | 28.37× | `low_coverage` |

Three conclusions:

1. **Breadth is worth having.** Genuine mitochondrial homology rises
   from ~0.30 to ~0.98.
2. **Raw RefSeq must not be used.** It erases the sibling margin on
   plastid contigs and silently re-opens task 25's false `ok`. A panel
   change that reverts two landed defect fixes is worse than no change
   ([rule 18](../CONSTITUTION.md)).
3. **Masking preserves both properties**, and does *not* rescue the
   `INT-PLANT-01-mt` fixture (28.37×, still under MIN) — confirming that
   fixture is genuinely under-sequenced rather than badly referenced.
   The replacement-fixture item in [`tasks/todo.md`](todo.md) stands
   independently of this task.

Panel composition, for scale:

| Panel | Sequences | Bases |
|---|---|---|
| `plant_mt` current (GetOrganelle seed) | 1 | 398 kb |
| RefSeq Viridiplantae mt, raw | 679 | 259.5 Mb |
| — of which plastid-derived | 624 genomes affected | 68.5 Mb (26.4 %) |
| RefSeq Viridiplantae mt, masked | 679 | 191.0 Mb |

## 2. Is masking a `principle 5` violation?

**No, but say so explicitly in the PR.**
[Principle 5](../CONSTITUTION.md) requires that host exclusion be
*positive recruitment*, not *negative depletion*: we recruit reads
toward a target panel, and never discard a read for resembling a host
reference.

Masking operates on the **reference**, not on reads. No read is
excluded for matching anything; we are removing sequence from the bait
that was never mitochondrial-specific to begin with, so that "aligns to
the plant mitochondrion panel" means what it claims. Every read is still
judged solely on positive similarity to target sequence.

This is materially different from the RECRUIT-side competitive filter
held open at
[§9 item 1](../spec/07-open-questions.md#9-fine-tuning-post-prototype-benchmarking),
which *does* discard reads on a comparison and *does* need sign-off.
This task does not pre-empt that item; if anything it reduces the case
for it, and the §5 measurements should be fed back into it.

## 3. Build changes

### 3.1 Masking step

Extend [`scripts/refdata/build_recruit.sh`](../scripts/refdata/build_recruit.sh).
`plant_pt` and `animal_mt` are unchanged — they keep copying the
GetOrganelle seed FASTA. Only `plant_mt` gains a new derivation:

```
# pseudocode — plant_mt only
input:  refseq_mt_viridiplantae.fa   (from §4.2's build, see §3.3)
        plant_pt.fa                  (already staged this run)

1. minimap2 -x asm20 plant_pt.fa refseq_mt_viridiplantae.fa  -> PAF
2. per mitogenome, merge the aligned QUERY intervals
3. rewrite the FASTA replacing each merged interval with N
4. minimap2 -d plant_mt.mmi  <masked fasta>
```

Step 2 is the same merged-interval routine C3 uses
([§3.7.2](../spec/03-organelles.md#372-aligned-fraction-is-merged-across-all-alignment-blocks));
do not write a third copy of it — see §4.

`-x asm20` is the preset used in the prototype (divergent
genome-to-genome). Confirm it during implementation and record the
choice; `asm10` would mask less and `map-ont` is wrong for
assembly-to-assembly comparison.

### 3.2 Masking thresholds are config, not magic numbers

Per [rule 18](../CONSTITUTION.md) and the standing to-do on hard-coded
script config, expose as script arguments with recorded defaults:

| Knob | Suggested default | Why it matters |
|---|---|---|
| minimap2 preset | `asm20` | Sets how diverged a NUPT can be and still be masked |
| min alignment length to mask | 200 bp | Below this, masking removes conserved-gene fragments shared by both organelles rather than genuine NUPT |
| max masked fraction per genome | 0.60 | A mitogenome that is >60 % plastid by this measure is more likely a misannotated RefSeq record than a real genome — warn and skip it rather than emitting a mostly-`N` reference |

The third is a data-quality guard, not a tuning knob. Emit the per-genome
masked fraction to the build log and the top offenders to
`manifest.json` so the decision is auditable.

### 3.3 Build ordering

`refseq_mt_viridiplantae.fa` is produced in scratch by
[`scripts/refdata/build_validate.sh`](../scripts/refdata/build_validate.sh)
(line ~112) before being turned into a BLAST DB — it is **not** retained
in the bundle today, only the BLAST DB is. The recruit build now needs
it.

Decide between:

- **(a) Sequence the builds** — `build_validate.sh` runs first and
  retains its scratch FASTA for `build_recruit.sh` to consume. Cheapest,
  but couples two scripts that are currently independent and makes the
  order load-bearing and undocumented.
- **(b) Retain the FASTA in the bundle** — `build_validate.sh` writes
  `validate/refseq_mt_viridiplantae.fa` alongside the BLAST DB.
  Recommended: it is already downloaded and split, the bundle gains a
  reproducibility artefact rather than a derived-only DB, and
  `build_recruit.sh` gets an explicit input path instead of an implicit
  ordering dependency. Costs ~250 Mb in the bundle tarball — check that
  against the Azure blob size budget
  ([spec §4a](../spec/04a-azure-blob-storage.md)) before committing.
- **(c) Re-download in `build_recruit.sh`** — rejected. Duplicates a
  large download and risks the two scripts pinning different RefSeq
  releases, which would break the provenance guarantee of
  [rule 10](../CONSTITUTION.md).

### 3.4 Bundle version and provenance

This is a new bundle: **`refs/v2026.08/`**. The existing `v2026.07` is
immutable ([rule 10](../CONSTITUTION.md)) — do not rebuild in place.

`manifest.json` must record, for `plant_mt`:

- source (RefSeq release number + download date, inherited from §4.2),
- the masking parameters actually applied (§3.2),
- masked bases and masked fraction, panel-wide and worst-per-genome,
- any genome skipped by the §3.2 guard, with its fraction,
- SHA256 of both the masked FASTA and the `.mmi`.

An auditor reading a result six months from now must be able to tell
that this panel was masked, by what rule, and how much was removed.

Also update:

- `scripts/refs-v2026.08.sha256` and the Azure upload
  ([spec §4a](../spec/04a-azure-blob-storage.md));
- `conf/integration.config`'s `organelle_refs` path;
- any P3 dataset config pointing at `v2026.07`.

## 4. Where the masking code lives

The masking logic is Python (interval merge + FASTA rewrite), and it is
**pipeline-adjacent, not a workflow stage** — it belongs beside
[`scripts/refdata/split_refseq_mt.py`](../scripts/refdata/split_refseq_mt.py),
not in `bin/`. It is not a C1–C7 component and does not run at pipeline
runtime, so [rule 14](../CONSTITUTION.md)'s 100 % branch-coverage
requirement does not formally bind it.

It should still be tested, for a specific reason: a silent bug here
produces a *plausible but wrong reference bundle*, which then produces
plausible but wrong results in every downstream run — the exact class of
failure [rule 18](../CONSTITUTION.md) exists to prevent, and the hardest
to notice.

Minimum test surface, as a new `scripts/tests/test_mask_panel.py`:

| # | Scenario | Expected |
|---|---|---|
| 1 | Interval merge: overlapping, adjacent, nested, disjoint, single, empty | Correct merged length in each case |
| 2 | Single interval masked mid-sequence | Exactly that span becomes `N`; flanks byte-identical |
| 3 | Two overlapping intervals | Union masked once, no double-counting in the reported fraction |
| 4 | Interval at sequence start / end | No off-by-one; length preserved |
| 5 | Alignment shorter than the min-length threshold | Not masked |
| 6 | Genome above the max-masked-fraction guard | Skipped, warned, absent from output, recorded in the report |
| 7 | Genome with no plastid hits at all | Passed through byte-identical |
| 8 | Masked FASTA length == input length for every record | Masking must not change coordinates |

Case 8 is the invariant that matters most — masking substitutes, never
deletes.

**Reuse, do not re-implement, the interval merge.** It now exists twice
([`bin/bin_target.py`](../bin/bin_target.py) and
[`bin/coverage_gate.py`](../bin/coverage_gate.py), both added by tasks 23
and 25). A third copy in `scripts/refdata/` is the point at which this
becomes a maintenance liability ([rule 19](../CONSTITUTION.md)). Decide
during implementation whether to factor it into a shared module — note
that `bin/` is staged onto the container `PATH` at runtime while
`scripts/refdata/` is not, so a shared import is not free. Recording
"three copies, deliberately, because X" is an acceptable outcome; doing
it by accident is not.

## 5. Re-validation — the point of the task

A panel change silently alters the behaviour of three stages. None of
this is optional.

### 5.1 Re-measure the signals tasks 23 and 25 depend on

Against the new bundle, on all three integration fixtures **and** the P3
datasets:

- per-contig merged aligned fraction against declared and sibling panels
  (`bin_metadata.json` `panel_aligned_frac`), and the resulting margin;
- `sibling_carryover.sibling_organelle_fraction` per sample;
- `coverage_basis`, `estimated_cov`, `target_assigned_bases` and
  `sibling_organelle_fraction` from `sample_status.json`.

Record before/after in the outcomes section.

### 5.2 Regression guards that must not move

- **No plastid contig may out-score the plastid panel.** The margins in
  §1 must stay strongly negative for `contig_1`/`contig_3`-equivalents.
  If any sibling margin lands inside ±0.05, stop — that is the raw-RefSeq
  failure mode reappearing, and no threshold retune fixes it.
- **`INT-PLANT-01-pt` must not regress.** Its sibling is `plant_mt`;
  a broader mt panel could plausibly pull plastid reads across. It sat
  at `sibling_organelle_fraction` 0.0193 and 260.43× before. Assert it
  stays `ok` and the fraction stays under
  `max_sibling_organelle_fraction` in
  `tests/integration/expected/plant_pt/coverage_bounds.json`.
- **`INT-ANIMAL-01` must be byte-identical.** It has no sibling panel
  and does not touch `plant_mt`. Any change is a bug in the build.

### 5.3 Threshold review, not retune

[§3.7.6](../spec/03-organelles.md#376-revised-per-target-criteria)'s
`plant_mt` `min_aligned_frac` of **0.15** was set against the
single-genome panel, where genuine mitochondrial contigs scored
0.17–0.35. Against a masked broad panel they score ~0.98, so 0.15 stops
being a meaningful floor — it would admit almost anything that hits the
panel at all.

Raise it, and record the evidence for the new value.
[Task 23 §2.3](completed/23_bin_target_recalibration.md) put these
thresholds in `nextflow.config` precisely so this is a config change.
**Do not tune until the fixtures pass** — derive the value from the
observed separation between the mito and plastid populations across the
P3 datasets, then check the fixtures.

The same review applies to `min_identity` (70.0) and to
`params.coverage_limits.plant_mt.min_cov` — the last only if §5.1 shows
the corrected estimate has moved materially for reasons other than this
fixture being under-sequenced.

## 6. Integration fixtures

No new fixture is required; this task changes reference data, not test
data. But note the interaction with the two open fixture items in
[`tasks/todo.md`](todo.md):

- `INT-PLANT-01-mt` stays soft-failed at ~28× (§1), so
  `ASSEMBLING_SAMPLES` in
  [`tests/integration/assertions.sh`](../tests/integration/assertions.sh)
  is unchanged and the `plant_mt` binning assertions stay dark. **This
  task does not restore them** — only a deeper fixture does.
- `tests/integration/expected/plant_mt/coverage_bounds.json` bounds
  `estimated_cov` to [15, 30]. The masked panel gives 28.37×, inside
  that band, so it should pass unchanged — but confirm rather than
  assume, and if it lands near 30 widen the *reasoning*, not just the
  bound.

If the replacement `plant_mt` fixture lands first, this task's §5.1
measurements should be taken on it too, and it becomes the primary
evidence for §5.3.

## 7. Decisions required

1. **§3.3 build ordering** — (a), (b) or (c). Recommendation (b), subject
   to the bundle-size check.
2. **§3.2 masking parameters** — confirm preset, min alignment length,
   and the max-masked-fraction guard, or supply alternatives.
3. **§4 interval-merge duplication** — factor into a shared module, or
   accept a third copy with the reason recorded.
4. **§5.3 new `min_aligned_frac`** — cannot be chosen before §5.1 runs;
   flag the value and its evidence for review rather than landing it
   silently.
5. **Sequencing against [task 26](26_binning_marker_genes.md)** — see
   the header. Recommendation: 28 first.

## 8. Exit criteria

- `refs/v2026.08/recruit/plant_mt.{fa,mmi}` built from masked RefSeq
  Viridiplantae mitogenomes; `plant_pt` and `animal_mt` unchanged from
  `v2026.07`.
- `manifest.json` records masking parameters, masked fraction, skipped
  genomes, and SHA256 for every artefact (§3.4).
- Bundle uploaded, `scripts/refs-v2026.08.sha256` committed,
  `conf/integration.config` repointed.
- `scripts/tests/test_mask_panel.py` covering §4 cases 1–8;
  `scripts/pytest.sh` no worse than before; `flake8` clean at 79 columns.
- §5.1 before/after table recorded in this task's outcomes.
- §5.2 regression guards all hold — in particular **no sibling margin
  inside ±0.05**.
- `-profile integration` green with `bash tests/integration/assertions.sh`
  passing; `-stub-run` green.
- [§4.1](../spec/04-reference-data.md#41-setup-task-source-reference-panel-from-getorganelle)
  updated: the `plant_mt` row no longer says "GetOrganelleDB
  `embplant_mt`", and the "re-evaluate after P1" clause is resolved with
  the §1 evidence.
- [§9 item 1](../spec/07-open-questions.md#9-fine-tuning-post-prototype-benchmarking)
  updated with §5.1's effect on recruitment — a better panel may reduce
  the case for a RECRUIT-side competitive filter.
- [Task 26](26_binning_marker_genes.md)'s §2 gate re-read against the new
  measurements, and its status recorded (proceed / close).
