# Task 28 — Plastid-masked `plant_mt` recruitment panel

> **Reconcile with [task 26](../26_binning_marker_genes.md) before starting
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
>   and the cheaper outcome under [rule 19](../../CONSTITUTION.md).
> - If task 26 lands first, its marker-gene criterion would be tuned
>   against a panel this task then replaces, and would need
>   recalibrating.
>
> **Recommended sequence: 28 before 26.** Task 26's §2 should consume
> this task's post-change measurements. Whoever picks up either task
> first must read both and record the sequencing decision.

**Phase:** P-tune-C ([spec §9 grouping](../../spec/07-open-questions.md#9-fine-tuning-post-prototype-benchmarking)).
Reference-data change; no workflow logic moves.
**Goal:** Replace the single-sequence `plant_mt` recruitment panel with a
broad, plastid-masked panel built from RefSeq, without re-opening the
sibling-organelle defects that
[task 23](23_bin_target_recalibration.md) and
[task 25](25_coverage_gate_carryover.md) closed.

**Prerequisite:** [task 25](25_coverage_gate_carryover.md)
landed (C2 gates on target-assigned bases; C3 records
`sibling_carryover`). This task's success metric is expressed in fields
those two tasks introduced.

**Spec basis:**
[§4.1](../../spec/04-reference-data.md#41-setup-task-source-reference-panel-from-getorganelle)
(the panel source, and its own "re-evaluate after P1" clause),
[§4.2](../../spec/04-reference-data.md#42-refseq-organelle-dbs-for-blast_validate)
(the RefSeq Viridiplantae mitochondrion set this reuses),
[§4.4](../../spec/04-reference-data.md#44-consolidated-build-script)
(bundle layout + `manifest.json`),
[§2.1.5](../../spec/02-stages.md#215-sibling-organelle-carry-over-in-the-estimate)
(the coverage split that consumes the panel),
[§3.7.1](../../spec/03-organelles.md#371-homology-is-measured-against-sibling-panels-too)
(the sibling discrimination that consumes it),
[rule 10](../../CONSTITUTION.md) (versioned, immutable, provenance-tracked
bundles), [rule 19](../../CONSTITUTION.md) (maintenance cost).

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
[task 25](25_coverage_gate_carryover.md), on the
`INT-PLANT-01-mt` fixture assembly. Merged aligned fraction against each
candidate panel; the bracketed figure is the **margin over the
`plant_pt` panel**, which is the quantity C3 actually decides on
([§3.7.1](../../spec/03-organelles.md#371-homology-is-measured-against-sibling-panels-too)).

| Contig | Truth | Vigna-only (current) | RefSeq raw | RefSeq pt-masked |
|---|---|---|---|---|
| contig_2 | mito | 0.278 (+0.226) | 0.999 (+0.946) | 0.996 (+0.944) |
| contig_5 | mito | 0.317 (+0.100) | 1.000 (+0.782) | 0.976 (+0.758) |
| contig_3 | plastid | 0.008 (−0.991) | 1.000 (**−0.000**) | 0.069 (−0.931) |
| contig_1 | plastid | 0.073 (−0.927) | 0.999 (**+0.000**) | 0.110 (−0.889) |

Read-level effect on the [§2.1.5](../../spec/02-stages.md#215-sibling-organelle-carry-over-in-the-estimate)
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
   ([rule 18](../../CONSTITUTION.md)).
3. **Masking preserves both properties**, and does *not* rescue the
   `INT-PLANT-01-mt` fixture (28.37×, still under MIN) — confirming that
   fixture is genuinely under-sequenced rather than badly referenced.
   The replacement-fixture item in [`tasks/todo.md`](../todo.md) stands
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
[Principle 5](../../CONSTITUTION.md) requires that host exclusion be
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
[§9 item 1](../../spec/07-open-questions.md#9-fine-tuning-post-prototype-benchmarking),
which *does* discard reads on a comparison and *does* need sign-off.
This task does not pre-empt that item; if anything it reduces the case
for it, and the §5 measurements should be fed back into it.

## 3. Build changes

### 3.1 Masking step

Extend [`scripts/refdata/build_recruit.sh`](../../scripts/refdata/build_recruit.sh).
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
([§3.7.2](../../spec/03-organelles.md#372-aligned-fraction-is-merged-across-all-alignment-blocks));
do not write a third copy of it — see §4.

`-x asm20` is the preset used in the prototype (divergent
genome-to-genome). Confirm it during implementation and record the
choice; `asm10` would mask less and `map-ont` is wrong for
assembly-to-assembly comparison.

### 3.2 Masking thresholds are config, not magic numbers

Per [rule 18](../../CONSTITUTION.md) and the standing to-do on hard-coded
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
[`scripts/refdata/build_validate.sh`](../../scripts/refdata/build_validate.sh)
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
  ([spec §4a](../../spec/04a-azure-blob-storage.md)) before committing.
- **(c) Re-download in `build_recruit.sh`** — rejected. Duplicates a
  large download and risks the two scripts pinning different RefSeq
  releases, which would break the provenance guarantee of
  [rule 10](../../CONSTITUTION.md).

### 3.4 Bundle version and provenance

This is a new bundle: **`refs/v2026.08/`**. The existing `v2026.07` is
immutable ([rule 10](../../CONSTITUTION.md)) — do not rebuild in place.

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
  ([spec §4a](../../spec/04a-azure-blob-storage.md));
- `conf/integration.config`'s `organelle_refs` path;
- any P3 dataset config pointing at `v2026.07`.

## 4. Where the masking code lives

The masking logic is Python (interval merge + FASTA rewrite), and it is
**pipeline-adjacent, not a workflow stage** — it belongs beside
[`scripts/refdata/split_refseq_mt.py`](../../scripts/refdata/split_refseq_mt.py),
not in `bin/`. It is not a C1–C7 component and does not run at pipeline
runtime, so [rule 14](../../CONSTITUTION.md)'s 100 % branch-coverage
requirement does not formally bind it.

It should still be tested, for a specific reason: a silent bug here
produces a *plausible but wrong reference bundle*, which then produces
plausible but wrong results in every downstream run — the exact class of
failure [rule 18](../../CONSTITUTION.md) exists to prevent, and the hardest
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
([`bin/bin_target.py`](../../bin/bin_target.py) and
[`bin/coverage_gate.py`](../../bin/coverage_gate.py), both added by tasks 23
and 25). A third copy in `scripts/refdata/` is the point at which this
becomes a maintenance liability ([rule 19](../../CONSTITUTION.md)). Decide
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

[§3.7.6](../../spec/03-organelles.md#376-revised-per-target-criteria)'s
`plant_mt` `min_aligned_frac` of **0.15** was set against the
single-genome panel, where genuine mitochondrial contigs scored
0.17–0.35. Against a masked broad panel they score ~0.98, so 0.15 stops
being a meaningful floor — it would admit almost anything that hits the
panel at all.

Raise it, and record the evidence for the new value.
[Task 23 §2.3](23_bin_target_recalibration.md) put these
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
[`tasks/todo.md`](../todo.md):

- `INT-PLANT-01-mt` stays soft-failed at ~28× (§1), so
  `ASSEMBLING_SAMPLES` in
  [`tests/integration/assertions.sh`](../../tests/integration/assertions.sh)
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
5. **Sequencing against [task 26](../26_binning_marker_genes.md)** — see
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
- [§4.1](../../spec/04-reference-data.md#41-setup-task-source-reference-panel-from-getorganelle)
  updated: the `plant_mt` row no longer says "GetOrganelleDB
  `embplant_mt`", and the "re-evaluate after P1" clause is resolved with
  the §1 evidence.
- [§9 item 1](../../spec/07-open-questions.md#9-fine-tuning-post-prototype-benchmarking)
  updated with §5.1's effect on recruitment — a better panel may reduce
  the case for a RECRUIT-side competitive filter.
- [Task 26](../26_binning_marker_genes.md)'s §2 gate re-read against the new
  measurements, and its status recorded (proceed / close).

## 9. Outcomes

Executed 2026-08-03. Bundle **`refs/v2026.08/`** built, integration green,
all §5.2 regression guards hold. One material regression found and
**deliberately not fixed here** — see §9.5.

### 9.1 Decisions (§7)

| # | Decision | Resolution |
|---|---|---|
| 1 | §3.3 build ordering | **(b)** — `build_validate.sh` retains `validate/refseq_mt_viridiplantae.fa` (262 MB); `build_recruit.sh` takes it as an explicit `--refseq-mt` path. Bundle grew 865 MB → **1.9 GB**, larger than the task's ~250 MB estimate because `plant_mt.mmi` is now 536 MB (was 1.5 MB). Accepted after review. |
| 2 | §3.2 masking parameters | Confirmed as suggested: `asm20`, min alignment 200 bp, max masked fraction 0.60. `asm20` was tested against `map-ont` and is clearly correct — map-ont masks **6.4 %** vs asm20's **19.7 %**, i.e. it is markedly *less* sensitive for genome-to-genome divergent homology. |
| 3 | §4 interval-merge duplication | **Factored out.** New [`bin/intervals.py`](../../bin/intervals.py) holds `merge_intervals` + `merged_length`; `bin_target.py` and `coverage_gate.py` import it as a bin/ sibling (the mechanism `bin_target.py` already used for `plastid_canonicalise`), and `mask_panel.py` runs on the host and `sys.path`-inserts `bin/`. Three copies became one. |
| 4 | §5.3 new `min_aligned_frac` | **Flagged, not landed** — see §9.4. |
| 5 | Sequencing vs task 26 | 28 first, as recommended. Task 26 re-read; stays open — see §9.6. |

### 9.2 Panel build

RefSeq release **236**, downloaded 2026-07-26 (the `v2026.07` scratch was
still on disk, so no re-download; the BLAST DBs and `proteins/` were
carried over unchanged from `v2026.07`, which is legitimate because the
RefSeq release did not move).

| | Value |
|---|---|
| Genomes in | 679 |
| Genomes written | 652 |
| Genomes skipped by the 0.60 guard | **27** (max 76.5 % plastid; the guard is doing real work) |
| Genomes with any plastid hit | 606 |
| Bases written | 249,343,119 |
| Bases masked | 49,071,472 (**19.7 %**) |

Masked fraction came in at 19.7 %, below §1's exploratory 26.4 %, because
the 200 bp floor and the 0.60 guard both remove sequence from the count
that the exploratory pass included.

`manifest.json` (new — spec §4.4 asked for one and no bundle had ever had
one) records RefSeq provenance, the masking parameters actually applied,
panel-wide and per-genome masked fractions, every skipped genome, and
SHA256 + byte count for all 54 artefacts. Generated by the new
[`scripts/refdata/build_manifest.py`](../../scripts/refdata/build_manifest.py).

### 9.3 §5.1 before/after

**Per-contig merged aligned fraction**, `INT-PLANT-01-mt` assembly, measured
with C3's own instrument (`mappy`, `preset='map-ont'`, `best_n=5`) rather
than the `asm20` used for §1's exploratory table. Bracketed figure is the
margin over `plant_pt`.

| Contig | Truth | v2026.07 (Vigna) | v2026.08 (masked) |
|---|---|---|---|
| contig_2 | mito | 0.200 (+0.143) | 0.514 (**+0.458**) |
| contig_5 | mito | 0.335 (+0.220) | 0.729 (**+0.614**) |
| contig_1 | plastid | 0.064 (−0.936) | 0.347 (**−0.653**) |
| contig_3 | plastid | 0.012 (−0.988) | 0.314 (**−0.686**) |

**§1's predicted margins (−0.889 to −0.931) do not reproduce under the
production metric.** §1 measured at `asm20`, where I reproduce it closely
(−0.920 / −0.910); C3 measures at `map-ont`/`best_n=5`, which is far more
permissive and lets plastid contigs pick up 31–35 % against 652 diverse
mitogenomes. The margin still widens for genuine mito contigs and stays
strongly negative for plastid ones, so the change is directionally right —
but the headroom is roughly a third of what §1 advertised. Probable cause
and a proposed fix are filed in [`tasks/todo.md`](../todo.md) (mask against the
full RefSeq plastid set, not the 101-genome `plant_pt` panel).

**Sample-level**, from `sample_status.json`:

| Sample | Metric | v2026.07 | v2026.08 |
|---|---|---|---|
| INT-ANIMAL-01 | everything | — | **byte-identical** ✓ |
| INT-PLANT-01-pt | `estimated_cov` | 260.43× | 246.94× |
| | `sibling_organelle_fraction` | 0.0193 | **0.0701** |
| | recruited reads | 2513 | 2513 (unchanged) |
| INT-PLANT-01-mt | `estimated_cov` | 27.29× | **21.01×** |
| | `sibling_organelle_fraction` | 0.7053 | **0.7485** |
| | `target_assigned_bases` | 10,916,501 | 8,403,865 |
| | `target_merged_aligned_bases` | 726,531 | **2,230,195** |
| | recruited reads | 1366 | **1253** |

`INT-PLANT-01-mt` landed at 21.01×, not §1's predicted 28.37× — same
root cause as §9.5.

### 9.4 §5.3 threshold review — flagged, not landed

`params.bin_target_thresholds.plant_mt.min_aligned_frac` stays at **0.15**.
Against the masked panel the observed populations are mito 0.514/0.729 and
plastid 0.314/0.347, so 0.15 is indeed no longer a meaningful floor and
something near **0.45** would separate them. It is not landed because that
value rests on two contigs of each class from a single fixture that does
not even reach BIN_TARGET — exactly the tuning §5.3 warns against ("derive
the value from the observed separation ... across the P3 datasets"). The
sibling test (criterion 3) separates these populations by >1.1 on its own,
so nothing is currently mis-binned by leaving the floor low. Recorded for
[§9 item 10](../../spec/07-open-questions.md#9-fine-tuning-post-prototype-benchmarking).

`min_identity` (70.0) and `coverage_limits.plant_mt.min_cov` unchanged.

### 9.5 Regression found: MAPQ collapse in RECRUIT

**The broader panel makes `INT-PLANT-01-mt` recruitment worse, not better.**
Recruited reads fell 1366 → 1253 and plastid carry-over rose 70.5 % → 74.9 %.

Measured cause: [`modules/local/recruit.nf`](../../modules/local/recruit.nf)
filters with `samtools view -F 4 -q 1`. Against one reference, every mapped
read cleared it. Against 652 near-redundant references:

| Panel | Primary alignments | Pass `-q 1` | Dropped at MAPQ 0 |
|---|---|---|---|
| v2026.07 (1 genome) | 1366 | 1366 | 0 |
| v2026.08 (652 genomes) | 1352 | 1074 | **278** |

`-q 1` is meant to mean "confidently placed"; on a redundant panel it means
"found in only one genome", which is the opposite of what breadth is for.

**Is it repairable? Probably yes — the discarded reads are the good ones.**
Composition of the two pools, each read assigned by C2's own rule (merged
aligned fraction vs `plant_mt` > vs `plant_pt`):

| Pool | Reads | Mito-assigned | Plastid-assigned |
|---|---|---|---|
| MAPQ 0 — discarded by `-q 1` | 278 | **150 (54.0 %)**, 5.05 Mb | 128 (46.0 %), 3.81 Mb |
| MAPQ ≥ 1 — kept | 1074 | 252 (**23.5 %**), 6.34 Mb | 822 (76.5 %), 21.40 Mb |

The discarded pool is **2.3× richer in mitochondrial reads than the pool
that survives the filter**. `-q 1` is not accidentally protective here; it
is preferentially deleting the signal the panel exists to find. That is the
best possible diagnosis for repairability — the fix is a filter change, not
a panel redesign.

Three candidate repairs, in increasing cost:

1. **Drop the MAPQ term** (`-F 4` alone). Recovers 5.05 Mb mito at the cost
   of 3.81 Mb plastid. `estimated_cov = target_assigned_bases / nominal`,
   so this fixture would move from 8.40 Mb / 21.01× to roughly 13.45 Mb /
   **33.6×**, and carry-over from 74.9 % to ~68 %. **Handle with care:** that
   crosses the 30× MIN and would flip `INT-PLANT-01-mt` back to `ok` — the
   exact false pass [task 25](25_coverage_gate_carryover.md) closed. Before
   adopting, confirm the recovered bases are real depth and not
   multi-mapping inflation.
2. **Replace MAPQ with a panel-size-independent length criterion** (the
   ptGAUL-style filter already named in §9 item 1). Measured on the same
   reads: merged aligned fraction ≥ 0.7 gives **95.9 % mitochondrial
   purity** but only 196 reads / 2.11 Mb; ≥ 0.5 gives 75.2 % / 2.63 Mb.
   Excellent purity, far too lossy to be the sole criterion at this depth.
   Likely useful in combination with (1).
3. **Reduce panel redundancy** — cluster the 652 genomes to a
   non-redundant representative set. Attacks the cause rather than the
   symptom and keeps MAPQ meaningful, but costs a build step and some of
   the breadth this task was for.

Caveat: this is one fixture, and the pathological one. The direction is
unambiguous; the magnitudes are not general.

**Deliberately not fixed in this task.** Task 28's own framing is
"reference-data change; no workflow logic moves", and RECRUIT filter
strictness is
[§9 item 1](../../spec/07-open-questions.md#9-fine-tuning-post-prototype-benchmarking),
which requires its own sign-off. §9 item 1 has been updated with this
measurement and the MAPQ term promoted to the urgent half of that item.
Homology *detection* did improve as intended — `target_merged_aligned_bases`
tripled, 726 kb → 2.23 Mb — so the panel is doing its job; the loss is
entirely in the filter downstream of it.

### 9.6 Task 26 status

**Re-read; stays open.** The new margins (mito +0.46/+0.61 vs plastid
−0.65/−0.69, a gap over 1.1 with no overlap) strengthen task 26's own §2
preference to close as unnecessary. They do not settle it: §2 requires
stability and margin measurements on the P3 datasets, and question 3 needs
a sample with known NUMT content. Neither exists, and `INT-PLANT-01-mt`
still soft-fails the gate so it cannot contribute. Task 26's premise is
unchanged by this task except that its close case is now better supported.

### 9.7 Verification

- `scripts/pytest.sh` — 143 passed, **100 % branch coverage** across all of
  `bin/` including the new `intervals.py`.
- `flake8` clean at 79 columns on every touched Python file.
- `-profile integration` green (40 tasks, 12m10s);
  `bash tests/integration/assertions.sh` — **all assertions passed**.
- `-profile stub,docker -stub-run` green.
- §5.2 guards: no sibling margin anywhere inside ±0.05 (closest is −0.653);
  `INT-PLANT-01-pt` stayed `ok` with sibling fraction 0.0701 < 0.15;
  `INT-ANIMAL-01` byte-identical across assembly, target, bin metadata,
  coverage status and barcodes.
- No fixture bound was retuned.

### 9.8 Follow-ups filed (superseded in part by §10)

In [`tasks/todo.md`](../todo.md), new "Reference data" section:

- mask against the full RefSeq plastid set rather than the 101-genome
  `plant_pt` panel (with the caveat that a first attempt was OOM-killed);
- confirm the 536 MB `plant_mt.mmi` fits CI and Azure memory budgets.

---

## 10. Follow-up — RECRUIT MAPQ repair (executed 2026-08-04)

Executed on request after §9.5, as a follow-up within this task rather
than a separate brief. **This moves workflow logic**, which §8's framing
("reference-data change; no workflow logic moves") excluded — the
sign-off [§9 item 1](../../spec/07-open-questions.md#9-fine-tuning-post-prototype-benchmarking)
requires for RECRUIT filter strictness is this instruction.

Two repairs were requested:

1. drop the MAPQ term;
2. replace it with a panel-size-independent length criterion.

**Repair 1 landed. Repair 2 is implemented but shipped inert, because
measurement showed every non-zero threshold breaks a fixture bound.**

### 10.1 Why MAPQ had to go, restated from the full distribution

`-q 1` discards MAPQ 0. Composition of the whole `INT-PLANT-01-mt`
primary-alignment set against the v2026.08 panel, each read assigned by
C2's rule:

| MAPQ | reads | mito | plastid | % mito |
|---|---|---|---|---|
| 0 | 278 | 150 | 128 | **54.0** |
| 1–9 | 188 | 133 | 55 | **70.7** |
| 10–29 | 67 | 36 | 31 | 53.7 |
| 30–59 | 175 | 27 | 148 | 15.4 |
| 60 (max) | 644 | 56 | 588 | **8.7** |

MAPQ is **anti-correlated with being the declared organelle**. Mean
merged aligned fraction is flat across the buckets (0.282 at MAPQ 0,
0.329 at MAPQ 60), confirming MAPQ is tracking placement uniqueness and
not match extent.

The mechanism: minimap2 derives MAPQ from the gap between best and
second-best chaining score. Against one reference, "second best" is
another locus in the same genome — a repeat — and discarding those is
ordinary hygiene. Against 652 near-identical mitogenomes, "second best"
is *the same locus in another species*. A read matching conserved
mitochondrial sequence present in hundreds of panel genomes has no score
gap and scores 0; a plastid read matching residual unmasked NUPT in one
genome is unique and scores 60. Recruitment asks "is this read
organellar?", not "which species' organelle is it?" — on a
single-reference panel those questions shared an answer, and on a broad
panel they invert.

### 10.2 Why repair 2 could not be landed hot

Threshold sweep, measured per fixture against its declared panel, with
C2's target/sibling split applied so the numbers are the coverage gate's
actual inputs. `cov` is `target_assigned_bases / nominal_size`.

| min_frac | min_bp | ANIMAL cov | PLANT-pt cov | PLANT-mt cov |
|---|---|---|---|---|
| **0.0** | **0** | **327.52** | **254.13** | **28.48** |
| 0.0 | 1000 | 162.37 | 169.47 | 10.24 |
| 0.1 | 0 | 149.51 | 154.56 | 7.34 |
| 0.2 | 0 | 129.05 | 148.95 | 6.35 |
| 0.3 | 0 | 118.95 | 147.16 | 6.10 |
| 0.5 | 0 | 113.17 | 145.49 | 5.35 |
| 0.7 | 0 | 104.95 | 143.55 | 5.20 |

Expected bounds are animal_mt [200, 500], plant_pt [150, 400],
plant_mt [15, 30]. **Only the 0.0 / 0 row sits inside all three.**

- `animal_mt` leaves its band immediately (327 → 150 at min_frac 0.1).
- `plant_pt` survives only to 0.1 (154.56, against a 150 floor).
- `plant_mt` collapses catastrophically: 28.48× → 7.34× at 0.1, and
  → 10.24× merely by requiring 1 kb aligned.

The `plant_mt` collapse is the substantive result. Its genuine reads are
long ONT reads that match the panel over only a *portion* of their
length, so a fractional floor removes real mitochondrial depth far
faster than it removes plastid. Purity does improve — frac ≥ 0.7 yields
95.9 % mitochondrial — but at 5.20×, which is not a usable sample.
Purity was never the objective; recovering genuine target depth is.

Per this repo's own rule, a criterion that moves fixtures outside their
bounds is **a finding, not a bound to retune**. So repair 2 ships as
`params.recruit_thresholds`, per-target, defaulting to 0/0 — a genuine
no-op today. The knobs exist so §9 item 1's benchmark is a config sweep
rather than a code change, exactly as task 23 did for
`bin_target_thresholds`.

**Read 0/0 as "not yet determined", not as "measured correct."** The
sweep above establishes only that the *fixtures* reject every non-zero
value, and those fixture bounds derive from the same under-powered data
— subsampled to a CI budget, pre-R10.4.1 chemistry, one accession
serving both plant targets. Tuning until the fixtures pass is circular
reasoning, which is now written up as
[spec §5.1](../../spec/05-test-data.md#51-the-integration-fixtures-are-correctness-tests-not-benchmarks).
A real value has to come from the P3 datasets. The purity data in the
table above is still useful input to that sweep; the *thresholds* it
implies are not.

### 10.3 Changes

- **[`bin/recruit_filter.py`](../../bin/recruit_filter.py)** (new) — reads
  the PAF, merges query intervals per read via
  [`bin/intervals.py`](../../bin/intervals.py), applies both floors,
  emits a read-ID list plus a `recruit_stats.json` sidecar.
- **[`modules/local/recruit.nf`](../../modules/local/recruit.nf)** —
  minimap2 now emits PAF (`-x map-ont`) instead of SAM; the
  `samtools view -F 4 -q 1 | cut -f1 | sort -u` pipeline is replaced by
  `recruit_filter.py`. `recruit_stats.json` added to `output:` and to the
  stub block.
- **[`nextflow.config`](../../nextflow.config)** —
  `params.recruit_thresholds`, per target, both floors 0 (rule 18).
- **[`conf/containers.config`](../../conf/containers.config)** — RECRUIT
  moves to the **already-SHA-pinned** coverage-gate image
  (`python3.12_seqkit2.13_seqtk1.5_minimap2-2.31`). PAF-based filtering
  removes the samtools dependency, and that image already carries
  python3 + minimap2 + seqtk, so no new container build and no rule 12
  exposure.
- **[`scripts/tests/test_recruit_filter.py`](../../scripts/tests/test_recruit_filter.py)**
  (new) — 14 cases, 100 % branch coverage of the selection path.

Note the interval-merge consolidation from §9.1 decision 3 paid off
immediately: `recruit_filter.py` is the fourth consumer of
`merge_intervals` and reuses it rather than adding an awk reimplementation
inside the module, which was the alternative given the old recruit image
had no python3.

### 10.5 Measured effect

Three states, so the panel change and the filter repair can be told
apart:

- **A** — `v2026.07` panel + `-q 1` (the pre-task-28 baseline)
- **B** — `v2026.08` masked panel + `-q 1` (task 28 §9 as shipped)
- **C** — `v2026.08` masked panel + merged-extent filter (this repair)

| Sample | Metric | A | B | **C** |
|---|---|---|---|---|
| INT-ANIMAL-01 | `estimated_cov` | 325.81 | 325.81 | **327.52** |
| | recruited reads | 887 | 887 | **889** |
| INT-PLANT-01-pt | `estimated_cov` | 260.43 | 246.94 | **254.13** |
| | `sibling_organelle_fraction` | 0.0193 | 0.0701 | **0.0822** |
| | recruited reads | 2513 | 2513 | **2605** |
| INT-PLANT-01-mt | `estimated_cov` | 27.29 | 21.01 | **28.48** |
| | `sibling_organelle_fraction` | 0.7053 | 0.7485 | **0.6887** |
| | `target_assigned_bases` | 10,916,501 | 8,403,865 | **11,390,928** |
| | recruited reads | 1366 | 1253 | **1352** |

Measured values match §10.2's predicted sweep row exactly (28.48 /
254.13 / 327.52), so the offline model of the filter is sound.

**The headline is B → C on `plant_mt`.** Coverage recovers 21.01× →
28.48× and plastid carry-over falls 0.7485 → 0.6887. Both are now
*better than the pre-task-28 baseline* (27.29×, 0.7053), which B was
not. The masked panel only pays for itself once the filter stops
selecting against it — §9.5's regression is fully reversed, and the
breadth this task bought is finally visible in the numbers.

`INT-PLANT-01-mt` remains `low_coverage` at 28.48× against the 30× MIN,
so [task 25](25_coverage_gate_carryover.md)'s verdict stands and no
false `ok` is re-opened. The fixture is genuinely under-sequenced; the
replacement item in [`tasks/todo.md`](../todo.md) is unaffected.

One benign observation: `INT-ANIMAL-01`'s assembly moved 348,580 →
215,078 bp. That sample exceeds `max_cov` and is subsampled, so a
slightly different recruited pool yields a different subsample and a
different assembly. The binned target is unchanged in substance —
16,799 bp, single contig, circular via `flye_circ` — and both figures
sit inside `assembly_bounds.json`. Not a regression, but it means
`animal_mt` assembly length is not a stable quantity across recruitment
changes.

### 10.6 Verification (post-change)

- `-profile integration` green (40 tasks, 9m58s);
  `bash tests/integration/assertions.sh` — **all assertions passed**,
  every fixture inside its existing bounds. **No bound was retuned.**
- `-profile stub,docker -stub-run` green.
- `scripts/pytest.sh` — 157 passed, 100 % branch coverage across `bin/`.
- `flake8` clean at 79 columns.
- `recruit_stats.json` published per sample, recording reads aligned vs
  recruited and the thresholds applied (rule 16 provenance). With the
  shipped 0/0 defaults, recruited == aligned on all three fixtures,
  which is the intended no-op.
