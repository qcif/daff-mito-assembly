# Task 68 — `animal_mt` recruitment inflated by microsatellites (stub)

**Phase:** P5, real-data tuning (spec/07-open-questions.md §9).

**Status:** scoped 2026-10-16 (see Scoped plan); awaiting go-ahead. Follows on
from task 58 (tasks/completed/58_animal_mt_bin_identity_floor.md, §5
Out of scope). Its §3.3 round 2/3 re-run (2026-10-15, `animal_mt`
identity floor 70%) is the baseline:

- `tests/manual/output/round_3_t58/` and
  `tests/manual/output/round_2_insect_t58/` (gitignored; also on the
  remote as `/mnt/data/wf5/output_round_{3,2_insect}_t58/`);
- `tests/integration/t58/output/` for `INT-ANIMAL-01` (gate 327.52×).

The baseline has gate estimates, statuses, bin selections and barcode
loci, but not recruited read counts or self-mapped depth. The remote
work dirs are deleted and intermediates were not published, so derive
self-mapped depth by mapping raw reads (`/mnt/data/wf5/data/`) back to
each sample's selected contigs as the first step.

## Overview

`RECRUIT` (stage 5) keeps any read that aligns anywhere to the target
organelle's reference panel. The panel is RefSeq mitogenomes, and some
of them carry simple-sequence repeats (microsatellites such as (CA)n
or (GA)n) in their control regions. Nuclear genomes are full of the
same repeats. A nuclear read that happens to contain a (CA)n run aligns
to that stretch of the panel and is recruited as "mitochondrial".

On client round 3, sample barcode20 (Cicadellidae), 78 % of recruited
reads came in on alignments under 300 bp to repeats like the (CA)n
tail of the *Indri indri* reference and the (GA)n of *Lepidodermella
squamata*. The coverage gate (stage 6) estimated 391× and passed the
sample on to assembly. Flye produced 94 contigs, mostly nuclear.
Mapping the reads back to the sample's own mitochondrial contigs puts
the real depth at low single digits up to about 20×.

That is two problems. The coverage gate's figure, which the report
shows and which decides `fail` / `low_coverage` / `ok`, is badly
overstated. And the assembler is handed mostly nuclear reads, which
costs time and raises the NUMT risk in binning.

`params.recruit_thresholds` already has floors for this
(`min_aligned_bp`, `min_aligned_frac`), but both ship at 0. Task 28
found that any non-zero value moves an integration fixture out of its
expected bounds. Task 58 measured that a 1000 bp floor alone still
leaves about 65× on barcode20, which is still overstated.

## Candidate approaches (to evaluate, not decided)

1. **Mask low-complexity sequence in the `animal_mt` recruit panel at
   build time,** in the same way the `plant_mt` panel has its
   plastid-derived regions masked (task 28, `scripts/refdata/mask_panel.py`).
   A tandem-repeat or low-complexity masker (for example `dustmasker`
   from the already-pinned BLAST image, or TRF) turns the repeats into
   N, so nothing aligns to them. This is a new bundle version with
   masking provenance in the manifest. It looks the most direct.
2. **A recruitment floor** (`min_aligned_bp`) on its own, or as well as
   (1). Measured to be insufficient alone.
3. **A complexity filter on the alignment** in `bin/recruit_filter.py`:
   drop hits whose aligned reference segment is mostly low-complexity.
   This is a runtime cost on every read.

## Scope to decide when this is picked up

- Evidence set: round 2 and round 3 `animal_mt` samples
  (`tests/manual/samples_round_2_insect.csv`,
  `tests/manual/samples_round_3.csv`), plus `INT-ANIMAL-01`. For each:
  the recruited read count, the gate estimate, and self-mapped
  mitochondrial depth, before and after.
- Acceptance: barcode20's gate estimate within a stated factor of its
  self-mapped depth. Barcode recovery on round 2 barcode05/06/13 and
  round 3 barcode13 unchanged. `INT-ANIMAL-01` stays within its
  coverage and assembly bounds, or a bounds change is justified as
  biology (spec §5.1: fixtures are correctness tests, not benchmarks).
- Check whether `plant_pt` and `plant_mt` panels have the same
  exposure (plant mitogenomes are repeat-rich). Masking all three in
  one bundle may be cheaper than three bundles.
- Unit tests for any `bin/` or `scripts/refdata/` change. Provenance
  for any masking (CONSTITUTION principle 10).

## Related tasks/todo.md items (fold in when scoped)

- Recruitment: "`params.recruit_thresholds` ship at 0/0".
- Recruitment: "`INT-ANIMAL-01`'s assembly length is not stable across
  recruitment changes". A recruitment change will move it again.

## Scoped plan (2026-10-16)

Decisions (user): approach 1, mask the `animal_mt` panel at build time;
measure-only for `plant_pt` / `plant_mt`; real-data runs on the
daff-admin remote (remote-run skill).

### Feasibility probe

`dustmasker` 1.0.0 is in the pinned BLAST 2.17.0 image. Run over
`refs/v2026.09_1/recruit/animal_mt.fa` (1853 records, 30.2 Mb):

| `-level` | masked bases | fraction |
|---|---|---|
| 20 (default) | 4.13 Mb | 13.7 % |
| 30 | 1.42 Mb | 4.7 % |

Default DUST probably masks too much. Insect mitogenomes are about
75-80 % A+T, and much of the 13.7 % is likely AT-rich control region
rather than microsatellite, which would cost real recruitment. The
masker and its settings must be chosen on evidence (step 4), not by
default.

### Things the plan must account for

- **The `animal_mt` gate counts whole reads.** `coverage_gate.py`
  divides total recruited read bases by a nominal genome size
  (`BASIS_TOTAL_RECRUITED`, no sibling panel). Even a clean sample's
  gate estimate sits above its self-mapped depth, so the acceptance
  test is a ratio judged against clean samples, not a fixed factor.
- **Self-mapped depth has the same flaw.** A sample's own mt contigs
  carry their own (CA)n and AT-rich stretches, so mapping raw reads to
  them also attracts nuclear repeat reads. Use median depth, and count
  only alignments of at least about 1 kb or covering most of the read.
  Record the unfiltered figure as well so the filter's effect shows.
- **Masking also changes binning.** `BIN_TARGET` aligns contigs to the
  same `${assembly_target}.mmi`. A masked control region lowers a real
  mt contig's aligned fraction and could drop it below the binning
  floor.
- **`INT-ANIMAL-01` is pre-recruited** (SRR8306868), so masking may
  barely move it. It checks for regressions; it is not evidence that
  the fix works.

### Steps

1. **Baseline** (remote, current bundle `v2026.09_1`). For each
   evidence sample (round 2 insect, round 3, `INT-ANIMAL-01`):
   - recruited reads and bases, from `minimap2` + `recruit_filter.py`
     run directly in the recruit image (there is no RECRUIT-only
     pipeline entry point);
   - the gate estimate (recruited bases / nominal size), checked
     against the task 58 baseline's reported figure;
   - self-mapped depth to the selected contigs, filtered median and
     unfiltered mean;
   - the gate-to-depth ratio. The clean samples (round 2
     barcode05/06/13, round 3 barcode13) set the normal range.

   For barcode20, also list which panel records and regions its short
   (< 300 bp) hits land on, and whether they are short-unit repeats
   (period 1-6 bp) or broader low complexity. This decides which
   masker candidates go into step 4.
2. **Masking code.** Factor the "intervals in, N-substituted FASTA
   plus JSON report out" logic out of `scripts/refdata/mask_panel.py`
   into a shared function, keeping its length-preserving contract and
   its max-masked-fraction guard. Add a small entry point that takes
   repeat intervals (dustmasker `-outfmt interval`, or a pure-Python
   tandem-repeat scan if step 1 favours it). Reuse `bin/intervals.py`.
   Unit tests in `scripts/tests/`; the existing `test_mask_panel.py`
   must still pass.
3. **Build.** `build_recruit.sh`: add the `animal_mt` masking step
   (BLAST image, `--user`, if dustmasker), with its settings exposed
   as knobs. Emit `animal_mt.masking.json`, fold it into
   `build_manifest.py` next to `plant_mt.masking.json`, and rebuild
   the `.mmi`.
4. **Tune** (remote, outside the pipeline). Build candidate panels for
   2-4 settings (e.g. dust level 20/30/40, and/or tandem-repeat
   minimum length), then map barcode20 and the clean samples to each.
   This takes minutes per setting. For each, record:
   - barcode20's recruited bases, gate estimate and gate-to-depth
     ratio;
   - for the clean samples, retention: the fraction of reads that the
     current panel recruits and that map to the sample's own mt
     contigs, still recruited by the masked panel.

   Also record barcode08/14's recruited bases. The target (revised at
   the progress review) is that barcode20's short-hit bases (best
   alignment < 300 bp, 76 % of its recruited bases) are mostly gone,
   and its gate moves toward the 92x ceiling. Pick the loosest setting
   that does this while clean-sample retention stays at least 98 %.
   Confirm the 98 % figure with the user when the numbers are in.
5. **Bundle and evidence run.** New bundle `v2026.10` (VERSIONING.md).
   Upload it to Azure, copy it to the remote under
   `/mnt/data/wf5/refs/`, update the default bundle version, and check
   how the integration CI fetches refs. Run round 2 insect, round 3
   and `INT-ANIMAL-01` through the full pipeline.
6. **Plant exposure (measure only).** Run the chosen masker over
   `plant_mt.fa` and `plant_pt.fa` and report the masked fraction and
   repeat classes. If cheap, check the composition of recruited reads
   on `round_2_plant_t70`. Record the result in Outcomes and raise a
   todo if the plant panels are exposed.
7. **Spec / docs.** Update the spec section on the recruit panel and
   spec/07 §9; fold in the two todo.md recruitment items. Render a new
   workflow report and give the user its path.

### Acceptance

Revised at the progress review (user decision: masking only; the gate
basis is a follow-up todo).

- Post-mask, barcode20's gate-to-depth ratio falls within the
  post-mask range of the clean samples (barcode05/06/13). The ceiling
  is about 4.0, so barcode20's status is expected to stay `ok`. That
  is a known limit of the whole-read gate, not a failure of this task.
- barcode08 and barcode14 (95-98 % repeat hits) drop to `fail`; this
  is the intended status change.
- Bin selections and barcode recovery unchanged on barcode05/06/13.
- Clean-sample mt-read retention at least 98 % (or the agreed figure).
- `INT-ANIMAL-01` stays within its coverage and assembly bounds, or a
  bounds change is justified as biology (spec §5.1).
- Masking provenance in the manifest (CONSTITUTION principle 10).

### Out of scope

Runtime complexity filter in `recruit_filter.py`; non-zero
`min_aligned_bp` defaults; masking the plant panels.

Changing the `animal_mt` gate basis from whole-read bases to aligned
bases is a todo in tasks/todo.md (Recruitment), raised at the progress
review.

### Checklist

- [x] Baseline table: recruited, gate, self-mapped depth, ratio
- [x] barcode20 hit-location breakdown
- [ ] `INT-ANIMAL-01` baseline
- [x] Shared masking function, new entry point and tests; flake8 clean
- [x] `build_recruit.sh` / `build_manifest.py` changes
- [x] Tuning sweep: barcode20 ratio, clean-sample retention (level 30 recommended, awaiting confirmation)
- [ ] `v2026.10` built, distributed (Azure, remote, default version)
- [ ] Evidence run; acceptance met
- [ ] Plant exposure measured
- [ ] Spec and todo.md updated; report rendered

## Outcomes

- 2026-10-16: scoped and plan reviewed; dustmasker probe above. No
  code changed.
- 2026-10-16, step 1 baseline (remote, bundle `v2026.09_1`; raw reads
  re-cleaned with the pipeline's chopper + filtlong settings, recruited
  with minimap2 + `recruit_filter.py`; scripts in
  `/mnt/data/wf5/t68/` on the remote). Recruited bases reproduce the
  gate's `total_recruited_bases` exactly. Self-mapped depth = cleaned
  reads mapped to the selected contig(s), alignments of at least 500 bp
  and at least 80 % identity (PAF matches / block), median per base.
  The unfiltered figure for barcode20 is 840x, so an identity filter
  is essential.

  | sample | recruited reads | gate (x) | filtered depth (x) | gate / depth |
  |---|---|---|---|---|
  | R2 barcode05 | 576 | 143.3 | 59 | 2.4 |
  | R2 barcode06 | 104 | 26.1 | 4 (6 alignments) | 6.5 |
  | R2/R3 barcode13 | 5035 | 1340.1 | 1160 | 1.2 |
  | R3 barcode20 | 1491 | 390.8 | 23 | 17 |

  Other samples recruit too few reads to assemble a contig (gate
  `fail`, 0.1-9x), so they have no self-mapped depth. barcode06's depth
  rests on six alignments, so its ratio is noisy; the clean range is
  roughly 1-6.5.
- barcode20 hit locations: 76 % of its recruited bases (1151 of 1491
  reads, best alignment under 300 bp) are short hits. 950 of the 1151
  land on *Indri indri* KJ944196.1, all in its last 300 bp
  (16780-17078), a minisatellite of about 30 bp units built on (CA)n;
  147 on *Lepidodermella squamata* NC_026985.1. A strict period of at
  most 6 bp covers only 162 of the 1151, so a masker must catch
  composite and interrupted repeats: dustmasker or TRF, not a
  simple-period scan.
- barcode20's selected "mt" contig is a 5.8 kb fragment at 77 %
  identity and 0.37 aligned fraction, i.e. barcode20 has no real mt
  assembly; its true mt depth is about 21-23x (Flye 21x).
- 2026-10-16, progress review:
  - R2 barcode13 and R3 barcode13 are the same file (identical md5),
    as are the two barcode18s. The clean set is three independent
    samples (barcode05, barcode06, barcode13), not four.
  - The "clean" samples carry the same artefact: short hits are 15 %
    (barcode05), 20 % (barcode06) and 26 % (barcode13) of recruited
    bases, and barcode13's land on the same *Indri* tail. Their
    baseline ratios are inflated too, so acceptance must compare
    post-mask ratios, not barcode20 against the pre-mask range.
  - R2 barcode08 (gate 35.6x, `ok`) and barcode14 (26.2x,
    `low_coverage`) are 98 % and 95 % short repeat hits, mostly the
    *Indri* tail. barcode08 went to assembly, produced one contig,
    selected nothing and recovered no barcode. Masking should drop
    both to `fail`. They belong in the evidence set as the intended
    status changes.
  - Ceiling, if masking removed every short-hit read: barcode20 gate
    92x / depth 23x = ratio 4.0, still `ok`, while its true depth
    would be `low_coverage`. Clean samples: barcode05 2.1, barcode06
    5.2, barcode13 0.85. Masking alone cannot correct barcode20's
    status, because the gate counts whole read bases.
  - `INT-ANIMAL-01` has no baseline yet (step 1 item outstanding).
- Benchmark samples: repeat content of recruited bases. Repeat content
  is the share of recruited bases from reads whose best panel
  alignment is under 300 bp. "After" is measured with the dustmasker
  level 30 panel (step 4 sweep, 2026-10-19), built with
  `mask_repeats.py`; the surviving short hits on clean samples are
  residual conserved-gene fragments, not microsatellites.

  | sample | class | recruited bases (before) | gate (before) | repeat content (before) | repeat content (after) |
  |---|---|---|---|---|---|
  | R2 barcode05 | clean | 2,436,939 | 143.3x `ok` | 15.5 % | 11.9 % (d30) |
  | R2 barcode06 | clean | 445,099 | 26.1x `low_coverage` | 20.5 % | 14.6 % (d30) |
  | R2/R3 barcode13 | clean | 22,782,397 | 1340.1x `ok` | 26.5 % | 0.2 % (d30) |
  | R3 barcode20 | affected | 6,644,227 | 390.8x `ok` | 76.3 % | 9.9 % (d30) |
  | R2 barcode08 | affected | 604,736 | 35.6x `ok` | 97.9 % | 0 % (4,457 bp recruited) |
  | R2 barcode14 | affected | 445,517 | 26.2x `low_coverage` | 95.4 % | 23.3 % (12,023 bp recruited) |
  | R2 barcode02 | fail | 2,161 | 0.1x `fail` | 100 % | 0 bp recruited |
  | R2 barcode03 | fail | 85,864 | 5.0x `fail` | 96.2 % | 0 % (3,229 bp recruited) |
  | R2 barcode04 | fail | 118,661 | 7.0x `fail` | 26.8 % | 33.7 % (57,936 bp recruited) |
  | R2/R3 barcode18 | fail | 109,551 | 6.4x `fail` | 93.6 % | 23.6 % (9,131 bp recruited) |
  | R3 barcode16 | fail | 152,788 | 9.0x `fail` | 100 % | 0 bp recruited |
  | `INT-ANIMAL-01` | regression | pending | 327.5x (task 58) | pending | pending (not yet run) |

- 2026-10-19: steps 2 and 3 done. `mask_panel.py` now delegates to a
  shared `mask_fasta`; new `scripts/refdata/mask_repeats.py` (dustmasker
  `-outfmt interval`, inclusive ends converted to half-open),
  `test_mask_repeats.py`, and `build_recruit.sh` / `build_manifest.py`
  emit and record `animal_mt.masking.json`. Level 30 is the provisional
  default pending confirmation. A full bundle build has not been run.
- 2026-10-19: step 4 sweep (remote, `/mnt/data/wf5/t68/sweep/`,
  summary via `sweep_summary.py`). Dust window 64, linker 1.

  | level | bases masked | records dropped by 60 % guard |
  |---|---|---|
  | 20 | 12.7 % | 32 |
  | 30 | 4.7 % | 0 |
  | 40 | 1.7 % | 0 |

  | sample | panel | recruited bp | gate | repeat share | retention | gate/depth |
  |---|---|---|---|---|---|---|
  | barcode20 | base | 6,644,227 | 390.8x | 76.3 % | - | 17.0 |
  | barcode20 | d40 | 1,608,166 | 94.6x | 8.2 % | 100 % | 4.1 |
  | barcode20 | d30 | 1,298,167 | 76.4x | 9.9 % | 100 % | 3.3 |
  | barcode20 | d20 | 508,452 | 29.9x | 20.2 % | 100 % | 1.3 |
  | barcode05 | base / d30 | 2,436,939 / 2,098,240 | 143.3x / 123.4x | 15.5 / 11.9 % | 100 % | 2.4 / 2.1 |
  | barcode06 | base / d30 | 445,099 / 221,790 | 26.2x / 13.0x | 20.5 / 14.6 % | 100 % | 6.5 / 3.3 |
  | barcode13 | base / d30 | 22,782,397 / 16,320,764 | 1340x / 960x | 26.5 / 0.2 % | 99.9 % | 1.2 / 0.8 |

  Retention = reads that both the baseline panel recruited and that
  align (>= 500 bp, >= 0.8 identity) to the sample's own contigs, still
  recruited. Level 30 meets the 98 % bar with 100 % (99.9 % barcode13),
  and brings barcode20 into the clean-sample ratio range (3.3 against
  2.1-3.3 post-mask) with status still `ok`. Level 20 is tighter but
  drops 32 real genomes from the panel, so it is rejected. barcode08
  and barcode14 fall to 0.3x and 0.7x (`fail`) as intended;
  barcode06 drops to 13.0x (still `low_coverage` territory, was 26.2x).
  Recommendation: level 30. Awaiting user confirmation of level and the
  98 % retention figure.
