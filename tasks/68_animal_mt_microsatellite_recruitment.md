# Task 68 — `animal_mt` recruitment inflated by microsatellites (stub)

**Phase:** P5, real-data tuning (spec/07-open-questions.md §9).

**Status:** stub. Scope it fully before starting. Follows on from task
58 (tasks/completed/58_animal_mt_bin_identity_floor.md, §5 Out of
scope), and depends on that task's unfinished §3.3 round 2/3 re-run for
its baseline.

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
