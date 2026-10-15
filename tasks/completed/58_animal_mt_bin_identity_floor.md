# Task 58 — Lower the `animal_mt` binning identity floor from 80% to 70%

## Overview

After Flye assembles a sample, `BIN_TARGET` (custom component C3,
`bin/bin_target.py`) decides which contigs are the target mitogenome.
A contig qualifies only if its alignment against the `animal_mt`
reference panel (`refs/<ver>/recruit/animal_mt.mmi`) clears two floors:

- **aligned fraction**: at least 30% of the contig aligns to the panel;
- **identity**: those aligned bases are at least 80% identical.

The 80% identity floor was calibrated on a single integration fixture
(`INT-ANIMAL-01`, an aphid, which is near-conspecific to a panel
reference). Real client insects now show it is set too high. Every
genuine insect mitogenome selected so far has only just cleared it, and
a divergent taxon fails it outright:

| Run | Sample | Selected contig | Identity to panel | Aligned frac |
|---|---|---|---|---|
| round 2 | barcode05 | contig_1, 15.9 kb, circular | 82.30% | 0.962 |
| round 2 | barcode06 | contig_1, 5.9 kb | 81.37% | 0.993 |
| round 2/3 | barcode13 (*Aedes*) | contig_1, 14.6 kb | 80.64% | 0.938 |
| round 3 | barcode20 (Cicadellidae) | none → `no_assembly` | contig_65: 76.96%; contig_87: 75.26% | 0.368; 0.644 |

Insect mitochondrial DNA evolves fast. A taxon from a family or genus
with no close relative in the panel will typically sit at 70–80%
nucleotide identity to its nearest reference, even if the assembly is
perfect. A BLAST of a round 3 barcode20 COX1 read against all of RefSeq
mitochondrion (`refs/v2026.09_1/validate/refseq_mt_metazoa`) gives a best
hit of only 79.4% (*Atkinsoniella*), and every top hit is a leafhopper.
That clade simply isn't in RefSeq at that depth. For biosecurity the
under-referenced taxon is exactly the case that matters most, and this
floor rejects it preferentially. tasks/todo.md already raises the same
argument for the barcode identity floor (C5).

The intent is to lower `bin_target_thresholds.animal_mt.min_identity` to
**70%**. This matches `plant_mt` and leaves the aligned-fraction floor
(0.30) as the main thing that separates target from off-target contigs.
That floor carries the separation on `INT-ANIMAL-01` (0.878/0.985 for
true contigs vs ≤ 0.069 for off-target), and the same holds in round 3:
the barcode20 contigs between 70% and 80% identity with aligned fraction
≥ 0.30 are exactly the two mito-like fragments.

**This does not fix recruitment.** The round 3 investigation found that
barcode20's 391× coverage-gate figure is mostly nuclear reads recruited
on < 300 bp microsatellite hits (for example the (CA)n tail of the
*Indri indri* reference). Self-mapping puts true mitochondrial depth in
the low single digits to ~20×. That is a separate problem (see Out of
scope). As a consequence, barcode20 is expected to move from
`no_assembly` to a selected short fragment, which is a *partial*
recovery, not a confident one. That is the correct outcome under
CONSTITUTION principle 7 ("a partial recovery is not a negative"),
provided it is reported honestly (§3.3).

## 1. Constraints

- **CONSTITUTION principle 2 / hard constraint:** the change is keyed by
  `assembly_target` and touches `animal_mt` only. `plant_pt` (75%) and
  `plant_mt` (70%) are unchanged.
- **Rule 18 / spec §3.7.6:** the value stays in versioned config
  (`params.bin_target_thresholds` in `nextflow.config`) and is echoed
  into `bin_metadata.json` as `thresholds_applied`. No constant moves
  into `bin/bin_target.py`.
- **Spec §5.1:** the integration fixtures are correctness tests, not
  benchmarks. 70% is still a *prototype default* under spec §9 item 10,
  justified here by real-data evidence, not by a fixture sweep. Say so
  in every place the value is documented.
- **Principle 7:** a newly admitted low-identity contig must stay
  distinguishable from a confident `ok` at every reporting surface.

## 2. Config change

- `nextflow.config`: `bin_target_thresholds.animal_mt.min_identity`
  80.0 → 70.0.
- Extend the comment block above `bin_target_thresholds` with one or two
  lines saying the `animal_mt` floor was lowered on round 2/3 real-data
  evidence (task 58). Real insect mitogenomes sit at 80.6–82.3% to the
  panel, and divergent clades at 75–79% even against all of RefSeq.
- Check that no other profile or config (`conf/integration.config`,
  `conf/stub.config`, `tests/manual/*.config`) overrides the value.

## 3. Verification

### 3.1 Unit tests (`scripts/tests/test_bin_target.py`)

No logic changes in `bin/bin_target.py`, and the tests pass explicit
thresholds, so no existing test should change. Confirm (and add if
absent) a boundary test on the identity floor:

- contig identity == `min_identity` → admitted as a target candidate;
- identity just below → `off_target`;
- a contig at ~75% identity with aligned fraction ≥ 0.30 is admitted
  under a 70 floor and rejected under an 80 floor (pseudocode:
  `classify(row(id=75, af=0.6), thresholds(min_identity=70))` →
  `target_candidate`; same row with 80 → `off_target`).

Keep 100% branch coverage of C3 (rule 14). Run via `scripts/pytest.sh`.

### 3.2 Integration fixture (`-profile integration`)

`INT-ANIMAL-01` must still select the same single circular contig, and
`tests/integration/expected/animal_mt/bin_bounds.json` must pass
unchanged (14–20 kb, 1 contig, circular). Before running, read the
fixture's `bin_metadata.json` from a previous run. If any off-target
contig there has aligned fraction ≥ 0.30 and identity 70–80%, that
contig will now become a candidate. Record it, and do not adjust the
bounds to accommodate it without discussion. No fixture data changes
are expected.

### 3.3 Real-data re-run (remote, via the `remote-run` skill)

Re-run round 3 (`tests/manual/samples_round_3.csv`, data already at
`/mnt/data/wf5/data/round_3`) and round 2 animal_mt
(`tests/manual/samples_round_2_insect.csv`) clean, since both remote
work dirs were deleted. Expected:

- round 2 barcode05/06/13 and round 3 barcode13: **unchanged**
  selection, status and barcodes.
- round 3 barcode20: `contig_65` (5.8 kb, 21×, 76.96%, af 0.368) is
  admitted and selected as the highest-coverage candidate, with
  `contig_87` logged as a secondary. Record the resulting
  `sample_status` and per-locus barcode outcome. `no_barcode` or a
  partial/low-locus `ok` are both plausible. If a barcode is emitted,
  check its BLAST_VALIDATE / `validation.tsv` evidence that it is
  Cicadellidae mitochondrial (top hits leafhoppers), not a NUMT or
  chimera.
- any other sample whose status changes: explain it.

Contig_65's coverage (21×) is indistinguishable from the nuclear
contigs in the same assembly (≤ 28×), so coverage cannot rule out a
NUMT (spec §3.7.3, §9 item 10a). Note this in the outcomes rather than
treating the selection as settled.

### 3.4 Report check

Render the round 3 barcode20 `report.html` and confirm that the
low identity (Reference axis) and short, non-circular target are
visible, so the sample does not read as a full-depth confident
recovery. If the report gives no such signal, record it as a follow-up
rather than expanding this task.

## 4. Spec updates

- `spec/03-organelles.md` §3.7.6: change the `animal_mt` row to 70%,
  and rewrite its Rationale cell to cite the real-data identities above
  and the aligned-fraction floor as the primary separator. Keep the
  "prototype defaults" paragraph and add that `animal_mt` was revised on
  real client data (task 58_animal_mt_bin_identity_floor.md).
- `spec/07-open-questions.md` §9 item 10: update the current-value
  column to 70/75/70%, and note that round 2/3 client insects are the
  first real `animal_mt` evidence for the sweep.
- `spec/02-stages.md`: check that the stage 9 row and the C3 row do
  not quote 80%. Update them if they do.

## 5. Out of scope (record in tasks/todo.md if not already there)

- **Recruitment inflation by microsatellites.** `recruit_thresholds`
  ship at 0/0, so any read with any panel alignment is recruited. On
  round 3 barcode20, 78% of recruited reads came in on < 300 bp hits to
  simple-sequence repeats in panel control regions (*Indri indri* (CA)n,
  *Lepidodermella squamata* (GA)n). This inflated the gate to 391× and
  sent 94 nuclear contigs through Flye. A `min_aligned_bp` floor alone
  leaves residual inflation (1000 bp → ~65× estimated, still
  overstated). Low-complexity masking of the `animal_mt` panel, similar
  to the `plant_mt` NUPT masking, looks like the more direct fix. Add
  this to the existing "Recruitment" item in tasks/todo.md with these
  figures. It should become its own task, using round 2/3 as evidence.
- `min_aligned_frac` for `animal_mt` (0.30) is unchanged. Contig_65
  clears it only narrowly (0.368). Leave any revision to the §9 item 10
  sweep.
- The C5 barcode identity floor (tasks/todo.md, "Barcode extraction")
  is a separate threshold with its own open item.

## 6. Outcomes

### Done (2026-10-13)

- `nextflow.config`: `animal_mt.min_identity` 80.0 → 70.0, comment added.
  No override in `conf/*.config` or `tests/manual/*.config`.
- Unit tests: added `test_identity_floor_boundary` (identity == floor
  admitted, +0.01 rejected) and `test_divergent_contig_admitted_under_lower_floor`
  (~75% identity contig: 70 floor → `target_candidate`, 80 → `off_target`).
  `scripts/pytest.sh`: 56 passed, `bin/bin_target.py` 100% branch coverage.
  flake8 clean. No logic change in `bin/bin_target.py`.
- Fixture pre-check (§3.2): from the previous `INT-ANIMAL-01`
  `bin_metadata.json` (`output/integration_20261012_190146_140425`), the
  off-target contigs with identity 70–80% (contig_1, 3, 2) have aligned
  fraction 0.079 / 0.060 / 0.019, all ≪ 0.30, so none becomes a candidate.
  Bounds unchanged. The full `-profile integration` run was not repeated.
- Spec: `03-organelles.md` §3.7.6 row + prototype-defaults paragraph,
  `07-open-questions.md` §9 item 10 updated. `02-stages.md` does not
  quote 80%; no change.

### Done (2026-10-15): remote re-runs and live integration run

- §3.2 live `-profile integration` (animal_mt): INT-ANIMAL-01 selects the
  same single circular contig (contig_10, aligned fraction 0.92);
  `thresholds_applied.min_identity` = 70.0; all off-target contigs have
  aligned fraction <= 0.113. No animal assertion failed. Bounds unchanged.
- §3.3 round 2 insect (`tests/manual/output/round_2_insect_t58/`):
  statuses and barcodes unchanged. barcode05 `ok` (contig_1, 82.30%),
  barcode06 `low_coverage` (81.37%; Flye renumbered the selected contig
  1 -> 2), barcode13 `ok` (80.64% -> 80.32%, assembly variation). 02/03/04/18
  `fail` and 08 `no_assembly`, as before. barcode14 METAFLYE assembled no
  contigs (ignored task failure, no metadata.json), as in the earlier run.
- §3.3 round 3 (`tests/manual/output/round_3_t58/`): barcode13 `ok`
  (same COX1 identity 0.9336), barcode16/18 `fail` (coverage gate),
  **barcode20 `no_assembly` -> `no_barcode`**. contig_65 (5,796 bp, 21x,
  76.96%, af 0.368, non-circular) is selected; the secondary is contig_80
  (3,705 bp, 16x, 75.26%, af 0.644; contig_87 in the earlier run, Flye
  renumbered). Annotation finds only ND4/ND5 (protein-coding completeness
  0.15), 0/6 barcode loci. BLAST validation top hit is *Atkinsoniella
  warpa* (leafhopper) at 80.40%, 37% coverage, so the fragment is
  Cicadellidae mitochondrial, not an obvious chimera. NUMT is **not ruled
  out**: 21x is indistinguishable from nuclear contigs (<= 28x).
- §3.4 report: barcode20 `report.html` shows the `no_barcode` outcome,
  target length, non-circularised, annotation 0/6 and the BLAST
  evidence. It does not show the 76.96% panel identity of the selected
  contig. **Follow-up** (not done here): surface bin identity / aligned
  fraction on the report.
- Remote work dirs removed. Task 68 (microsatellite recruitment) can use
  these runs as its baseline.
