# Task 75 — Collapse circular `animal_mt` contigs that hold the genome twice

**Phase:** P5 — production deployment (spec/06-phases.md).

**Status:** complete (2026-10-22). No preconditions (task 76, rotation, is complete and
this task builds on it).

**Related spec:** spec/01-pipeline-flow.md §1 (stage 9a),
spec/02-stages.md stage table row 9a and §2.2 C13,
spec/03-organelles.md §3.1 (typical sizes), §3.7.4 (circularity comes
from Flye) and §3.8 (animal mitogenome rotation), spec/06a-reports.md
§6a.2. CONSTITUTION principles 7 and 8, rules 14, 16, 18 and 19.

## Overview

metaFlye sometimes closes a circular organelle contig that holds the
genome **twice in a row** (a tandem duplicate, or concatemer) instead
of once. The contig is flagged circular, aligns well to the organelle
panel and carries every gene, so BIN_TARGET selects it as the target.
But it is about twice as long as the real genome, and everything
downstream that treats it as one genome is wrong: the deliverable
FASTA, the length, the gene order, the organelle map, and the barcodes.

This has happened once, in task 68 (2026-10-20). On the `INT-ANIMAL-01`
fixture (*Acyrthosiphon pisum*, mitogenome about 17 kb), after
recruitment changed the read pool that Flye subsampled, the target was
one **33,905 bp** circular contig (Flye `circ.=Y`, 58x, multiplicity
9). The sample still came out `ok`, "clean recovery", and its
`barcodes.fasta` held **COX2 twice**, once from each copy
(`COX2_contig_1_8624_9294` and `COX2_contig_1_25577_26244`). Only the
integration test's target-length bound caught it. On a client sample
nothing would have flagged it.

The cause is the assembler on a particular read set, not anything the
pipeline did before it. Flye has no seed option on this path (task
76), so a doubled contig can turn up on any sample. The fix belongs
after assembly.

This task makes C13 (`bin/rotate_target.py`, stage 9a) detect a
circular `animal_mt` target that repeats itself end to end, and
**collapse it to a single copy** (decided with the user, 2026-10-22).
The collapsed copy starts at the rotation anchor (COX1), so the
deliverable is the same as for a normal assembly. The collapse is
recorded in `rotation.json` and `metadata.json` and shown in the
report, because the pipeline has changed the assembler's sequence and
the auditor must be able to see that (rule 16, rule 18).

### Evidence gathered while scoping (2026-10-22)

The doubled contig from task 68 survived in an old scratchpad. It is
now kept as a unit-test fixture, with a normal contig as a control,
in `scripts/tests/fixtures/tandem_duplicate/`:

| File | What it is |
|---|---|
| `INT-ANIMAL-01.doubled.fasta` | The 33,905 bp BIN_TARGET target (`contig_1`), pre-rotation |
| `INT-ANIMAL-01.doubled.bin_metadata.json` | Its `bin_metadata.json` (`contig_1`: circular, `flye_circ`, 58x) |
| `INT-ANIMAL-01.single.fasta` | A normal 16,960 bp target (`contig_10`, rotated) from task 76's remote run 1 |

A circular 31-mer test separates the two cleanly. For each 31-mer seen
exactly twice, take the circular distance between the two copies:

| Contig | Length | Dominant offset | Positions whose 31-mer recurs at that offset (±50 bp) |
|---|---|---|---|
| Doubled | 33,905 bp | 16,951 bp | 0.77 |
| Single (`INT-ANIMAL-01`) | 16,960 bp | none (top offset 124 bp, 31 k-mers) | 0.002 |
| `INT-PLANT-01-pt` | 155,276 bp | none (9 k-mers) | 0.000 |

In the doubled contig, the offsets cluster at 16,947–16,952 bp,
because the two copies differ by a few indels. 2 × 16,951 = 33,902,
3 bp short of the contig length, so the period tiles the circle. The
contig's leading 31-mer occurs again at exactly position 16,951
(0-based), which gives an exact cut point.

**Frequency.** None of the 18 `animal_mt` targets in the local outputs
(`tests/manual/output/`, `tests/integration/output/`) is longer than
17 kb. Remote and Azure outputs are surveyed in §3.0.

## 1. Constraints

- **Collapse, don't just flag.** A flag that left both copies in place
  would still ship duplicate barcodes and a wrong map. The collapsed
  copy is the deliverable.
- **Fail closed (rule 18).** Collapse only when the evidence is
  unambiguous (§2.2). Record a contig that looks repetitive but does
  not tile cleanly as `unresolved`, leave its sequence unchanged, and
  warn. Never guess a cut.
- **No new failure mode (principle 8).** Detection and collapse never
  fail a sample. C13 still exits 0 except on unreadable input.
- **`animal_mt` only.** Detection runs only where it is configured
  (§3.1). `plant_pt` is rebuilt from the graph by C4 (spec §3.6).
  `plant_mt` has genuine large repeats and no single circle, and
  collapsing it would destroy real biology.
- **Circular contigs only.** The collapse relies on the copies joining
  end to start around the circle. Linear contigs are `not_tested`.
- **Config, not code.** The k-mer size and thresholds are pipeline
  params in `nextflow.config`, validated in `validateParams()` and
  passed to C13 as explicit CLI args. They are not constants in
  `rotate_target.py`.
- **Rule 14.** Stdlib only (the k-mer index is a dict over a contig of
  at most a few tens of kb). 100% branch coverage on the new code.
- **Rule 19.** Integration assertions check status values, not exact
  lengths or offsets.

## 2. Definitions

### 2.1 Period and support

For a circular contig `seq` of length `L`, with `k = params.tandem_collapse.<target>.kmer`
(default 31):

```
index every k-mer of seq read around the circle (seq + seq[:k-1])
for each k-mer occurring 2..MAX_COPIES times:
    for each pair of adjacent occurrences:
        d = circular distance, min(delta, L - delta)
        histogram[d] += 1
P = histogram's highest-count offset            (the candidate period)
support = fraction of the L positions whose k-mer also occurs at
          circular offset P ± tolerance_bp
copies  = round(L / P)
```

`MAX_COPIES` can be a small module constant (for example 4). It is a
bound on the search, not a biological threshold.

### 2.2 Decision

| Condition | `tandem.status` | Sequence |
|---|---|---|
| Not configured for this target, not circular, or empty | `not_tested` | unchanged |
| `support < min_support`, or `P < min_period_bp`, or no repeated k-mers | `single` | unchanged |
| `support ≥ min_support`, `P ≥ min_period_bp`, `copies ≥ 2`, and `abs(L − copies × P) ≤ copies × tolerance_bp` | `collapsed` | one copy (§2.3) |
| `support ≥ min_support` and `P ≥ min_period_bp`, but the period does not tile the circle | `unresolved` | unchanged, report warning |

The comparison is `support ≥ min_support` (the boundary counts as
repetitive). Defaults: `kmer` 31, `min_support` 0.5,
`tolerance_bp` 50, `min_period_bp` 1,000. The measured gap is
0.002 versus 0.77, so 0.5 sits well inside it. `min_period_bp` keeps a
long microsatellite run from counting as a genome period. Mark all of
these **provisional**, calibrated on one doubled and two single
contigs.

A genuine partial duplication (a duplicated control region or gene
block) gives a short `P` with low support, or a `P` that does not tile
the contig. It is `single` or `unresolved`, never `collapsed`.

### 2.3 Which copy is kept, and where the cut falls

Detection does not depend on rotation, so it runs on the input
sequence. The cut is made **after** rotation, so that the kept copy
starts at the anchor gene:

```
rotated = rotate_contig(seq, anchor)        # task 76, unchanged; or seq if no_anchor
lead    = rotated[:k]
cut     = first occurrence of lead in rotated within [P - tol, P + tol]
          (fallback: cut = P, and record cut_method "period")
output  = rotated[:cut]
```

The kept copy is the first one in the rotated frame. It contains the
anchor at position 1, and because the copies join end to start, the
kept copy is itself a closed circle. With a `no_anchor` contig the cut
is made in Flye's frame, and the existing `no_anchor` warning still
applies.

Picking the first copy is a convention, not a quality judgement. Both
copies are consensus from the same reads. Choosing between them by
annotation quality is out of scope (§7).

## 3. Scope

### 3.0 Survey of past runs (record in Outcomes, not blocking)

Run the §2.1 test (C13's new function, called directly) over every
`animal_mt` target available:

- local `tests/manual/output/` and `tests/integration/output/` (done
  during scoping: 18 targets, none doubled; re-check with the final
  function);
- remote outputs on daff-admin under `/mnt/data/wf5/output_*` (see the
  `remote-run` skill for access);
- Azure outputs already downloaded under `tests/manual/output/az_*`.

Record per target: sample, length, `P`, support, decision. If any real
sample is `collapsed` or `unresolved`, say so to the user before going
further. It is either a second real case, or a false positive that
moves the thresholds.

### 3.1 Config (`nextflow.config`, `main.nf`)

Add, next to `rotation_anchors`:

```
tandem_collapse = [
    animal_mt: [ kmer: 31, min_support: 0.5, tolerance_bp: 50,
                 min_period_bp: 1000 ],
    plant_pt:  null,
    plant_mt:  null,
]
```

The comment cites spec §3.8, this task, and marks the values
provisional. `validateParams()` checks that every assembly target has an
entry (`null` allowed), `kmer` is an integer in [15, 63], `min_support`
is in (0, 1], and `tolerance_bp` and `min_period_bp` are positive
integers. Each failure names the param.

### 3.2 C13 `bin/rotate_target.py`

- New pure functions: one that returns `(P, support, copies)` for a
  sequence (§2.1), one that applies the §2.2 decision, and one that
  finds the cut (§2.3).
- `process_contig()` runs detection on every circular contig when
  tandem config is given. The order is: detect, rotate (unchanged),
  then cut. Rotation status values and logic are unchanged. A doubled
  contig is still `rotated` (or `no_anchor`), and the collapse is
  reported separately.
- New CLI args: `--tandem-kmer`, `--tandem-min-support`,
  `--tandem-tolerance-bp`, `--tandem-min-period-bp`. Detection is
  disabled when they are absent, so `plant_pt` and `plant_mt` get
  `not_tested`.
- Each contig record in `rotation.json` gains:

  ```
  "length_bp": 16951,            // output length (was input length)
  "input_length_bp": 33905,
  "tandem": {
    "status": "collapsed",       // not_tested | single | collapsed | unresolved
    "period_bp": 16951,
    "copies": 2,
    "support": 0.768,
    "cut_bp": 16951,             // null unless collapsed
    "cut_method": "leading_kmer",// leading_kmer | period | null
    "params": { "kmer": 31, "min_support": 0.5,
                "tolerance_bp": 50, "min_period_bp": 1000 }
  }
  ```

- Bump the schema id to `wf5/rotation/v2`. C13 is no longer a pure
  permutation, and `length_bp` changes meaning. Update the module
  docstring to match.
- The FASTA header line is still written back as read. The length
  change is recorded in `rotation.json`, not in the header.

### 3.3 Wiring

- `modules/local/rotate_target.nf`: pass the four args when the
  target's `tandem_collapse` entry is not null. Build the arg string in
  `main.nf` per sample beside the existing `rotation_anchors` lookup,
  the same way the anchor list is passed (a `val` in the input tuple).
- Nothing downstream changes. `BLAST_VALIDATE`, `MINIPROT_CDS`,
  annotation, extraction, the map and COLLATE already consume
  `ROTATE_TARGET.out.binned`, so they all see the collapsed copy.
  Task 69's target length (`fasta_total_length`) reads the same FASTA
  and needs no change.
- `bin_metadata.json` is not rewritten. It records what BIN_TARGET saw
  (the 33,905 bp contig), and `rotation.json` records the change.

### 3.4 Metadata and report

- `assets/sample_metadata.schema.json`: update the `rotation`
  description to `wf5/rotation/v2` and mention the `tandem` record.
  The `contigs` array is untyped, so no structural change is needed
  unless C6 is changed to validate it.
- `bin/report/report.py` `rotation_lines()`: add one line per contig
  whose tandem status is `collapsed` or `unresolved`:
  - `collapsed`, severity `info`: "Flye assembled this contig as 2
    tandem copies of the genome (33,905 bp). One copy was kept
    (16,951 bp), starting at COX1."
  - `unresolved`, severity `warning`: "This contig repeats itself
    (period 16,951 bp, 77 % of 31-mers) but the repeat does not tile
    the circle, so it was not collapsed. Length, gene order and
    barcodes may count some genes twice."
  - `single` and `not_tested`: no line.
- The `unresolved` warning goes into the warnings mirror (Assembly
  tab), like the `no_anchor` warning. Check that `collapsed` shows on
  the Assembly tab and appears nowhere as a failure.

### 3.5 Flye

No Flye parameter change in this task. Flye's circular-contig
duplication is not configurable on this path, and the subsample seed
is already pinned (task 76). Record this in spec §3.8, so the next
person does not go looking for a seed.

## 4. Testing

### 4.1 Unit tests (`scripts/tests/test_rotate_target.py`)

Run with `scripts/pytest.sh`. The new code needs 100% branch coverage.

Real fixtures (`scripts/tests/fixtures/tandem_duplicate/`):

- Doubled: `collapsed`, `copies` 2, `period_bp` within 16,951 ± 5,
  `support` > 0.7, `cut_method` `leading_kmer`. The output is
  `input[:cut]` when there is no anchor. With a synthetic COX1 hit
  placed in copy 2 (a hand-written anchor GFF line in the test), the
  output is rotated so that hit is at position 1, and its length is the
  cut.
- Single: `single`, `support` < 0.01, sequence unchanged.

Synthetic sequences (seeded `random`, built in the test, not files):

- A random 15 kb genome ×2 with about 0.5 % substitutions gives
  `collapsed`. The same genome ×3 gives `collapsed` with `copies` 3,
  keeping one copy.
- A segmental duplication (5 kb of a 15 kb genome repeated once) gives
  `single`.
- Two copies plus a 3 kb unique insert, so the period does not tile,
  gives `unresolved` with the sequence unchanged.
- A long `(AT)n` run gives `single`, because of `min_period_bp`.
- `support` exactly equal to `min_support` counts as repetitive.
- The leading k-mer is mutated in copy 2, so the cut falls back to
  `P` with `cut_method` `period`.
- Linear contig, empty target, and no tandem args all give
  `not_tested`, with the sequence unchanged.
- A `no_anchor` doubled contig is collapsed in Flye's frame, and the
  rotation status stays `no_anchor`.
- CLI: the four args parse, and leaving them out disables detection.
  `rotation.json` carries `wf5/rotation/v2`, `input_length_bp` and the
  `tandem` block.

Report (`scripts/tests/test_report.py`): one test per tandem status
for `rotation_lines()`, covering text and severity. `unresolved`
appears in the warnings mirror. Update any existing test that asserts
`wf5/rotation/v1` or `length_bp` equal to the input length.

Collate (`scripts/tests/test_collate.py`): a `rotation.json` with a
`tandem` block passes through into `metadata.json`, and the result
validates against the schema.

### 4.2 Stub run

The stub `ROTATE_TARGET` is unchanged (it touches outputs). Confirm the
stub run passes with the new `val` in the input tuple.

### 4.3 Integration (`tests/integration/assertions.sh`)

No fixture reproduces the doubled contig on demand. Flye is not
deterministic, and building one by assembly would be fragile (rule
19). The real doubled contig is covered by the §4.1 unit tests.
Record that gap in the assertions' header comment.

- `INT-ANIMAL-01`: every target contig's `tandem.status` is `single`
  or `collapsed`. Fail on `unresolved` or `not_tested`. If it is
  `collapsed`, print the input and output lengths as `OK:` lines, so a
  nightly occurrence is visible in the log. The existing
  `bin_bounds.json` target-length bound (14,000–20,000 bp) now applies
  to the output length. Check it reads the published
  `organelle_assembly.fasta`, not `bin_metadata.json`.
- `INT-PLANT-01-pt` and `INT-PLANT-01-mt`: every contig's
  `tandem.status` is `not_tested`.

### 4.4 Real-data acceptance

1. Run `ROTATE_ANCHORS`' miniprot by hand (pinned container, `docker
   run --user`) over the COX1/COX2/COX3/CYTB panels from
   `refs/v2026.10_1` against `INT-ANIMAL-01.doubled.fasta`, under
   tables 2 and 5. Then run `rotate_target.py` with the tandem args.
   Expected: `rotated`, COX1, `collapsed`, output about 16,951 bp.
2. Run miniprot over the full `animal_mt` panel on the collapsed
   output. Each anchor gene must have its best hit once (not twice),
   with COX1 starting at 1 on `+`.
3. Align the collapsed copy to `INT-ANIMAL-01.single.fasta` (minimap2,
   `asm5`). Expect a single near-full-length alignment at about 99 %
   identity or better. They are different assemblies of the same
   reads, so the match will not be exact.
4. Run the full `-profile integration` once. All assertions must pass.
   Render the run report and give the user its path, per CLAUDE.md.

Record the commands and numbers in Outcomes.

## 5. Spec and doc updates

- spec/03-organelles.md §3.8: add a "Tandem duplicates" paragraph. Cover
  what Flye does, the task 68 evidence, the §2.2 decision table, the
  first-copy convention, why it is `animal_mt` only, and that Flye has
  no setting to prevent it.
- spec/03-organelles.md §3.1: in the `animal_mt` "Assembly pitfall"
  cell, add "can close as a tandem duplicate; collapsed by C13 (§3.8)".
- spec/02-stages.md row 9a and C13 row: replace "Rotation is a pure
  permutation — same length and bases" with the collapse behaviour.
  Add `tandem` statuses and `wf5/rotation/v2`, and add this task to
  the C13 task list.
- spec/01-pipeline-flow.md §1: one more line on stage 9a ("collapses a
  circle that holds the genome twice to one copy").
- spec/06a-reports.md §6a.2: the two new Assembly-tab lines.
- `scripts/tests/README.md`: one line on the `tandem_duplicate/`
  fixtures and where they came from.
- tasks/todo.md: update task 76's follow-up (c) to point at the
  completed task, and add §7's items.
- CONSTITUTION: no amendment needed. Principle 7 and rule 18 already
  require that an output which is wrong in a way the pipeline can
  detect must not ship looking normal.

## 6. Design decisions

### 6.1 Collapse rather than flag

With flag-only, the sequence would ship with two copies of every gene.
On the real case, that put COX2 twice in `barcodes.fasta`. It would
also give a double-length organelle map and a target length that hides
a partial recovery (task 69). The collapsed copy is what a correct
assembly would have produced. Changing the assembler's output is the
same kind of step as rotation and C4's plastid canonicalisation, and it
is recorded the same way.

### 6.2 K-mer self-similarity rather than self-alignment or gene counts

- **Self-alignment** (minimap2 / mappy) would need a new dependency in
  the `wf5/scripts` image (rule 19). The k-mer index is stdlib, runs in
  milliseconds on 40 kb, and the measured separation (0.77 against
  0.002) leaves plenty of margin.
- **Counting anchor genes** (COX1 found twice) depends on annotation
  quality and on how divergent the panel is. A doubled contig whose
  second COX1 copy is poorly aligned would be missed, and a genuine
  gene duplication would be collapsed by mistake. The k-mer test
  checks the property itself, a whole-sequence period.

### 6.3 In C13, not BIN_TARGET

C13 already handles circular `animal_mt` targets. The cut must be made
in the anchor frame, so it has to come after the rotation decision, and
C13 already writes a provenance record that collates into
`metadata.json`. In BIN_TARGET, the collapse would have to happen
before rotation, in Flye's frame, and the copy kept could split a gene.

## 7. Out of scope (add to tasks/todo.md on completion)

- A generic "target longer than `expected_target_size.max_bp`" finding.
  Task 69 left `max_bp` unused. A too-long target that is not a tandem
  duplicate (a genuine expansion, or a chimera) is a separate question.
- Choosing which copy to keep by annotation quality, or building a
  consensus of the copies.
- Tandem detection for `plant_pt` (C4 rebuilds the circle from the
  graph) and `plant_mt` (genuine repeats).
- Correcting Flye's per-contig coverage for the collapse. The reported
  contig coverage is Flye's figure for the doubled contig.
- Run-report marker for collapsed samples. The per-sample report is
  enough for now.
- EXTRACT_BARCODES emitting one locus twice from a single contig (the
  COX2 ×2 seen here). After this task, that can only come from a
  genuine duplication, but C5 has no guard of its own.

## 8. Acceptance criteria

- [x] §3.0 survey recorded in Outcomes, and the user told about any
      real `collapsed` or `unresolved` target.
- [x] C13 detects, collapses and records tandem duplicates as §2–§3.2
      describe. `rotation.json` is `wf5/rotation/v2`. 100% branch
      coverage, flake8 clean.
- [x] `tandem_collapse` is in `nextflow.config`, validated in
      `validateParams()`, passed per target, and marked provisional.
- [x] The report shows `collapsed` (info) and `unresolved` (warning)
      lines, with `unresolved` in the warnings mirror.
- [x] The integration assertions check `tandem.status` on all three
      fixtures, and the no-fixture gap is recorded.
- [x] §4.4 hand run on the doubled fixture: collapsed to about
      16,951 bp, COX1 at 1, each anchor gene once, and about 99 %
      identity or better to the single contig.
- [x] A full `-profile integration` run passes, and the run-report path
      is given to the user.
- [x] The spec, README and todo updates in §5 are made.

## Outcomes

### §3.0 Survey (2026-10-22)

`detect_tandem()` (params: k 31, min_support 0.5, tolerance 50,
min_period 1000) was run directly over every `animal_mt` and plant target
FASTA available. **No real sample was `collapsed` or `unresolved`.** Only
the task 68 doubled fixture is `collapsed`.

| Source | Targets | Result |
|---|---|---|
| Fixtures | doubled 33,905 bp: P 16,951, support 0.773; single 16,960 bp: P 124, support 0.023 | `collapsed`; `single` |
| Local `tests/integration/output`, `tests/manual/output` (incl. `az_*`) | 13 `animal_mt` targets (5,340–16,961 bp), 12 plant contigs | all `single`; animal support 0.000–0.041 |
| daff-admin `/mnt/data/wf5/output`, `output_round_3`, `output_round_3_t68`, `output_round_2_insect_t68`, `output_t76_barcode05` | R2/R3 insect samples (5,771–15,908 bp) | all `single`; support 0.000–0.041 |

The fixed `single` control measured 0.023, not the brief's 0.002. The brief's
figure counted only the dominant-offset positions once; this function counts
every position whose k-mer recurs within ±50 bp of the period. The gap to
0.5 is still more than 20×. The unit test bounds the control at < 0.03.

### Implementation notes

- `process_contig()` is now `rotate_record()` (unchanged rotation logic,
  split out) plus tandem detection and cut. `rotation.json` also carries a
  top-level `tandem_params` (null when disabled).
- Support counts a position when any other occurrence of its k-mer lies at
  circular distance P ± tolerance. K-mers occurring more than `MAX_COPIES`
  (4) times are excluded, which is what makes `(AT)n` runs `single`.
- Deviation from the brief's unit test list (no behaviour change): the
  segmental duplication test uses a 30 kb unique genome plus one repeated
  5 kb segment, not 5 kb of 15 kb. At 15 kb, the repeated k-mers cover
  enough of the circle that support is at or above 0.5 and the result
  depends on the exact fraction (it is `unresolved`, not `single`). The
  test shows that a partial duplication is never `collapsed`.
- The BIN_TARGET target-length bound in `assertions.sh` now measures the
  published `organelle_assembly.fasta`, so a collapsed contig is bounded
  by its output length. The rotation length check is now `input_length_bp`
  against `bin_metadata.json`, and any length change requires `collapsed`.
- Local Nextflow 26.10.3 cannot parse `nextflow.config` (task 72 strict
  syntax, shelved). The stub and integration runs used `NXF_VER=25.04.0`,
  as CI does.

### §4.4 Real-data hand run (doubled fixture)

Miniprot 0.18 biocontainer (pinned digest, `docker run --user`), COX1/COX2/
COX3/CYTB panels from `refs/v2026.10_1` under tables 2 and 5, then
`rotate_target.py` with the default tandem args:

- `rotated`, COX1, offset 6,992, `collapsed`, `copies` 2, `period_bp`
  16,951, support 0.773, `cut_method` `leading_kmer`, **cut 16,957 bp**
  (copy 2 carries a few indels, so the next leading k-mer lies 6 bp past
  P; within tolerance).
- Full `animal_mt` panel on the collapsed output (tables 2 and 5): every
  gene's best hit occurs once; no hit scoring ≥ 80 % of a gene's best lies
  more than 1 kb from it. COX1 is at 1–1563 on `+` (table 5, score 2,491).
- minimap2 `asm5` of the collapsed copy against `INT-ANIMAL-01.single.fasta`
  (both start at COX1): two co-linear alignments, 0-5,062 and 5,817-16,957
  of the collapsed copy, at 99.37 % identity (5,044/5,076 and 11,105/11,175
  matches). About 750 bp between them (5,062-5,817) did not align at
  `asm5`; this is a divergent stretch, not a rotation or duplication
  effect (the criterion of "a single near-full-length alignment" was not
  met literally). The aligned blocks are at the expected identity.

### Verification

- `scripts/pytest.sh`: full suite passes, 100 % branch coverage
  (`bin/rotate_target.py` 219 statements, 72 branches). flake8 clean on
  `bin/` and `scripts/tests/`.
- Stub run (`NXF_VER=25.04.0 nextflow run . -profile stub -stub-run`)
  passes with the new `val(tandem)` input.
- Full `-profile integration` run (with `-c tests/manual/local_resources.config`,
  because `process_medium` asks for 16 GB and this host has 15 GB; the
  first attempt stalled for that reason, unrelated to this task):
  `assertions.sh` exits 0, 151 `OK:` lines, no `FAIL:`. INT-ANIMAL-01:
  `tandem.status` `single`, rotated at COX1, 16,957 bp. INT-PLANT-01-pt:
  `not_tested`. INT-PLANT-01-mt is not in `ASSEMBLING_SAMPLES` (it stays
  a `low_coverage` sample, see task 35), so its tandem status is not
  asserted.
- Run report: `tests/integration/output/run-report.html` (per-sample
  reports alongside, e.g. `tests/integration/output/INT-ANIMAL-01/report.html`).
  The integration fixtures do not produce a doubled contig, so the
  `collapsed` report line is covered by the unit tests only.
