# Task 76 — Rotate circular `animal_mt` assemblies to start at COX1

**Phase:** P5 — production deployment (spec/06-phases.md).

**Status:** complete. No preconditions. Pick up ahead of the other open
tasks: the nightly integration run fails intermittently until this lands
(see Overview).

**Related spec:** spec/01-pipeline-flow.md §1 (stage diagram) and §1a,
spec/02-stages.md stage table and §2.2 (custom logic components),
spec/03-organelles.md §3.6 (plastid canonicalisation, the `plant_pt`
precedent) and §3.7.4 (circularity comes from Flye),
spec/04-reference-data.md §4.3 (protein panel), spec/05-test-data.md
§5a/§5b, spec/06a-reports.md §6a.2. CONSTITUTION principles 8 and 10,
rules 11, 12, 14, 16, 17, 18 and 19.

## Overview

An animal mitochondrial genome is a circle of about 15–20 kb. Flye
closes the circle and writes it out as a linear FASTA record, but where
it cuts the circle to do that (the *origin*, position 1) is arbitrary.
The cut can also change between two runs on the same reads. Flye runs
multithreaded and is not bit-reproducible, and nothing downstream of it
in this pipeline puts the sequence back into a standard frame.

That causes two problems:

1. **A gene can be cut in two.** If the origin lands inside a gene, half
   the gene sits at the end of the record and half at the start. Every
   tool after assembly (miniprot, MITOS2, the barcode extractor, the
   organelle map) treats the record as linear, so each tool sees two
   partial genes or misses the gene.
2. **Coordinates are not comparable between runs.** The same genome
   comes out with every feature shifted by thousands of bp, so outputs
   from two runs, or two samples of the same species, can't be compared
   by position.

This has already happened. The nightly integration run on 2026-10-20
(GitHub Actions run 39385848346) failed `INT-ANIMAL-01 cds_crosscheck
conflicts 1 exceeds max 0`. The fixture and code were the same as in a
run that passed earlier that day (run 39365421908). Comparing the two
outputs:

| | Passing run | Failing run |
|---|---|---|
| Target contig | 16,957 bp, circular | 16,960 bp, circular |
| COX1 (miniprot) | 13766–15293 | 8694–10221 |
| ND6 (miniprot) | 4976–5465, whole | 3–389, at the contig start |
| nad6 (MITOS2) | 4978–5040 | 16862–16924, at the contig end |

The origin moved by about 5 kb and landed inside ND6. MITOS2 called the
piece of ND6 at the end of the contig and miniprot called the piece at
the start. `bin/annotate_summary.py` (C8) found no overlap between the
two calls and recorded a coordinate conflict. The biology was the same
in both runs. Only the coordinate frame differed.

Pinning a random seed does not fix this. Flye has no seed option on this
path. The one random subsample before assembly (`seqtk sample` in
COVERAGE_GATE) is already seeded with `params.seqtk_seed`. The variation
comes from Flye itself.

The standard fix is to **rotate** each circular mitogenome to a fixed
start, as GenBank records and mitogenome tools do. Most invertebrate
mitogenomes, which make up this pipeline's intake, start at the first
base of COX1 on the forward strand. This task adds a step after
BIN_TARGET that finds COX1 on each circular `animal_mt` target contig,
reverse-complements the contig if COX1 is on the minus strand, and
rotates it so COX1 starts at position 1. Every later stage then works on
a sequence whose position 1 is a gene boundary, and every run lands in
the same frame. If COX1 can't be found, the step falls back to the start
of another conserved gene (§3.3). If no anchor gene is found at all, it
leaves the contig as it is and says so in the provenance and the report.

## 1. Constraints

- **`animal_mt` only.** `plant_pt` already gets a deterministic frame
  from C4's canonical `path1` (spec §3.6). `plant_mt` assemblies are
  rarely one circle and have no conserved gene order. Both pass through
  unchanged, with the reason recorded (§3.4). This is driven by config
  (§3.1), not by a branch on `assembly_target` in `main.nf`.
- **Circular contigs only.** Rotating a linear contig is meaningless: its
  ends are real breaks, not an arbitrary cut. Circularity comes from
  `bin_metadata.json` (`circular_method` is `flye_circ` or `end_overlap`;
  spec §3.7.4), not from a fresh test.
- **Rotation is a pure permutation.** Same length, same bases. The only
  operations are reverse-complement and cyclic shift. No trimming,
  padding or re-joining. The doubled contigs of task
  75_circular_tandem_duplicate_contigs.md are explicitly left to that
  task. Rotating such a contig is harmless and does not hide the
  problem.
- **Contig IDs are unchanged.** Every later join keys on the contig name
  (`contig_10`). The FASTA header line is written back exactly as read.
  Provenance lives in `rotation.json`, not in the header.
- **Rotate before annotation, not after.** Shifting `cds.gff`
  coordinates after the fact can't repair a gene the origin already
  split, because miniprot only aligned the part it could see. The step
  must run before BLAST_VALIDATE and MINIPROT_CDS.
- **Never fail a sample over rotation.** Anything short of a
  configuration error gives a passthrough and a recorded status
  (principle 8). A configured anchor gene with no panel file in the
  reference bundle *is* a configuration error, and fails loudly (rule
  18).
- **Containers.** Follow the MINIPROT_CDS → SELECT_GENETIC_CODE pattern:
  miniprot runs in the existing pinned miniprot biocontainer, and the
  Python logic runs in the existing `wf5/scripts` image (rules 11, 12).
  No new image.
- **Panel data comes from the versioned reference bundle**
  (`proteins/animal_mt/<GENE>.faa`, spec §4.3; principle 10). No new
  reference data is added.

## 2. Definitions

- **Anchor gene:** a protein-coding gene whose start defines position 1.
  The anchors are an ordered list per target in `nextflow.config`, with
  COX1 first.
- **Anchor hit:** the best miniprot alignment of an anchor gene's panel
  proteins to a given contig. Best means highest miniprot score, with
  ties broken by lowest start coordinate, then by table, so the choice
  is deterministic.
- **Start estimate:** the first base of the anchor's CDS on the anchor's
  strand. miniprot often starts aligning a few residues into the
  protein, because panel N-termini vary and insect COX1 start codons are
  often non-canonical. So the start estimate extends upstream by
  `3 × (query_start − 1)` bp, where `query_start` is the first aligned
  residue of the winning panel protein (from the `Target=` attribute).
  The extension wraps around the circle, and the raw and adjusted
  positions are both recorded.

## 3. Scope

### 3.1 Config (`nextflow.config`)

Add a per-target anchor list beside `genetic_code_tables`:

```
rotation_anchors = [
    animal_mt: ['COX1', 'COX2', 'COX3', 'CYTB'],
    plant_pt:  [],     // C4 path1 frame (spec §3.6)
    plant_mt:  [],     // no single circle, no conserved order
]
```

The fallback order shown is a starting point: long, conserved genes that
miniprot finds reliably. ATP6 is left out because of its frameshift
over-extension (tasks/todo.md, Annotation). Confirm the order against
`INT-ANIMAL-01` with COX1 removed from the panel (§4.3), and record the
evidence.

### 3.2 New process `ROTATE_ANCHORS` (stage 9a, miniprot container)

Input `(meta, target.fasta)` plus `protein_panel`. For each gene in
`params.rotation_anchors[meta.assembly_target]`, check that
`${protein_panel}/<target>/<GENE>.faa` exists. If one is missing, exit
non-zero with a message naming the gene and the bundle path. Then run
miniprot `--gff` against the concatenated anchor panels under each table
in `params.genetic_code_tables[target]`, writing
`anchors/anchors.<T>.gff`. An empty anchor list or an empty
`target.fasta` writes empty GFFs and succeeds. The total work is about
4,400 queries against about 17 kb, which is negligible next to
MINIPROT_CDS.

The panel is passed as an `input: path` value channel, never as
`${params.x}` (rule 13). Add a `withName` entry in
`conf/containers.config` with the same digest as `MINIPROT_CDS`. The
`stub:` block creates every declared output.

### 3.3 New process `ROTATE_TARGET` (stage 9a, C13 `bin/rotate_target.py`, `wf5/scripts` image)

Input: `(meta, target.fasta, secondaries.tsv, bin_metadata.json,
anchors/*.gff)` plus the anchor list (passed as a `val`, not read from
`params` inside the script). Output: `(meta, target.fasta,
secondaries.tsv)` with the same shape as `BIN_TARGET.out.binned`, plus
`(meta, rotation.json)`.

Pseudocode per contig in `target.fasta`:

```
if no anchors configured:                 status = not_configured
elif target.fasta empty:                  status = no_target
elif contig not circular in bin_metadata: status = not_circular
else:
    for gene in anchors (in order):
        hit = best hit of gene on this contig (across tables)
        if hit: break
    if no hit:                            status = no_anchor
    else:
        seq = revcomp(seq) if hit.strand == '-'   # map hit coords too
        start = adjusted start estimate (§2), mod length
        seq = seq[start:] + seq[:start]
        status = rotated  (fallback = gene != anchors[0])
write contig with original header line
```

Contigs are rotated independently. A secondary circular contig with no
anchor hit passes through with `no_anchor`. `secondaries.tsv` passes
through untouched. The script always exits 0 except on unreadable input.

`rotation.json` (schema `wf5/rotation/v1`) has one record per contig:
`contig_id`, `status`, `anchor_gene`, `fallback` (bool),
`anchor_table`, `reverse_complemented`, `raw_start`, `query_start`,
`offset` (bases moved to the front), `length_bp`, and the `anchors`
list as configured. It also has a top-level `status`: the primary
contig's status. The primary contig is the longest one, the same rule
the report already uses.

### 3.4 Wiring (`main.nf`)

```
BIN_TARGET.out.binned ─► ROTATE_ANCHORS ─┐
BIN_TARGET.out.binned + .metadata ───────┴► ROTATE_TARGET ─► BLAST_VALIDATE ─► …
```

- Every consumer of `BIN_TARGET.out.binned` downstream of stage 9 takes
  `ROTATE_TARGET.out.binned` instead. That includes `BLAST_VALIDATE` and
  the COLLATE input join (currently `main.nf`'s
  `.join(BIN_TARGET.out.binned, by: 0)`), so the deliverable
  `organelle_assembly.fasta` is the rotated sequence.
- `BANDAGE_NG` keeps reading Flye's graph, which has no linear
  coordinates, so nothing changes there.
- `ORGANELLE_MAP` already draws from the GFF made on the rotated target.
  COX1 now sits at 12 o'clock.
- `rotation.json` joins into COLLATE.

### 3.5 Provenance and report

- **C6 (`bin/collate.py`):** copy `rotation.json`'s per-contig records
  into `metadata.json` under `rotation`, and extend the sample metadata
  schema to match (rule 16). The bundle also carries `rotation.json`
  under `diagnostics/`.
- **Per-sample `report.html`, Assembly tab (spec §6a.2):** one line per
  primary contig. "Rotated to start at COX1 (reverse-complemented)."
  "Rotated to start at COX2 — COX1 not found." "Not rotated: no anchor
  gene found; coordinates follow Flye's arbitrary origin." "Not rotated:
  contig is linear." Show nothing for `not_configured` targets.
  `no_anchor` on a circular contig gets the same warning styling as
  other assembly caveats. It isn't a failure, but a reader comparing
  coordinates needs to know.

## 4. Testing

### 4.1 Unit tests (`scripts/tests/test_rotate_target.py`)

C13 needs 100 % branch coverage (rule 14). Use small synthetic FASTA and
GFF fixtures written in the test, not reference data (rule 19). Cover:

- Plus-strand anchor: correct offset, and the anchor CDS starts at 1.
- Minus-strand anchor: reverse-complemented, and the anchor is on the
  plus strand at 1.
- Anchor already at position 1: offset 0, `rotated`.
- `query_start > 1` adjustment, including an adjustment that wraps
  across the origin.
- COX1 absent with COX2 present gives `fallback: true`. COX1 present
  wins over a higher-scoring COX2 (order beats score).
- Several hits of one gene: highest score wins, then lowest start, then
  table.
- No hits gives `no_anchor`. A linear contig gives `not_circular`. An
  empty FASTA gives `no_target`. An empty anchor list gives
  `not_configured`.
- Multi-contig FASTA: each contig is handled independently. A hit on a
  contig absent from the FASTA is ignored.
- Invariants as a property check: same length, same base composition,
  and a rotated sequence (after reverse-complement where applied) is a
  cyclic shift of the original. Header lines are byte-identical.
- The `rotation.json` shape validates against `wf5/rotation/v1`.

C6 and report tests: extend `test_collate.py` and `test_report.py` for
the `rotation` block and each Assembly-tab message.

### 4.2 Stub run

Add both processes to the stub path. Extend the layout assertions in
`.github/workflows/tests.yml` so `STUB-01` carries
`diagnostics/rotation.json`. The container-coverage and bare-params lint
checks must pass.

### 4.3 Integration (`tests/integration/assertions.sh`)

`INT-ANIMAL-01` exercises the COX1 path. No new fixture is needed. Add:

- `rotation.json`: `status == rotated`, `anchor_gene == COX1`,
  `fallback == false`.
- `length_bp` equals the target contig's `length_bp` in
  `bin_metadata.json`.
- In `cds.gff`, the winning COX1 feature is on `+` and starts within a
  small tolerance of position 1. Set the tolerance from the measured
  value, with the evidence recorded, and expect tens of bp at most.
  Don't assert `reverse_complemented`, because that depends on Flye's
  output strand.
- The existing `max_cds_crosscheck_conflicts: 0` stays and should now be
  stable.

The plant fixtures assert `status == not_configured`.

Fallback evidence: run C13 by hand on the `INT-ANIMAL-01` target with
COX1 removed from the anchor list, and record which gene won and where
position 1 landed. This is evidence only. Don't add a permanent fixture.

### 4.4 Real-data acceptance

Because the original failure was nondeterministic, one green run proves
little. Run `-profile integration` twice, using the `remote-run` skill,
and record:

- Flye's raw origin for `INT-ANIMAL-01` in each run.
- The rotated COX1 start in each run.

COX1 must land at about position 1 in both runs, even if the raw
origins differ. Also run one round-2/3 client `animal_mt` sample, and
record the anchor used and the conflict count.

## 5. Spec and doc updates

- spec/01-pipeline-flow.md §1: add stage `[9a] ROTATE_ANCHORS →
  ROTATE_TARGET` on BIN_TARGET's branch, before BLAST_VALIDATE. Keep the
  existing stage numbers.
- spec/02-stages.md: add a stage-table row for 9a, and a C13 row in
  §2.2.
- spec/03-organelles.md: add a new §3.8, "Animal mitogenome rotation".
  It covers why the frame is arbitrary (with the 2026-10-20 evidence
  above), the COX1 / forward-strand convention, the fallback list, the
  passthrough statuses, and why rotation runs before annotation.
- spec/06a-reports.md §6a.2: add the Assembly-tab rotation line.
- docs/testing.md: add the new integration assertions.
- tasks/todo.md, Recruitment, the item "INT-ANIMAL-01's assembly length
  is not stable…": the closing sentence says the subsample seed needs
  pinning first. Correct it. The seed is already pinned, and Flye varies
  between runs regardless. Rotation stabilises coordinates, not
  assembly length.

## 6. Design decisions

### 6.1 Fallback: next anchor gene, not a gene-offset matrix

One alternative was a matrix of each gene's typical offset from the
COX1 start. With it, if COX1 isn't found but, say, CYTB is, position 1
could be set at "CYTB start minus the usual CYTB-to-COX1 distance". It
was rejected for this task:

- **It assumes a conserved gene order and spacing.** The arthropod
  ground pattern holds for most insects. But rearrangements are common
  in groups that matter for biosecurity (Hymenoptera, Thysanoptera,
  Phthiraptera, many Hemiptera), and vertebrate order differs outright.
  `animal_mt` covers all of Metazoa (spec §4.3), so the matrix would
  need per-clade tables that drift as the panel grows (rule 19).
- **An approximate origin can split a gene.** That reintroduces the bug
  this task fixes, and hides it behind a confident-looking frame (rule
  18). Rotating to the start of a gene that was actually found
  guarantees position 1 is a gene boundary.
- **The fallback is visible.** "Rotated to start at COX2 — COX1 not
  found" tells the reader exactly which frame they are in. An
  approximated COX1 start reads as if COX1 had been found.

The cost is that a COX1-less sample isn't in the COX1 frame. That's
acceptable: such a sample is already flagged by protein-coding
completeness. Revisit if real data shows COX1 failing often while order
is conserved. Even then, the estimate should snap to the nearest
annotated gene boundary rather than land wherever the offset puts it.

### 6.2 A separate miniprot pass, not reuse of `cds.gff`

Rotation has to come before MINIPROT_CDS (§1), so it can't read that
stage's output. Running MINIPROT_CDS twice (once to find COX1, again
after rotating) would double the most expensive annotation step. A
pass over only the anchor panels costs almost nothing.

## 7. Out of scope (add to tasks/todo.md on completion)

- **A circular-aware C8 crosscheck.** This would let an annotator call
  and a miniprot call on either side of the origin count as agreeing. It
  is still worth having as a safety net for `no_anchor` and `plant_mt`
  circles, but it isn't needed once rotation lands.
- **A conventional start for `plant_pt`** (`trnH-GUG` / start of LSC).
  C4's `path1` is deterministic given the graph, but not in the GenBank
  frame.
- **Detecting and trimming doubled contigs.** That is task
  75_circular_tandem_duplicate_contigs.md.

## 8. Acceptance criteria

- [x] `rotation_anchors` is in `nextflow.config`. `ROTATE_ANCHORS` and
      `ROTATE_TARGET` are wired as in §3.4, with containers pinned and
      stub blocks complete.
- [x] C13 has 100 % branch coverage. flake8 is clean.
- [x] The deliverable FASTA, `cds.gff`, annotation, barcodes and map are
      all in the rotated frame. `metadata.json` and the report show the
      rotation.
- [x] The integration assertions in §4.3 pass. Two remote integration
      runs show COX1 at about position 1 in both (§4.4), with evidence
      recorded in Outcomes.
- [x] The fallback evidence from §4.3 is recorded.
- [x] The spec and todo updates in §5 are made. The follow-ups in §7 are
      added to tasks/todo.md.

## Outcomes

### Implemented (2026-10-21)

- `bin/rotate_target.py` (C13): 25 unit tests, 100 % branch coverage,
  flake8 clean. Full suite 920 passed.
- `modules/local/rotate_anchors.nf`, `modules/local/rotate_target.nf`,
  `rotation_anchors` in `nextflow.config`, container pins in
  `conf/containers.config`, wiring in `main.nf`. `ROTATE_TARGET.out.binned`
  now feeds `BLAST_VALIDATE` and the COLLATE join.
- C6 (`bin/collate.py`) copies `rotation.json` into `metadata.json` under
  `rotation` and into `diagnostics/`; schema, Assembly-tab lines and the
  `no_anchor` warning are in `bin/report/report.py` and
  `scripts/report/templates/components/assembly.html`, with tests.
- Stub run passes (Nextflow 25.10.2); `tests.yml` asserts
  `STUB-01/diagnostics/rotation.json`. The stub COLLATE now copies
  `rotation.json` into `diagnostics/` for that purpose.
- Integration assertions added to `tests/integration/assertions.sh`
  (`ROTATION_COX1_TOLERANCE_BP`).
- Spec §1, §2 (stage 9a, C13), §3.8, §6a.2, docs/testing.md and todo.md
  updated; §7 follow-ups added to todo.md.

### Design notes / deviations

- `ROTATE_TARGET` takes the anchor list as a `val` element inside the
  input tuple (built per sample in `main.nf`), not a separate value input,
  because the list depends on `meta.assembly_target`.
- `ROTATE_TARGET` writes `rotated/target.fasta` (the input is already named
  `target.fasta`); consumers stage it by basename, so downstream is
  unchanged.
- Locally installed Nextflow 26.04 cannot parse `nextflow.config` (existing
  `if` inside `plugins {}`); runs here used `NXF_VER=25.10.2`.

### Evidence (local, existing INT-ANIMAL-01 target, refs v2026.10)

- Ran miniprot (pinned biocontainer) over the COX1/COX2/COX3/CYTB panels
  under tables 2 and 5 and C13 by hand on the target (contig_10, 16,799 bp,
  `flye_circ`). COX1 won (table 5), reverse-complemented, offset 12,378.
  Re-running miniprot on the rotated contig puts the best COX1 hit at
  1–1559.. on `+`, i.e. starting at position 1.
- Fallback: with COX1 removed COX2 won (offset 14,007) and re-annotation
  puts COX2 at position 1 `+`; with COX1 and COX2 removed COX3 won (offset
  15,602), COX3 at position 1 `+`.

### Evidence (local `-profile integration`, INT-ANIMAL-01 only, 2026-10-21)

- First run failed `ROTATE_TARGET` with exit 126: `bin/rotate_target.py`
  was not executable. Fixed with `chmod +x`; the resumed run completed.
- Target contig_10, 16,957 bp, `flye_circ`. `rotation.json`: `rotated`,
  COX1, `fallback: false`, table 5, `reverse_complemented: true`,
  `raw_start` 1628, `query_start` 1, `offset` 1627, `length_bp` 16957.
- Winning COX1 in `cds.gff`: 1–1559 on `+` (starts at position 1).
- `cds_crosscheck.coordinate_conflicts`: 0.
- The new rotation assertions pass. `ROTATION_COX1_TOLERANCE_BP` set to 30
  (measured start offset 0); the full assertions script was not run
  because only the animal sample was run locally.

### Evidence (remote, round-2 client sample, 2026-10-21)

- R2-INSECT-barcode05 (`animal_mt`, refs v2026.10) on daff-admin, 12m38s,
  all 26 tasks succeeded. Previously `ok` with 1 `coordinate_conflicts`.
- contig_4, 15,908 bp, circular: `rotated`, anchor COX1, no fallback,
  not reverse-complemented, `raw_start` 1576, offset 1575.
- COX1 in `organelle_annotation.gff`: 1–1533 on `+`. Coordinate conflicts
  now 0; status `ok`. The report carries the rotation line.
- Output: tests/manual/output/t76_barcode05/ (gitignored).

### Two remote integration runs (§4.4, 2026-10-21)

INT-ANIMAL-01 on daff-admin, two independent runs, refs v2026.10, 13 min
each, all 26 tasks succeeded, `ok`, 0 coordinate conflicts both times.

| | Run 1 | Run 2 |
|---|---|---|
| Flye target contig (pre-rotation) | contig_10, 16,960 bp | contig_10, 16,961 bp |
| Raw first 40 bp | `GGTATTATAAATATATTTAAT...` | `TTTTTTCATTAATTATTTTTTG...` |
| Rotation | COX1, reverse-complemented, offset 6702 | COX1, reverse-complemented, offset 6589 |
| Rotated first 40 bp | `ATTTACAATTTACAACCTATT...` | `ATTTACAATTTACAACCTATT...` |

Flye's raw origin differed between runs (different starting bases and a 1 bp
length difference), and after rotation both start identically with COX1 at
position 1 on `+` (best COX1 hit in `cds.gff` starts at 1 in both runs).
Output: tests/manual/output/t76_int1/ and t76_int2/ (gitignored).

The earlier "Deviation from §4.4" is superseded: both remote runs are done.
