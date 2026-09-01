# Task 38 — `MINIPROT_CDS`: align under the target's genetic code

**Phase:** P4a (a correctness repair to the stage 12 pass built in task 30).
**Goal:** Stop miniprot translating organelle DNA with the standard genetic
code. Select the correct NCBI translation table per sample, pass it to
miniprot, and record the choice as provenance.

**Status:** implementation. This is a defect fix, not a design change —
except for §2, where `animal_mt` genuinely has no single correct answer
and a selection mechanism has to be chosen.

**Related:** repairs `modules/local/miniprot_cds.nf` as built in
task 30_unified_locus_pass.md. Should land **before**
task 39_annotate_rescue_annotator_only.md, which operates on the
cross-check buckets this task changes, and before
task 40_annotation_confidence_scores.md.

---

## 1. Overview

Stage 12 runs one miniprot pass over the organellar protein panel. Its
output, `cds.gff`, is the single source of protein-coding features for
both barcode extraction (stage 13) and the merged annotation (stage 13a)
— so anything wrong with that alignment is wrong everywhere downstream.

miniprot's `-T` option selects the NCBI translation table it uses to
conceptually translate the target DNA. It **defaults to table 1**, the
standard genetic code. `modules/local/miniprot_cds.nf` never passes
`-T`. Every organelle this pipeline assembles uses a different code:
`animal_mt` is table 2 or 5, `plant_pt` is table 11, `plant_mt` is
table 1 (the one case where the default happens to be right).

The practical damage on `animal_mt` is that **TGA codes for tryptophan,
not STOP**. Under table 1 miniprot reads every internal Trp as a
premature stop, so it truncates alignments, invents frameshifts to work
around the artefact, and in the worst case abandons a gene entirely.

Measured on a client `animal_mt` intercept (`CLIENT-BC05`, a
*Bactrocera dorsalis* mitogenome — the BLAST validation places it at
99.1 % nucleotide identity to `NC_008748.1`), re-running the identical
panel with `-T 5`:

| | table 1 (current) | table 5 |
|---|---|---|
| alignments emitted | 88 | 97 |
| ATP6 identity (best ref) | 0.7853 | **0.8540** |
| COX1 identity | 0.7988 | **0.8750** |
| COX3 identity | 0.8023 | **0.8588** |
| ND4 identity | 0.6440 | **0.7087** |
| ND6 identity | 0.5398 | **0.6087** |
| `Frameshift=2` flags on ATP6 | on 6 references | **none** |

The frameshift row is the one that matters most for
[principle 18](../CONSTITUTION.md). Those flags are the pipeline
asserting a disruptive mutation in a protein-coding gene of a
high-priority plant pest. They are false. The assembly is intact and the
ATP6 coordinates are identical under both tables — only the reported
quality of the alignment changes. A reviewer acting on that flag would
be acting on an artefact of a missing command-line argument.

Every `Identity=` value currently in a shipped `cds.gff`,
`annotation_summary.json` and `validation.tsv` is depressed by this, and
the `coordinate_conflicts` bucket of the CDS cross-check should be
re-read once the table is right.

**Not fixed by this task:** `ATP8` is still missed on this sample under
table 5 (zero alignments either way). That is a seeding limit, not a
translation problem, and it is task 39's subject.

**Exit criteria:**

- `MINIPROT_CDS` passes an explicit `-T` on every target.
- The table used is chosen per-sample for `animal_mt` (§2), emitted
  from the process, and recorded in `metadata.json`
  ([rule 16](../CONSTITUTION.md)).
- No stage downstream infers the genetic code from config when the real
  value is now available (§4).
- Integration bounds re-baselined against the corrected alignments (§6).

---

## 2. Which table — `animal_mt` cannot take a static value

`params.genetic_code_tables` in [nextflow.config](../nextflow.config)
already records the answer for two of three targets:

```groovy
genetic_code_tables = [
    animal_mt: [2, 5],   // vertebrate / invertebrate — trial both
    plant_pt:  [11],
    plant_mt:  [1],
]
```

`plant_pt` and `plant_mt` are single-valued: pass the one entry as `-T`
and this task is done for those branches.

`animal_mt` is a list because the kingdom gate
([constraint 2](../CONSTITUTION.md)) declares the organelle, not the
clade. Tables 2 and 5 differ at exactly one pair of codons:

| codon | table 2 (vertebrate) | table 5 (invertebrate) |
|---|---|---|
| `AGA` / `AGG` | **STOP** | Ser |
| `ATA` | Met | Met |
| `TGA` | Trp | Trp |

So a static choice is wrong for half the submissions, in both
directions. Table 2 on an invertebrate truncates at every Ser — the same
failure mode this task is fixing. Table 5 on a vertebrate reads straight
through genuine stop codons and over-extends the gene. Neither is
acceptable for an intercept workflow that does not know the clade in
advance.

**Recommended mechanism: trial both, keep the better pass.** Run
miniprot once per configured table and select between the results. This
mirrors what C5 already does per-locus in `bin/validate_barcodes.py`,
and miniprot is cheap enough that a second pass over a ~16 kb target and
a 130-sequence panel is not a meaningful cost.

**Selection criterion** (the implementer's call, but state it in the
code and the summary): sum the best per-gene alignment score across the
panel, tie-break on the number of distinct genes recovered. Both tables
score the same query set against the same target with the same matrix,
so the totals are comparable. On `CLIENT-BC05` this separates cleanly —
table 5 wins on both counts (97 vs 88 alignments, uniformly higher
identity). Do not use `Frameshift`/`StopCodon` counts alone; they are
symptoms, and using them as the criterion couples the selector to
miniprot's flag semantics.

**Rejected alternative — infer the clade from `BLAST_VALIDATE`.** Stage
11 already identifies the source organism to species level and could
name the right table directly. Rejected: it makes CDS annotation depend
on a validation stage it currently does not consume, and it fails
exactly where a biosecurity workflow most needs to succeed — a novel or
under-referenced intercept with no close BLAST hit. The trial has no
such blind spot.

**Rejected alternative — reuse `params.annotate[target].genetic_code`.**
That key exists only for `animal_mt` (it configures MITOS2) and hard-codes
5, so it reintroduces the vertebrate failure above and is undefined for
both plant branches.

## 3. Where the selection logic lives

Trial-and-select is more than ~10 lines of shell, so under
[rule 14](../CONSTITUTION.md) it belongs in `bin/`, not in the process
`script:` block. Add a small component — suggested
`bin/select_genetic_code.py` — that takes the candidate GFFs and emits
the winning one plus a JSON record of the decision:

```
# pseudocode
for each candidate table T:
    run miniprot --gff -T T  ->  cds.T.gff
select_genetic_code.py --candidate 2:cds.2.gff --candidate 5:cds.5.gff \
    --out-gff <sample>.cds.gff --out-json genetic_code.json
```

Single-entry targets skip the selector entirely — one miniprot call,
`-T` from config, no decision to record beyond the table used. Keep that
path free of the selector so `plant_pt`/`plant_mt` gain no new failure
surface from a change that does not concern them.

`genetic_code.json` should carry, at minimum: the chosen table, the
candidates considered, and the score/gene-count for each. That is the
audit trail for a decision made under uncertainty
([rule 16](../CONSTITUTION.md)) and it is what §4 consumes.

Per [rule 14](../CONSTITUTION.md) the new component needs
`scripts/tests/test_select_genetic_code.py` at **100 % branch
coverage**: single-candidate passthrough, two-candidate selection either
way, exact tie, empty candidate GFF (a target where miniprot found
nothing under one table), and malformed input. Fixtures should be
hand-written minimal GFF fragments, not captured tool output — drift in
test data is a standing maintenance cost
([rule 19](../CONSTITUTION.md)).

## 4. Downstream: retire the config proxy in `ANNOTATE`

`modules/local/annotate.nf` currently guesses the barcode-side genetic
code, with a comment saying so:

```groovy
// Config proxy for C5's per-run clade-trial choice — ANNOTATE has
// no dependency on EXTRACT_BARCODES, so it cannot see what table
// C5 actually picked (task 31 §3). Last entry of the configured
// trial order (params.genetic_code_tables) is the representative
// value recorded for comparison.
def barcodeCode = params.genetic_code_tables[meta.assembly_target][-1]
```

Once stage 12 selects and emits a real table, that guess is unnecessary:
`MINIPROT_CDS` is already upstream of both consumers, so the chosen
value can travel in `meta` or as a staged file.

**This changes what `genetic_code_agreement` means, and the task must
resolve it rather than leave a field whose semantics have quietly
shifted.** Today it compares MITOS2's configured table against a config
constant — an assertion about `nextflow.config`, not about the sample.
After this change the honest comparison is between the table stage 12
*selected* and the table MITOS2 was *configured* with (`params.annotate
[animal_mt].genetic_code`, fixed at 5). A vertebrate `animal_mt`
submission would then correctly report `false` — miniprot chose 2,
MITOS2 ran under 5 — which is a real and currently invisible QC signal,
and exactly the kind of disagreement
`tasks/todo.md`'s standing `RUN_REPORT` item expects to surface.

Consider whether MITOS2's `--code` should simply follow the selected
table too. That would make the pipeline self-consistent but removes the
independence that gives the cross-check its value; note the decision
either way.

## 5. Unit tests

Beyond §3's new module, `scripts/tests/test_annotate_summary.py` needs
its `genetic_code_agreement` cases reworked for the new comparison, and
should gain a case where the two sources genuinely disagree.

`MINIPROT_CDS` itself stays a tool wrapper — `-stub-run` for wiring,
integration for behaviour ([rule 15](../CONSTITUTION.md)). The stub must
emit `genetic_code.json` alongside the GFF or the wiring check will not
exercise the new output.

## 6. Integration reconciliation

Correcting the table moves real numbers, so `-profile integration` will
shift on `INT-ANIMAL-01` and `INT-PLANT-01-pt`:

- **Identity rises across the board.** Anything asserting an identity
  floor gets easier, not harder — but re-baseline rather than assume.
- **More alignments survive**, so `feature_counts.CDS` and recovered
  barcode loci may increase.
  `tests/integration/expected/animal_mt/annotation_bounds.json` sets
  floors (`min_cds: 8`), so it should not break, but the floors are now
  slack and should be tightened to keep their regression value.
- **`coordinate_conflicts` should shrink.** `max_cds_crosscheck_conflicts`
  is 2; re-measure and tighten.
- **The barcode identity floor interacts.** `EXTRACT_BARCODES` hard-codes
  `--min-identity 60`, and `tasks/todo.md`'s standing item records only
  1/6 `animal_mt` loci clearing it on `INT-ANIMAL-01`. Some of that
  shortfall is this bug. Re-measure the six loci and update the
  todo item with corrected figures — but **do not sweep the threshold
  here**; that is still blocked on the benchmark set
  ([spec §5.1](../spec/05-test-data.md)).

`plant_mt` is unaffected (table 1 either way) and is not currently
asserted downstream of the coverage gate anyway.

## 7. Spec updates

- [spec §2 stage 12](../spec/02-stages.md) — record that miniprot runs
  under an explicit table, and how it is selected on `animal_mt`.
- [spec §2.2](../spec/02-stages.md) — add the §3 component to the
  custom-logic table with its container and test surface.
- [spec §2.2 C7](../spec/02-stages.md) / provenance — add the selected
  table to the `metadata.json` contract.

## 8. Acceptance criteria

- No miniprot invocation anywhere in `modules/` runs without `-T`.
- `plant_pt` uses 11, `plant_mt` uses 1, `animal_mt` selects between 2
  and 5 per sample and records why.
- `metadata.json` carries the selected table and the candidates
  considered.
- `genetic_code_agreement` compares two independently-derived values and
  its meaning is documented where it is written.
- `scripts/tests/test_select_genetic_code.py` at 100 % branch coverage;
  `scripts/pytest.sh` green; `flake8` clean.
- `-profile stub -stub-run` green.
- `-profile integration` green with re-baselined bounds, and the
  `INT-ANIMAL-01` ATP6 `Frameshift` flags gone.

## 9. Out of scope

- **Recovering ATP8** — task 39. Table 5 does not fix it.
- **Sweeping `--min-identity`** — blocked on benchmark data
  ([spec §9](../spec/07-open-questions.md)).
- **Broadening the protein panel.** The panel's 10 reps/gene contain no
  tephritid, which is why `CLIENT-BC05` identities top out around 0.87
  even corrected. That is the standing panel-breadth item in
  `tasks/todo.md`, not this task.
- **Changing MITOS2's `--code`** beyond noting the decision in §4.

## Outcomes

Implemented as specified in §1–§5, with one structural deviation from
the §3 pseudocode and several bugs found and fixed during
`-profile integration` testing.

**§2 selection.** `bin/select_genetic_code.py` trials both configured
tables and sums the best per-gene miniprot alignment score (`mRNA`
column 6) across the panel, tie-breaking on distinct gene count, then
trial order — as specified. `Frameshift`/`StopCodon` counts are not
part of the criterion. `scripts/tests/test_select_genetic_code.py`
covers single-candidate passthrough, both selection directions, an
exact tie, an empty candidate, and malformed input, at 100% branch
coverage.

**Deviation from §3 — two Nextflow processes, not one.** The task's
pseudocode runs `select_genetic_code.py` from inside `MINIPROT_CDS`'s
own `script:` block. That doesn't work: the `miniprot` biocontainer
(`quay.io/biocontainers/miniprot`) carries no Python, discovered only
when `-profile integration` hit `select_genetic_code.py: No such file
or directory` / `env: can't execute 'python3'`. Per rule 14 ("custom
logic ... runs in its own container"), the fix is a second process,
`SELECT_GENETIC_CODE` (`modules/local/select_genetic_code.nf`), sharing
`EXTRACT_BARCODES`'s already-pinned `neoformit/daff-wf5-scripts` image
rather than a new build (rule 12). `MINIPROT_CDS` now only runs
miniprot, once per configured table, into `candidates/cds.<table>.gff`;
for single-table targets it also resolves `<sample>.cds.gff` /
`genetic_code.json` directly in plain shell (no Python needed, no new
failure surface for `plant_pt`/`plant_mt`). `main.nf` branches
`MINIPROT_CDS.out.candidates` on
`params.genetic_code_tables[meta.assembly_target].size() > 1`, routes
the `animal_mt` branch through `SELECT_GENETIC_CODE`, and rejoins both
branches (keyed on `meta`) into one `ch_cds` channel feeding both
`EXTRACT_BARCODES` and `ANNOTATE` — preserving the "one broad pass, two
consumers" invariant from spec §2 stage 12.

**§4 downstream.** `annotate.nf`'s config-proxy `barcodeCode` line is
gone. `genetic_code.json` is passed straight through to
`annotate_summary.py`'s new `--genetic-code-json` flag (a staged
`Path`), which reads `selected_table` itself — not parsed in Groovy;
an earlier attempt to do `new groovy.json.JsonSlurper().parseText(
genetic_code_json.text)` inside the `.nf` script block threw
`java.nio.file.ProviderMismatchException`, and pushing the parse into
Python sidesteps it entirely. `genetic_code_barcodes` was renamed to
`genetic_code_cds` throughout (`annotate_summary.py`'s arg, `run()`
param, and the `annotation_summary.json` field) to name what it now
actually is: the table `MINIPROT_CDS`/C9 selected for this sample, not
a config-constant proxy for EXTRACT_BARCODES's independent trial.
**Decision on MITOS2's `--code`:** left as-is (fixed at `params.annotate
.animal_mt.genetic_code`, currently 5) rather than following C9's
selection — a vertebrate sample where miniprot selects table 2 while
MITOS2 stays configured for 5 now surfaces as `genetic_code_agreement:
false`, a real cross-check signal per §4's framing, which following C9
would erase.

**§6 integration reconciliation**, measured on `INT-ANIMAL-01`
(*Acyrthosiphon pisum*, invertebrate — miniprot correctly selected
table 5, `total_score` 8170 vs. 7043 for table 2):

- `annotation_bounds.json` tightened to measured values: `min_cds` 8→9,
  `min_trna` 15→20, `min_rrna` 1→2, `max_cds_crosscheck_conflicts` 2→0.
- Direct table-1-vs-table-5 comparison on the same target+panel: 62→74
  miniprot alignments, `StopCodon` flags 123→26, `Frameshift` flags
  125→148. The false-stop artifact this task exists to fix is
  substantially cleared (not zeroed — a handful of very-low-identity
  panel hits still trip it). **`Frameshift` did not disappear on
  `ATP6`** the way it did on the task's `CLIENT-BC05` anecdote: even
  under the now-correct table 5, `ATP6`'s best panel hit is only 48.4%
  identity (no aphid in the panel — same standing panel-breadth gap as
  `tasks/todo.md`'s barcode-floor item), and a `Frameshift=3` at that
  divergence is plausibly a genuine indel, not a translation artifact.
  What *did* clear on that exact hit: its co-occurring `StopCodon=3`
  flag is gone under table 5. Recorded in `tasks/todo.md` rather than
  claimed as met, since it doesn't hold on this fixture as stated.
- `tasks/todo.md`'s barcode-floor item re-measured: 2/6 `animal_mt`
  loci now clear the 60% identity floor on `INT-ANIMAL-01` (was 1/6),
  identities now 48.4–72.8% (was 42–65%). Updated in place with the
  full per-locus breakdown; the underlying panel-breadth diagnosis is
  otherwise unchanged, since it is not this task's fix.
- `plant_pt`/`plant_mt` unaffected as predicted (`plant_mt` still
  soft-fails the coverage gate before reaching this stage; `plant_pt`
  ran under its single configured table 11 throughout).

**Not implemented — `metadata.json`.** `COLLATE` (C6) is still the P0
stub (`touch metadata.json`); there is no real per-sample bundle logic
to add the selected table to yet. Recorded the contract in spec §2.2's
C6 row instead (`genetic_code.json` → `metadata.json`'s provenance,
rule 16) for whichever task implements `COLLATE` for real to pick up.

**Bugs found and fixed during implementation/integration testing** (in
the order hit, five full `-profile integration` cycles):

1. `bin/select_genetic_code.py` was missing its executable bit —
   `chmod +x`.
2. `annotate.nf` reading `genetic_code_json.text` in Groovy threw
   `ProviderMismatchException` — fixed by passing the path to
   `annotate_summary.py` instead (see §4 above).
3. The `miniprot` biocontainer has no Python — fixed by splitting
   `MINIPROT_CDS`/`SELECT_GENETIC_CODE` into two processes (see
   structural deviation above).
4. `MINIPROT_CDS`'s single-table `genetic_code.json` heredoc emitted a
   literal `\$schema` (invalid JSON) — a Groovy triple-quoted string
   needs `\$schema` (one backslash) to produce a literal `$`, not the
   `\\\$schema` first written.
5. **Nextflow quirk:** `path("name", optional: true)` inside a `tuple`
   output does not suppress the missing-file error per-path in Nextflow
   25.10.2 — confirmed with a minimal standalone repro. `optional: true`
   has to sit at the tuple level (after the last `path(...)`, before
   `emit:`). `modules/local/miniprot_cds.nf`'s `resolved` output uses
   the corrected form; `modules/local/collate.nf` has the same
   (currently dormant, since `COLLATE`'s stub always creates every
   file) latent bug, flagged in `tasks/todo.md` for whoever implements
   `COLLATE` for real.

**Verification:** `flake8 bin/ scripts/` clean; `scripts/pytest.sh`
green, 254 tests, 100% branch coverage repo-wide; `-profile stub
-stub-run` green (exercises both the single-table and multi-table
`MINIPROT_CDS`/`SELECT_GENETIC_CODE` wiring); `-profile integration`
green, `tests/integration/assertions.sh` all-pass.
