# Task 42 — Stage 15 `COLLATE`: per-sample bundle + `metadata.json` (C6)

**Phase:** P4a (spec 06-phases.md).
**Goal:** Replace the P0 stub in `modules/local/collate.nf` with a real
`bin/collate.py` (C6) that classifies each sample, assembles the
per-sample output bundle, and emits `metadata.json`.

**Scope boundary:** this task ships the **data layer only**. The
`report.html` render is task 43_per_sample_report.md; `COLLATE` continues
to `touch` an empty `report.html` until then, so the process contract
does not change shape twice.

**Related tasks:** 43_per_sample_report.md and 45_run_report.md both
consume `metadata.json` and should not start before this lands.
44_organelle_map.md is independent.

---

## 1. Overview

Fifteen stages of this pipeline each drop their own outputs into their
own work directory: a coverage estimate here, a binned contig there, a
GFF, a barcode FASTA, two NanoPlot report directories, a BLAST table. A
biosecurity officer does not want fifteen directories. They want one
folder per sample containing the assembled organelle, its annotation, the
barcodes destined for Taxodactyl, a human-readable report, and a machine-
readable record of how all of it was produced.

`COLLATE` is the stage that produces that folder. It is the pipeline's
join point for a single sample — the first and only place that sees every
upstream output at once — and therefore the only place that can answer
the question the whole workflow exists to answer: **what happened to this
sample?**

That question has more than two answers. A sample can be a clean success;
it can be a real but incomplete success (assembled under the warn
coverage floor); it can have failed at the coverage gate before assembly
was ever attempted; it can have assembled nothing recognisable; or it can
have assembled an organelle from which no barcode locus could be
extracted. CONSTITUTION.md principle 7 requires that these be five
*distinguishable* signals, not one boolean — a degraded sample must never
produce output indistinguishable from a confident negative, and equally,
a partial recovery must never be reported as a failure merely for being
incomplete.

`COLLATE` today does none of this. It is the P0 scaffold stub: it
`touch`es five empty files unconditionally, for every sample, whatever
happened upstream. Everything in this task is the work of replacing that
`touch` with a decision.

The second half of the task is `metadata.json`. This is simultaneously
the Taxodactyl handoff record (brief.md §5), the biosecurity audit trail
(CONSTITUTION.md rule 16), and the input contract for both report tiers.
Building it first — before either renderer — is deliberate: it means
tasks 43 and 45 are written against a fixed, schema-validated structure
rather than co-evolving with it.

---

## 2. What `COLLATE` receives today — and what is missing

`main.nf` already joins ten upstream channels into `COLLATE`. Before any
of the logic below can be written, four plumbing defects have to be
fixed, because the stage cannot currently see data it is specified to
report on.

### 2.1 Current input tuple

`main.nf`'s `ch_ok_inputs` supplies, per sample: `status_json`,
`coverage_json`, `nanoplot_raw`, `nanoplot_clean`, `target_fasta`,
`blast_tsv`, `barcodes_fasta`, `coords_gff`, `validation_tsv`,
`annotation_gff`, `annotation_summary`, `organelle_map_svg`, and a
literal `[]` placeholder for `graph_png`.

`ch_failed_inputs` supplies the first five and nine `[]` placeholders.

### 2.2 Four defects to fix in this task

1. **`BIN_TARGET.out.metadata` carries no `meta`.**
   `modules/local/bin_target.nf` declares it as a bare
   `path("bin_metadata.json")`, so it cannot be `join`ed per-sample and
   is not wired into `COLLATE` at all. `bin_metadata.json` holds
   `contigs_selected`, `target_source`, `sibling_carryover`,
   `ref_identity_pct`, `circular`, and the whole
   `plastid_canonicalisation` block — all of it specified as report
   content (spec 06a-reports.md §6a.2, spec 03-organelles.md §3.6).
   Add `tuple val(meta),` to the emit and join it.
2. **`BIN_TARGET.out.isoforms` likewise carries no `meta`.** Same fix.
   Task 44_organelle_map.md needs it too.
3. **`BANDAGE_NG`'s graph PNG is discarded.**
   `BANDAGE_NG.out.assembly` does carry
   `path("${meta.sample_id}.graph.png")`, but `main.nf` passes a literal
   `[]` into `COLLATE`'s `graph_png` slot. Wire the real channel.
4. **`METAFLYE`'s `assembly_info.txt` never reaches `COLLATE`.** It is
   the source of the per-contig length / coverage / circularity table
   that spec 06a-reports.md §6a.2's assembly section specifies. Thread
   it through.

Fixing 1–4 widens `COLLATE`'s input tuple. Do it in this task rather
than task 43 — the renderer should not be the thing that discovers its
inputs are unreachable.

---

## 3. Sample classification: one status, derived from three sources

This is the core logic and the part most worth getting right.

Three upstream files each carry their own status field, with overlapping
but **non-identical** vocabularies. `COLLATE` must not simply copy one of
them forward.

| Source file | Its `status` values | What it describes |
|---|---|---|
| `sample_status.json` (C2) | `ok`, `low_coverage`, `fail` | the coverage gate's decision |
| `annotation_summary.json` (C8) | `no_assembly`, `ok_cds_only`, `annotator_failed`, `ok`, `no_features` | the annotation outcome |
| `validation.tsv` (C5) | per-locus `pass` / `fail` / `not_found` | barcode recovery |

`COLLATE` derives a single **sample-level** status from all three. Keep
the derived field distinctly named (e.g. `sample_status`) so it can never
be confused with the gate's `status` when both appear in one JSON
document.

### 3.1 The dispatch matrix

| Derived `sample_status` | Trigger | Bundle |
|---|---|---|
| `fail` | `sample_status.json.status == "fail"` | **minimal** |
| `no_assembly` | gate passed, but no target contig was selected | **minimal** |
| `no_barcode` | assembly exists, but zero loci passed validation | **full** |
| `low_coverage` | `sample_status.json.status == "low_coverage"`, assembly and ≥1 locus exist | **full** |
| `ok` | gate `ok`, assembly and ≥1 locus exist | **full** |

Two rules that are easy to get backwards:

- **`low_coverage` is a full bundle, not a failure path.** Assembly,
  annotation and barcodes all genuinely exist; the warning rides along in
  `metadata.json` so both report tiers and the Taxodactyl handoff can
  record that the result is partial (spec 02-stages.md §2.1.3, spec
  06a-reports.md §6a.4). Routing it to the minimal bundle would discard a
  real result — CONSTITUTION.md principle 7's second half.
- **`no_barcode` is also a full bundle.** The organelle assembly and its
  annotation are valuable standalone outputs (CONSTITUTION.md principle
  6) even when no panel locus survived validation. Only the *barcode*
  claim is negative, and it is `validation.tsv` — shipped in the bundle —
  that records per-locus why.

Dispatch on a **parsed, explicit allowlist**, never a substring or prefix
test, and **fail closed** on an unrecognised value. This mirrors the
branch already in `main.nf` and the reasoning in spec 02-stages.md §2.1.3
(CONSTITUTION.md rule 18): a status added later must not land on
whichever side its spelling happens to match.

### 3.2 Deriving `no_assembly`

Do not infer this from a missing file — every upstream stage runs and
emits something even when the assembly is empty. The two genuine triggers
are:

- `bin_metadata.json`'s `contigs_selected` is empty, i.e. C3 selected no
  target contig; and/or
- `target.fasta` is present but zero-length.

The second case is specifically the **withheld plastid substitution**
(spec 03-organelles.md §3.6, task 24_plastid_substitution_guard.md §3.2,
carried in tasks/todo.md). A `plant_pt` sample with
`plastid_canonicalisation.substitution_applied: false` and
`substitution_withheld_reason: "no_c3_selection"` emits an **empty**
`target.fasta` alongside a **populated** `plastid_isoforms/` directory.

This is the trap: the sample looks superficially like a success — there
are isoform files on disk, and every downstream stage ran without error —
but C3 found no recognisable plastid, so the canonical path was withheld
deliberately to avoid emitting a confident ~150 kb false positive. It
must route to `no_assembly`, and the isoforms must be shipped as
**diagnostics**, with `substitution_withheld_reason` surfaced in
`metadata.json`. A full bundle full of empty downstream outputs would be
exactly the "output indistinguishable from a confident negative" that
principle 7 forbids.

### 3.3 Reconciling with spec 06a-reports.md §6a.4

Spec §6a.4 lists the run-level status vocabulary as `ok` /
`low_coverage` / `no_recovery` / `fail` / `error`, whereas
CONSTITUTION.md principle 7 names `no_assembly` and `no_barcode`
separately. The constitution wins.

**Decision for this task:** `metadata.json` carries the fine-grained
five-value `sample_status` above. `no_recovery` is a *display grouping*
belonging to C7, not a stored value, and `error` is a C7-only state (no
`metadata.json` exists for a sample that crashed). Record this in the
spec update (§10) so task 45_run_report.md inherits it unambiguously.

---

## 4. The bundle

Per spec 00-overview.md, `outdir/<sample_id>/` contains
`organelle_assembly.fasta`, `organelle_annotation.gff`, `barcodes.fasta`,
`metadata.json`, `report.html`.

- **Full bundle:** all five. `organelle_assembly.fasta` from
  `target.fasta`, `organelle_annotation.gff` from the scored annotation
  GFF (`ANNOTATION_SCORING`'s output, not `ANNOTATE`'s raw one),
  `barcodes.fasta` verbatim from C5.
- **Minimal bundle:** `metadata.json` + `report.html` only. No empty
  placeholder FASTA/GFF — an empty file that looks like a result is worse
  than an absent one (CONSTITUTION.md rule 18).
- **Diagnostics** (`validation.tsv`, `secondaries.tsv`,
  `bin_metadata.json`, `annotation_summary.json`, the graph PNG, the
  organelle map SVG, `plastid_isoforms/`) go to a subdirectory rather
  than the bundle root, so the five specified files stay the obvious
  deliverable. Confirm the subdirectory name against spec 00-overview.md
  and extend that spec's tree if this task adds to it.

### 4.1 The `optional: true` placement bug — fix before anything else

`modules/local/collate.nf` marks `organelle_assembly.fasta`,
`organelle_annotation.gff` and `barcodes.fasta` `optional: true`
**per path, inside the `tuple`**. That placement is a **no-op** in
Nextflow 25.10.2 — confirmed with a minimal standalone repro during task
38 — and the process still errors with "Missing output file(s)" when one
of those paths is not created. `optional: true` must sit at the **tuple**
level, after the last `path(...)` and before `emit:`; see
`modules/local/miniprot_cds.nf`'s `resolved` output for the corrected
form.

This has never fired because the stub `touch`es every file
unconditionally. **It fires the moment the minimal bundle is
implemented** — which is this task. Fix it first, or every soft-failed
sample crashes the process and breaks cross-sample isolation
(CONSTITUTION.md principle 8).

---

## 5. `metadata.json` — the contract

### 5.1 Contents

A single JSON document, `$schema: "wf5/sample-metadata/v1"`, carrying at
minimum:

- **Sample identity** — `sample_id`, `assembly_target`, derived `kingdom`
  and `organelle`, plus every submitter-supplied optional column
  (`sample_info`, `sample_type`, `sample_receipt_date`,
  `storage_location`) threaded from `meta` (spec 00-overview.md).
- **Derived `sample_status`** (§3) with a human-readable `reason`.
- **Coverage** — the whole of `coverage.json` and `sample_status.json`:
  `estimated_cov`, both floors, `coverage_basis`, `subsampled`,
  `fraction`, `seed`, and the sibling-organelle split. Where
  `coverage_basis` is `total_recruited`, that fact must survive into the
  report — it means the sibling split was unavailable and the estimate
  may over-state target depth.
- **Assembly + binning** — from `assembly_info.txt` and
  `bin_metadata.json`: contig count, N50, total bp, `contigs_selected`,
  `target_source`, `sibling_carryover`, `circular`, and the
  `plastid_canonicalisation` block where present.
- **Homology** — top BLAST hit per target contig from the C5-adjacent
  BLAST TSV (`qaccver saccver pident length qcovs evalue bitscore
  stitle`).
- **Barcodes** — per-locus outcome from `validation.tsv`, including the
  `not_found` / `invalid_length` / `identity_below_floor` /
  `internal_stop_codon` reasons for loci that dropped out.
- **Annotation** — the whole `annotation_summary.json` payload, or a
  reference to it (see §5.4).
- **Provenance (CONSTITUTION.md rule 16)** — tool versions, pipeline
  commit, reference-bundle version, and `genetic_code.json` in full: the
  selected table *and*, on `animal_mt`, every candidate considered with
  its score and gene count. Recording only the winner would lose the
  audit trail for a decision made under uncertainty (task 38 §4).

Author the JSON Schema as `assets/sample_metadata.schema.json`, alongside
the existing `assets/samplesheet.schema.json`, and validate against it in
the unit tests (§8). P4a's acceptance criterion in spec 06-phases.md
requires the Taxodactyl handoff bundle to validate against a schema, and
no such schema exists yet.

### 5.2 The reference-bundle version is currently unreachable

CONSTITUTION.md rule 16 and principle 10 both require the reference-
bundle version in `metadata.json`. It cannot currently be obtained.

`nextflow.config` exposes four *independent subdirectory* params —
`organelle_refs`, `protein_panel`, `blast_db`, `annotate_refs` — each
pointing inside `refs/v<version>/`. The `manifest.json` that records the
bundle version sits at the bundle **root**, and no param points there.

Add a `params.refs_manifest` (or a bundle-root param the four subdirs
derive from), materialise it as a value channel in `main.nf`, and stage
it into `COLLATE` as `input: path` — pattern 1 of spec
01-pipeline-flow.md §1a, and required by CONSTITUTION.md rule 13, since a
bare `${params.x}` string cannot be staged on a remote executor. Extend
`validateParams()`'s existing existence check to cover it, and add it to
`conf/integration.config` and `conf/stub.config`.

The manifest's shape is `{$schema, version, generated_at, refseq,
recruit_panels, protein_panel, annotate_refs, artefacts}`. Copy
`version` and `generated_at` into `metadata.json`; do not inline the
whole manifest — that is `run_manifest.json`'s job in task
45_run_report.md.

### 5.3 Decide one gene-name convention

tasks/todo.md (task 39 carry-forward) explicitly defers this decision to
the task that consumes these fields — this one.

`annotation_summary.json`'s `cds_crosscheck.agreed` and `annotator_only`
report **raw** annotator names (`atp6`, `cob_1`, `cox1_0`, `nad2_1` —
MITOS2's own fragment-suffixed strings), while `cds_rescued` reports
**canonical** ones (`ATP8`, `ND3`). Both are defensible alone — raw names
preserve exactly what the annotator said, canonical names are what a
reader looks up — but mixing them in one document means a consumer cannot
join the two lists without re-normalising.

**Recommended:** emit both, as `{"raw": [...], "canonical": [...]}`.
It costs little, loses no information, and satisfies CONSTITUTION.md rule
18 — the auditor reading the report in six months can see both what the
tool said and what it means. Whichever is chosen, apply it uniformly and
delete the todo.md entry.

### 5.4 Open question to settle during implementation

Does `metadata.json` **inline** the full `annotation_summary.json` and
`bin_metadata.json` payloads, or reference them as sibling files in the
diagnostics subdirectory?

Inlining makes `metadata.json` a single self-contained handoff artifact —
better for Taxodactyl and for archival — at the cost of duplication and
size. Referencing keeps it small but means a consumer needs the whole
directory. Recommend **inlining**, on rule 18 grounds: the handoff record
should survive being copied on its own. Record the decision either way.

---

## 6. Container

`conf/containers.config` currently pins `COLLATE` to `python:3.12-slim`
with a `TODO P4` comment. Replace it with the shared report image — the
one C6, C7 and the Jinja renderer all use (spec 02-stages.md §2.2, spec
06a-reports.md §6a.5) — pinned by SHA, never `latest` (CONSTITUTION.md
rule 11).

The image needs Python plus Jinja2 (for task 43) and a JSON-schema
validator. Per CONSTITUTION.md rule 12 it carries **deps only** —
`bin/collate.py` is auto-staged by Nextflow at runtime, not baked in.
Build it in this task even though only the JSON half is used, so task 43
inherits a working image rather than building one mid-task; add its deps
to `scripts/requirements.txt` and let `images.yml` rebuild.

---

## 7. `bin/collate.py` shape

Pseudocode only — this is a brief, not an implementation.

```
parse args: sample meta fields, every input path, refs manifest
read sample_status.json, coverage.json
read bin_metadata.json, annotation_summary.json, validation.tsv,
     genetic_code.json, blast tsv, assembly_info.txt   (each may be
     absent on the failed branch — absence is data, not an error)

sample_status = classify(gate_status, contigs_selected,
                         target_fasta_size, loci_passed)   # §3.1

if bundle_is_full(sample_status):
    copy target.fasta      -> organelle_assembly.fasta
    copy scored annotation -> organelle_annotation.gff
    copy barcodes.fasta    -> barcodes.fasta
copy diagnostics -> <diagnostics subdir>/

metadata = assemble_metadata(...)          # §5.1
validate(metadata, assets/sample_metadata.schema.json)
write metadata.json
touch report.html                          # task 43 fills this in

exit 0 ALWAYS
```

**Always exit 0.** A collation failure must not abort a batch
(CONSTITUTION.md principle 8). A sample whose inputs are unreadable gets
a `metadata.json` recording that fact, not a stack trace that takes its
siblings down. This mirrors C8's existing always-exit-0 contract.

---

## 8. Unit tests — `scripts/tests/test_collate.py`

C6 is custom logic, so CONSTITUTION.md rule 14 applies: **100% branch
coverage**, run via `scripts/pytest.sh`.

Mock at the tool boundary only, per task
27_unit_test_boundary_mocking.md — build real temporary input files and
assert on real output files. Cases to cover:

- Each of the five `sample_status` values dispatches to the right bundle.
- An **unrecognised** gate status fails closed rather than defaulting to
  `ok`.
- `low_coverage` yields a **full** bundle with the warning present in
  `metadata.json` — the specific regression §3.1 warns about.
- `no_barcode` yields a **full** bundle (assembly + annotation shipped,
  `barcodes.fasta` empty or absent per the §4 decision).
- The withheld-substitution `plant_pt` case: empty `target.fasta` +
  populated `plastid_isoforms/` → `no_assembly`, isoforms shipped as
  diagnostics, `substitution_withheld_reason` present.
- Minimal bundle emits **no** empty placeholder FASTA/GFF.
- Emitted `metadata.json` validates against the schema in every branch.
- Provenance round-trip: `animal_mt`'s full `genetic_code.json` candidate
  list survives into `metadata.json`.
- Missing / malformed / zero-byte inputs still exit 0 and still produce
  a schema-valid `metadata.json`.

---

## 9. Integration reconciliation

- **Stub run.** `conf/stub.config` gains the new refs-manifest param;
  `-stub-run` must stay green, including the widened input tuple.
- **`tests/integration/assertions.sh`.** It currently asserts only
  `sample_status.json .status == "ok"`. Add per-sample bundle assertions:
  the five specified files exist for a full bundle, `metadata.json`
  validates against the schema, and the derived `sample_status` matches
  the fixture's expectation.
- **`INT-PLANT-01-mt` is the `low_coverage` fixture.** After task 35 it
  sits at 28.48x — between the 10x hard floor and the 30x warn floor —
  so it is the only real-data exercise of the `low_coverage` full-bundle
  path. Assert that it produces a **full** bundle. Note it is still off
  `assertions.sh`'s `ASSEMBLING_SAMPLES` list (task
  36_plant_mt_integration_coverage.md, shelved); adding the bundle
  assertion does not require adding it there, and this task should not
  absorb task 36.
- **No fixture currently exercises `fail`.** The two-floor gate means no
  integration fixture soft-fails any more. Either construct a minimal
  synthetic low-depth fixture, or cover `fail` in the stub profile and
  record that real-data coverage of the minimal bundle is absent — do
  not leave it silently untested (CONSTITUTION.md rule 18). Recording the
  gap explicitly is acceptable; pretending it is covered is not.
- Keep fixture drift in mind (CONSTITUTION.md rule 19): assert on
  bundle *structure* and status, not on exact coverage figures that
  shift with reference-data rebuilds.

---

## 10. Spec updates

- **spec/02-stages.md** — expand stage 15's row and the C6 entry in §2.2
  with the five-value classification and the full/minimal dispatch.
- **spec/06a-reports.md** — record the §3.3 reconciliation: `no_recovery`
  is a C7 display grouping, not a stored `metadata.json` value.
- **spec/00-overview.md** — extend the per-sample output tree with the
  diagnostics subdirectory.
- **spec/04-reference-data.md** — document the new refs-manifest param
  and why `COLLATE` needs the bundle root.
- **tasks/todo.md** — delete the three `COLLATE` carry-forward items and
  the task 39 gene-name item, now resolved here. Leave the `RUN_REPORT`
  and per-sample-report items for tasks 45 and 43.

Per the task-authoring convention: when editing spec files, reference
tasks as plain text (`task 42_collate_bundle_metadata.md`), never as
markdown links, and never with line anchors.

---

## 11. Acceptance criteria

1. `bin/collate.py` exists, is flake8-clean at 79 columns, and has 100%
   branch coverage under `scripts/pytest.sh`.
2. All five `sample_status` values dispatch correctly; unrecognised
   statuses fail closed.
3. `low_coverage` and `no_barcode` both produce full bundles; `fail` and
   `no_assembly` produce minimal bundles with no empty placeholder
   FASTA/GFF.
4. The withheld-substitution `plant_pt` case routes to `no_assembly`
   with its isoforms shipped as diagnostics.
5. `optional: true` sits at tuple level in `modules/local/collate.nf`,
   and a minimal-bundle sample completes without a "Missing output
   file(s)" error.
6. `assets/sample_metadata.schema.json` exists; every emitted
   `metadata.json` validates against it.
7. `metadata.json` records the reference-bundle version, the pipeline
   commit, tool versions, and `animal_mt`'s full genetic-code candidate
   list.
8. `bin_metadata.json`, `plastid_isoforms/`, the Bandage PNG and
   `assembly_info.txt` all reach `COLLATE` (the §2.2 defects are fixed).
9. `COLLATE` exits 0 on every input, including malformed ones.
10. `-stub-run` green; `assertions.sh` extended and passing on the
    nightly integration profile.

---

## 12. Out of scope

- **`report.html` content** — task 43_per_sample_report.md. This task
  emits an empty placeholder.
- **`run_manifest.json` / `run-report.html`** — task 45_run_report.md.
- **The organelle map SVG** — task 44_organelle_map.md. `COLLATE` copies
  whatever it is given, stub or real.
- **Adding `INT-PLANT-01-mt` to `ASSEMBLING_SAMPLES`** — task
  36_plant_mt_integration_coverage.md.
- **Sourcing a deeper `plant_mt` fixture** — tasks/todo.md, unchanged.
- **Trimming over-extended miniprot CDS spans** — tasks/todo.md (task 39
  carry-forward, item b), unchanged.

---

## 13. Outcomes

To be completed when the task is executed — deviations from this brief,
measurements taken, and any new carry-forward items for tasks/todo.md.
