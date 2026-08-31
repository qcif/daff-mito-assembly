# Task 41 — `ANNOTATE`: a status for "the annotator ran and returned nothing"

**Phase:** P4a (a correction to the stage 13a merge built in task 31).
**Prerequisite:** none. Independent of tasks 38–40; touches the status
vocabulary and the summary writer, not the merge or the scoring.

**Goal:** Make a failed specialist annotator a reported state rather
than a silent one, so an annotation missing every tRNA and rRNA can
never again be published as `status: "ok"`.

---

## 1. Overview

Stage 13a merges miniprot CDS calls with a specialist annotator's
non-CDS features. `annotate_summary.json` carries four statuses
([spec §2 C8](../spec/02-stages.md#22-custom-logic-components)):

| Condition | `status` |
|---|---|
| `cds.gff` empty | `no_assembly` |
| CDS merged, no non-CDS annotator for this target | `ok_cds_only` |
| CDS merged + non-CDS features present | `ok` |
| Nothing parsed from either source | `no_features` |

There is a fifth condition, and it is the one that actually occurred:
**a non-CDS annotator was configured, was invoked, and produced
nothing.** `bin/annotate_summary.py` `run()` reaches its `else` branch
whenever `mitos_dir is not None`, regardless of whether anything was
parsed out of it, and writes `status: "ok"`.

### 1.1 The incident

The nightly integration test failed on every run from 2026-08-07
(`c9bad169`, the commit that introduced the MITOS2 branch) to
2026-08-16 — nine consecutive nights — on:

```
FAIL: INT-ANIMAL-01 tRNA count 0 below floor 15
FAIL: INT-ANIMAL-01 rRNA count 0 below floor 1
```

The cause was environmental: `scripts/refdata/build_annotate.sh` was run
on a workstation on 2026-08-07, but `refs-v2026.09.tar.gz` in blob
storage was never rebuilt. Its 182 entries are `recruit/`, `validate/`,
`proteins/` and `manifest.json` — no `annotate/`. In CI,
`refs/v2026.09/annotate` therefore did not exist, Nextflow staged the
directory input as a **dangling symlink** (`ln -s` succeeds on a
non-existent target), and MITOS2 exited 2 with
`no such directory annotate/mitos/refseq89m`.

Two guards have since closed the specific hole:

- `main.nf` `validateParams()` now existence-checks every refdata param
  before any process runs.
- `modules/local/annotate.nf` records the annotator's exit code and
  output to a published `mitos.log` instead of discarding it with
  `|| true`, and `.github/workflows/integration.yml` collects the work
  directory on failure.

Neither addresses what let the defect survive nine nights. The run
published, for a metazoan mitogenome, an annotation asserting 9 genes
and **zero** tRNAs — a biological impossibility, since 22 of the
canonical 37 metazoan mitochondrial genes are tRNAs
([spec §2 stage 13a](../spec/02-stages.md#2-stage-detail)) — and
labelled it `ok`. Only `assertions.sh`'s numeric floors caught it. A
production run has no such floors, and under
[principle 18](../CONSTITUTION.md) a summary that reports a tool
failure as a clean result is the failure mode to close.

**This is a status-vocabulary defect, not a MITOS2 defect.** The same
hole is open for every future annotator: task 32 (`plant_pt`) and task
33 (`plant_mt`) each add a non-CDS tool and inherit it.

**Exit criteria:**

- An invoked-but-unproductive non-CDS annotator produces a distinct,
  self-describing status, never `ok`.
- The annotator's exit code reaches `annotation_summary.json`.
- `ok` means what the table says: non-CDS features are present.

---

## 2. The status to add

Add `annotator_failed` as a fifth value:

| Condition | `status` |
|---|---|
| `cds.gff` empty (empty `target.fasta` upstream) | `no_assembly` |
| CDS merged, no non-CDS annotator for this target | `ok_cds_only` |
| **CDS merged, annotator configured but yielded no features** | **`annotator_failed`** |
| CDS merged + non-CDS features present | `ok` |
| Nothing parsed from either source | `no_features` |

Three properties matter:

1. **It is distinguishable from `ok_cds_only`.** Both ship a CDS-only
   GFF, but they mean opposite things: `ok_cds_only` is the honest,
   provisional state of a target with no annotator yet
   ([spec §8 item 3](../spec/07-open-questions.md#8-remaining-open-questions));
   `annotator_failed` is a tool that was expected to deliver and did
   not. Collapsing them would hide this defect just as effectively as
   `ok` did.
2. **It is not a sample failure.** Annotation is supplementary and
   never gates a sample ([spec §2 stage 13a](../spec/02-stages.md#2-stage-detail)).
   `COLLATE`'s negative-clarity classification
   ([principle 7](../CONSTITUTION.md)) is untouched — this status
   describes the annotation, not the specimen.
3. **C8 still exits 0** ([principle 8](../CONSTITUTION.md)). The status
   is the signal; the exit code stays 0.

### 2.1 Zero non-CDS features is the trigger, not the exit code

The test is on the *output*, not on how the tool terminated. MITOS2
exits 0 on some no-output paths, and a tool that dies after writing a
partial `result.gff` should be judged on what it wrote. So:

```
annotator configured AND parsed non-CDS feature count == 0
    -> annotator_failed
```

The exit code is recorded as evidence, not used as the test — see §3.

### 2.2 `reason` must name the cause

`reason` is currently `null` for `ok`. For `annotator_failed` it
carries the annotator, its exit code where known, and the directory
that was searched, e.g.:

```
"mitos2 exited 2 and wrote no result.gff under mitos_out/"
```

Enough for a report reader to act without opening `mitos.log`.

---

## 3. Plumbing the exit code

`modules/local/annotate.nf` already captures `MITOS_EXIT`. Pass it to
C8 as `--annotator-exit <n>`, defaulting to `null` when no annotator
ran, and surface it in the summary alongside the existing tool
versions:

```json
"non_cds_source": "mitos2",
"annotator_exit_code": 2,
```

`null` when no annotator was configured. This is provenance
([principle 16](../CONSTITUTION.md)) — the number that would have
saved nine nights of debugging belongs in the machine-readable output,
not only in a log the CI discards on success.

---

## 4. Bundle completeness — the upstream half

The dangling symlink was reachable because nothing checks that a
fetched bundle contains what its `manifest.json` claims.
`scripts/fetch_refs.sh` verifies the *tarball's* SHA256 and stops
there, so a bundle built before `annotate/` existed passes cleanly and
fails 12 minutes later inside a container.

[Principle 10](../CONSTITUTION.md) requires the bundle be
provenance-tracked; verifying only the archive checksum satisfies the
letter and not the intent. Add a post-unpack check that every
top-level component named in `manifest.json` is present on disk, and
fail `fetch_refs.sh` with the missing names if not.

This is in scope here because it is the same defect one layer up:
a missing input reported as a successful fetch.

---

## 5. Unit tests — `scripts/tests/test_annotate_summary.py`

Extending task 31 §5's suite; 100 % branch coverage maintained.

1. **`annotator_failed`** — `--mitos-dir` supplied, directory contains
   no `result.gff` → status `annotator_failed`, `reason` names the
   tool, CDS features still emitted in full, exit 0.
2. **`annotator_failed` with a parseable but empty `result.gff`** —
   file present, no tRNA/rRNA rows → same status. Distinguishes "no
   features" from "no file".
3. **Not confused with `ok_cds_only`** — no `--mitos-dir` at all still
   gives `ok_cds_only` with `non_cds_source: null`.
4. **Exit code recorded** — `--annotator-exit 2` appears as
   `annotator_exit_code: 2`; omitted → `null`.
5. **`ok` unchanged** — the task 31 happy-path fixture still yields
   `ok`; this task must not reclassify a working annotation.
6. **CDS-only genome is not `annotator_failed` by accident** — an
   annotator that legitimately returns only rRNA (no tRNA) is `ok`.
   The trigger is *zero* non-CDS features, not a missing type.

---

## 6. Integration reconciliation

`tests/integration/assertions.sh` currently asserts
`INT-ANIMAL-01 annotation status='ok'` and then separately floors the
tRNA/rRNA counts. Once this task lands the status assertion alone
catches the regression, but keep both: the floors test biology
(a metazoan mitogenome has ~22 tRNAs), the status tests plumbing.

Add the negative assertion that the status is never
`annotator_failed` on either fixture, with the published `mitos.log`
echoed on failure — already wired in `integration.yml`.

**This task does not make the nightly green.** That needs
`refs-v2026.09.tar.gz` rebuilt with `annotate/mitos/refseq89m` and
`scripts/refs-v2026.09.sha256` updated. This task ensures that if it
regresses, the failure names itself.

---

## 7. Spec updates

- [spec §2 stage 13a](../spec/02-stages.md#2-stage-detail) and the C8
  row in [§3](../spec/02-stages.md#22-custom-logic-components) — "Four statuses"
  becomes five; add `annotator_failed` with its meaning.
- [spec §6a](../spec/06a-reports.md) — the per-sample report must
  render `annotator_failed` distinctly from `ok_cds_only`; a reader
  must not read a tool failure as a provisional-by-design state.
  Pairs with the existing task 31 carry-forward in
  [todo.md](todo.md) about surfacing `cds_crosscheck` disagreements.
- [spec §4](../spec/04-reference-data.md) — record the bundle
  completeness check from §4.

---

## 8. Acceptance criteria

1. A configured annotator yielding zero non-CDS features produces
   `status: "annotator_failed"` with a `reason` naming the tool and
   exit code; C8 exits 0.
2. `annotator_exit_code` is present in `annotation_summary.json`,
   `null` when no annotator was configured.
3. `ok_cds_only`, `ok`, `no_assembly` and `no_features` are unchanged
   on their task 31 fixtures.
4. `fetch_refs.sh` fails with the missing component names when the
   unpacked bundle does not match its `manifest.json`.
5. Unit tests in §5 pass at 100 % branch coverage; `flake8` clean.
6. Spec §2, §4 and §6a updated as in §7.

---

## 9. Out of scope

- **Retrying or repairing a failed annotator.** The status reports; it
  does not recover.
- **Gating the sample on annotation.** Annotation stays supplementary
  and never influences `COLLATE`'s classification.
- **Rebuilding or republishing the refdata bundle.** An operational
  action, tracked separately.
- **The `ok_cds_only` plant arms.** Tasks 32 and 33 own those; they
  inherit this status vocabulary when their annotators land.

---

## 10. Outcomes

Implemented as specified, no deviations from the design. Notes:

- `bin/annotate_summary.py`: added `STATUS_ANNOTATOR_FAILED`, an
  `annotator_exit` parameter threaded through `run()`/`main()`
  (`--annotator-exit`, `type=int`, `None` when omitted), and
  `annotator_exit_code` in every summary branch including
  `no_assembly`. The new status sits in the existing `elif` chain
  between `ok_cds_only` and `ok`, keyed on `not non_cds` (zero
  non-CDS features), so it can never fire when no annotator was
  configured and never fires on a legitimate rRNA-only annotation.
  `reason` distinguishes "no result.gff found under `<dir>/`" from
  "result.gff present but no tRNA/rRNA features" using
  `find_mitos_result_gffs()`, and names the annotator + exit code
  (falling back to "exit code unknown" if none was supplied).
- `modules/local/annotate.nf`: `NONCDS` now appends
  `--annotator-exit ${MITOS_EXIT}` inside the branch that actually
  runs MITOS2; the no-annotator branch passes nothing, so the script
  sees `None` as designed.
- `scripts/fetch_refs.sh`: after unpacking, derives the set of
  top-level bundle components from `manifest.json`'s `artefacts` keys
  (first path segment of each — `manifest.json` has no separate
  literal "components" list) and fails listing any missing directory
  names. Manually replayed the task's own incident (an `annotate/`-less
  bundle) against this logic — it correctly reports `missing: annotate`.
  Verified with `shellcheck` (clean) rather than executing a real
  Azure fetch.
- `tests/integration/assertions.sh`: added an explicit
  `annotator_failed` tripwire on both `INT-ANIMAL-01` (real MITOS2
  annotator; echoes `mitos.log` on failure, per §6) and
  `INT-PLANT-01-pt` (no annotator exists yet, so this cannot fire
  today — kept as a forward guard for task 32).
- Unit tests: added `TestAnnotatorFailed` (6 cases per §5) plus one
  CLI-arg-parsing case in `TestExitCodeZero`. Full suite:
  `bash scripts/pytest.sh` — 234 passed, 100% branch coverage across
  all of `bin/*.py`, including `annotate_summary.py`.
  `flake8 bin/annotate_summary.py scripts/tests/test_annotate_summary.py`
  clean.
- Did not run a full `nextflow run . -profile integration` — the
  refdata/MITOS2 path was validated at the unit level (task brief's
  own §6 note: this task cannot make the nightly green on its own,
  only ensure a regression names itself). No new tool/container was
  introduced, so the brief's "confirm the tool works" light-touch
  check didn't apply here.
- Spec updated: §2 stage-detail C8 row (five statuses,
  `annotator_exit_code`), a new §4.5 "Fetching and verifying a
  bundle" section, and the annotated-genome-map row in §6a
  distinguishing `annotator_failed` from `ok_cds_only` in the report.
