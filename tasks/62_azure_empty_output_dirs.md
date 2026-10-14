# Task 62 — Azure: directory outputs must never be empty

**Phase:** P5 — production deployment (spec/06-phases.md).

**Status:** open. **No preconditions.** Independent of task
completed/61_intake_integrity_and_failure_visibility.md. Without task 61, a
sample hit by this bug vanishes silently. With it, the sample shows as
`error`. This task stops the failure from happening at all.

## Overview

Azure Blob Storage has no real directories. A "directory" is only a
shared prefix of the files inside it, so an empty directory does not
exist once it is uploaded. On Azure Batch, Nextflow's work directory
is in blob storage. A process that declares a directory output and
leaves it empty therefore fails there with "missing output", even
though the same run passes locally.

Task 57's Azure stub run (§5.2) hit this. `TRANSLATE_ANNOTATION_CDS`
declares `path("queries")`, and its stub only runs `mkdir -p queries`.
The real script branch can leave the directory empty too:
`bin/translate_annotation_cds.py` writes one `<gene>.faa` per gene
with a CDS feature, so a sample with zero CDS features writes nothing.
That case is plausible for an assembly that passed binning but is
short or divergent. With `errorStrategy = 'ignore'` the task is
dropped. `SCORE_ANNOTATION_BLASTP`, `JOIN_ANNOTATION_SCORES`,
`ORGANELLE_MAP` and `COLLATE` then never run for that sample, which is
exactly the silent drop task 61 describes.

The fix is a small placeholder file in every directory output that can
be empty, plus a rule in the spec so new processes do the same.

## 1. Constraints

- **Rule 11/13 (containers, stageable paths):** no new tools. A
  placeholder is a `touch`.
- **Behaviour on local and HPC executors is unchanged.** A consumer
  must ignore the placeholder: globs must not match it, and nothing may
  count it as content.
- **Placeholder name: `.keep`.** Use one name across the pipeline.
  It's hidden, so `ls`-style listings and `*.faa` globs skip it.

## 2. Scope

### 2.1 `TRANSLATE_ANNOTATION_CDS`

`touch queries/.keep` in both `script:` and `stub:`. Write it in the
`.nf` block rather than in the Python, because the stub never runs the
Python.

Check the consumer: `SCORE_ANNOTATION_BLASTP` globs
`${queries_dir}/*.faa` under `nullglob`. The placeholder is not
matched, and zero `.faa` files give an empty `blast.tsv`. Then follow
a zero-CDS sample through `JOIN_ANNOTATION_SCORES`, `ORGANELLE_MAP` and
`COLLATE`, and confirm it ends as a normal bundle (no scored features)
with no crash. Fix any stage that assumes at least one CDS. If one
needs a `bin/` change, add a unit test for the zero-feature input.

### 2.2 Audit every directory output and input

List every `path("<dir>")` or `path("<dir>/")` output in
`modules/local/*.nf`. For each, record in Outcomes whether it can be
empty on a real run or in the stub, and what was done. Known
candidates:

| Process | Output | Notes to check |
|---|---|---|
| `TRANSLATE_ANNOTATION_CDS` | `queries` | Fixed in §2.1. |
| `NANOPLOT_RAW` / `NANOPLOT_CLEAN` | `nanoplot_*/` | Always populated on success. The stub touches files. Confirm. |
| `BIN_TARGET` | `plastid_isoforms/` (`optional: true`) | If the script can create the directory empty, an absent output is still fine for Nextflow but changes meaning. Make sure "absent" only ever means "not the canonical branch". |
| `COLLATE` | `diagnostics` (`optional: true`) | Same check. |
| `MINIPROT_CDS` | `candidates/cds.*.gff` | A file glob, not a directory. The GFFs can be empty files, which blob storage stores fine. No change expected. |

Also check directory **inputs** staged from blob storage. The
reference bundle subdirectories are uploaded by
`scripts/publish_refs_extracted.sh`. An empty subdirectory in a bundle
would vanish the same way. Check `v2026.09_1` and record the result.
Task 65_plastid_noncds_barrnap.md already stages a placeholder file
rather than an empty directory, so it follows this rule.

### 2.3 A lint check

Add a check to the `lint` job in `.github/workflows/tests.yml`, in the
same style as the container-coverage check. For every module that
declares a directory output (a `path("x")` or `path("x/")` output with
no file extension and no glob), at least one of `script:`/`stub:` must
write into it. If both write `x/.keep`, the check passes. Otherwise the
module is listed on an explicit allow-list with a one-line reason (for
example "NanoPlot always writes its report"). Keep it simple. The goal
is that a new empty-directory output is a deliberate decision, not
something found on Azure.

## 3. Testing

- **Unit:** only if §2.1 needs a `bin/` change for zero CDS. If so,
  add a zero-feature test there.
  `scripts/tests/test_translate_annotation_cds.py` should already
  cover "no CDS → no `.faa` written". Add the test if it doesn't.
- **Stub (local):** passes, and `queries/.keep` exists in the
  `TRANSLATE_ANNOTATION_CDS` work directory.
- **Lint:** the new check passes on the tree and fails on a scratch
  module with an unmarked empty directory output. Revert the scratch
  change.
- **Integration (local):** passes unchanged. The placeholder appears
  under `annotation/scoring/queries/` when `publish_intermediates` is
  on, and nothing else changes.
- **Azure (task 57 launcher, `--stub` on the `wf5_test` pool):**
  completes with **zero ignored tasks** (task 57 recorded one). Check
  the end-of-run summary or `nextflow log`. The user approved this
  pool run on 2026-10-14.

## 4. Spec updates

- spec/01-pipeline-flow.md §1a: a new engineering-constraint bullet,
  "Directory outputs are never empty". Explain the blob-storage reason,
  the `.keep` convention, that consumers must not match it, and the
  lint check.
- spec/04a-azure-blob-storage.md: one line pointing at that rule.

## 5. Out of scope

- The other Azure items from task 57 (`ps` in images, which is task
  67_procps_in_images.md, and resume ergonomics, which stays in
  tasks/todo.md).

## 6. Checklist

- [x] §2.1 `queries/.keep` in script and stub; zero-CDS path traced to
      COLLATE
- [x] §2.2 audit table recorded, including reference bundle subdirs
- [x] §2.3 lint check added; proven to fail on a scratch module
- [ ] §3 stub and integration green; Azure stub with zero ignored tasks
- [x] §4 spec updated
- [ ] Outcomes section written

## Outcomes

Status: **partly done.** The integration run and the Azure stub run are
not done, so the task stays in tasks/.

- §2.1: `touch queries/.keep` added to script and stub of
  `TRANSLATE_ANNOTATION_CDS`. `SCORE_ANNOTATION_BLASTP` uses a nullglob
  over `*.faa`, so it skips `.keep` and writes an empty `blast.tsv`.
  `annotation_scores.py` and its tests already cover empty/missing
  blast input (explicit null scores), and
  `test_translate_annotation_cds.py::test_no_features_no_queries_written`
  already covers zero CDS. No `bin/` change needed. I did not run a
  zero-CDS sample through ORGANELLE_MAP/COLLATE end to end; that is
  covered only by reading the code.
- §2.2 audit of directory outputs in `modules/local/*.nf`:
  - `queries`: could be empty (zero CDS); fixed.
  - `nanoplot_raw/`, `nanoplot_clean/`: always populated on success; the
    stubs touch files. Allow-listed.
  - `plastid_isoforms/` (BIN_TARGET, optional): `bin_target.py` removes
    the dir off the canonical branch, so absent means "not canonical";
    C4 writes files on it. The stub makes none. Allow-listed.
  - `diagnostics` (COLLATE, optional): `collate.py` creates it only when
    there is content. Allow-listed.
  - `MINIPROT_CDS` `candidates/cds.*.gff`: file glob, no change.
  - Reference bundle `refs/v2026.09_1`: `find -type d -empty` finds no
    empty directories.
- §2.3: "Directory-output placeholder check" added to the `lint` job in
  `.github/workflows/tests.yml`. It passes on the tree (5 outputs) and
  fails with the `queries` fix reverted.
- §3: local stub run passes (Nextflow 25.04.0) and `queries/.keep` is in
  the work dirs. The local default Nextflow (26.04.6) fails to parse
  `nextflow.config` (`if` in the `plugins` block), which is task 60
  territory and unrelated. Not run: local integration, and the Azure
  `--stub` run on `wf5_test` (zero ignored tasks).
- §4: spec/01-pipeline-flow.md §1a and spec/04a-azure-blob-storage.md
  updated.
