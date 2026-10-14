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

- [ ] §2.1 `queries/.keep` in script and stub; zero-CDS path traced to
      COLLATE
- [ ] §2.2 audit table recorded, including reference bundle subdirs
- [ ] §2.3 lint check added; proven to fail on a scratch module
- [ ] §3 stub and integration green; Azure stub with zero ignored tasks
- [ ] §4 spec updated
- [ ] Outcomes section written
