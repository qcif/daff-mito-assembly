# Task 60 — CI: run the `bin/*.py` unit tests and gate on coverage

**Phase:** P5 — production deployment (spec/06-phases.md).

**Status:** open. **No preconditions.** Do this first. It is small,
and every task after it changes `bin/*.py`, so the tests should gate
those changes.

## Overview

All custom pipeline logic lives in `bin/*.py`: samplesheet validation,
the coverage gate, binning, barcode validation, annotation merge, and
the reports. CONSTITUTION rule 14 requires each of these scripts to be
unit-tested at 100 % branch coverage, and spec §5b says CI runs those
tests on every push.

CI does not run them. The `python-tests` job in
`.github/workflows/tests.yml` runs `flake8` and then an "Unit tests
(placeholder)" step. That step looks for `tests/unit/test_*.py`, a
directory task 27 removed, finds nothing, prints "No unit tests yet",
and passes. The real suite in `scripts/tests/` has 778 tests. It runs
only when someone runs `scripts/pytest.sh` by hand.

On 2026-10-14 the suite passes with 100 % branch coverage on every
`bin/` module. Wiring it into CI now costs nothing in fixes and
protects that state while the plastid annotation and deployment
changes land.

## 1. Constraints

- **Same environment as local.** spec §5b item 4: the tests run
  inside the `neoformit/daff-wf5-scripts:test` image via
  `scripts/pytest.sh`, not host Python. `test_bin_target.py` needs
  `mappy` pinned to the `BIN_TARGET` container's minimap2 version,
  and the image is where that pin lives (`scripts/requirements-test.txt`).
- **Runnable locally.** spec §5b "Local re-run": no CI-only tooling.
  CI calls the same `scripts/pytest.sh` a developer does.
- **Rule 14:** 100 % branch coverage is the target for `bin/` modules.
  CI enforces it. Lowering the bar needs a constitution amendment, not
  a CI tweak.

## 2. Scope

### 2.1 Build the test image in CI

The `:test` image is not published (`images.yml` pushes only the
runtime image). Build it in the `python-tests` job from
`scripts/Dockerfile` with `--build-arg TEST=1`, tag it
`neoformit/daff-wf5-scripts:test`, and load it into the local Docker
daemon so `scripts/pytest.sh` finds it and skips its own build.

Use the GitHub Actions layer cache (`docker/build-push-action` with
`cache-from/to: type=gha`, `push: false`, `load: true`). The `mappy`
compile is the slow layer, and it only changes when
`scripts/requirements*.txt` or `scripts/Dockerfile` change.

### 2.2 Run the suite and gate on coverage

Replace the placeholder step with `bash scripts/pytest.sh`. Running
the whole of `scripts/tests/` is the script's default.

Add a coverage floor that applies in CI but not to ad-hoc local runs
of a single test file, which would otherwise always fail. Sketch:

```
# scripts/pytest.sh
FAIL_UNDER=${COVERAGE_FAIL_UNDER:-}        # unset locally → report only
coverage report -m --include=... ${FAIL_UNDER:+--fail-under=$FAIL_UNDER}

# tests.yml
env: COVERAGE_FAIL_UNDER: 100
```

Document the variable in the script's header comment.

`scripts/pytest.sh` runs `docker run --user "$(id -u):$(id -g)"` and
mounts `/tmp`. Check that both work on the hosted runner.

### 2.3 Tidy the job

- Rename the step and job so the check name says what it does (for
  example "Python unit tests (pytest, 100 % branch) + flake8"). Keep
  the job id stable.
- Leave the `flake8` step as it is.

### 2.4 Branch protection: leave as-is (decided 2026-10-14)

spec §5b's "`main` requires `tests.yml` green before merge" is the
policy and stays. Pull requests are rare, and most commits go straight
to `main`, where a red `tests.yml` run on the pushed commit is the
signal to fix. **Do not change repository settings.**

## 3. Testing

No unit tests are needed. This task changes CI and one shell script.
Verify it end to end:

1. **Locally:** `bash scripts/pytest.sh` behaves as before (report
   only). `COVERAGE_FAIL_UNDER=100 bash scripts/pytest.sh` passes on
   the current tree. With one branch-covering test temporarily skipped,
   it fails on coverage. That shows the floor works. Revert the skip.
2. **In CI, on a throwaway branch** (the user approved pushing one on
   2026-10-14):
   - The branch as committed gives a green job with 778+ tests reported
     in the log.
   - A commit adding `assert False` to one test gives a red job.
   - A commit adding an untested branch to a `bin/` module gives a red
     job on coverage.
   Delete the branch afterwards. Record the job's wall time with a cold
   cache and with a warm one in Outcomes.
3. The stub and lint jobs are unchanged and still green.

## 4. Spec updates

- spec/05-test-data.md §5b `tests.yml` item 4: state that CI builds
  the `:test` image with the GHA layer cache and that it enforces
  `--fail-under=100` through `COVERAGE_FAIL_UNDER`.
- spec/05-test-data.md §5b "Local re-run": the list names
  `python -m unittest` and `nf-test`, neither of which the project
  uses. Replace them with `scripts/pytest.sh` and `-stub-run`.

## 5. Out of scope

- A coverage gate for `scripts/refdata/*.py`. Their tests
  (`test_build_proteins.py`, `test_mask_panel.py`) run in the suite,
  but coverage is measured only for `bin/`. Task
  64_protein_panel_label_check.md keeps `build_proteins.py` coverage
  where it is. Widening the gate is a separate decision.
- Pulling the inline lint heredocs out into scripts (tasks/todo.md,
  Test + CI backlog).
- Pinning `flake8`/`flake8-bugbear` versions.

## 6. Checklist

- [ ] §2.1 test image built in CI with layer cache
- [ ] §2.2 placeholder replaced; `COVERAGE_FAIL_UNDER` added and
      documented
- [ ] §2.3 job and step renamed, job id stable
- [ ] §3 local and throwaway-branch checks done; branch deleted; wall
      times recorded
- [ ] §4 spec updated
- [ ] Outcomes section written
