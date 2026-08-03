# Task 27 — Fake the external-tool boundary; bring C1 + C2 under coverage

**Phase:** cross-cutting test infrastructure. No pipeline stage, module,
or `bin/*.py` behaviour changes — this task restores the
[spec §5a](../../spec/05-test-data.md#5a-tests) pytest surface for C1 and
C2, which is currently both **failing** and **unmeasured**.

**Goal:** Make `scripts/pytest.sh` green across the whole suite, and
bring [`bin/coverage_gate.py`](../../bin/coverage_gate.py) (C2) and
[`bin/parse_samplesheet.py`](../../bin/parse_samplesheet.py) (C1) to the
100 % branch coverage [rule 14](../../CONSTITUTION.md) requires of every
custom-logic component.

**Prerequisite:** none. Independent of
[task 24](24_plastid_substitution_guard.md) and
[task 25](25_coverage_gate_carryover.md), though it touches the same
component as 25 — sequence 27 before 25 so that task inherits a green,
measured C2 to work against.

**Spec basis:**
[§5a](../../spec/05-test-data.md#5a-tests) (the pytest surface, and the
100 % branch-coverage rule of thumb for C1–C7),
[§5b item 4](../../spec/05-test-data.md#5b-ci) (where the Python tests run
— see §5 below, it currently contradicts two other documents),
[§2.2](../../spec/02-stages.md#22-custom-logic-components) (C2's declared
tool contract: "Python stdlib only + `seqkit`/`seqtk` in `PATH`"),
[rule 14](../../CONSTITUTION.md), [rule 15](../../CONSTITUTION.md) (three test
surfaces, three cadences — this task must not create a fourth).

## 1. Two defects, one cause

### 1.1 The suite is red

`scripts/tests/test_coverage_gate.py` fails all 4 tests inside the
`neoformit/daff-wf5-scripts:test` image.
[`bin/coverage_gate.py`](../../bin/coverage_gate.py) shells out to `seqkit`
(unconditionally, via `total_bases()`) and to `bash` + `seqtk` (the
subsample branch only). Neither binary is in that image: they live in
the production `COVERAGE_GATE` container, which is built out-of-band via
`mulled-build` and is not derived from any Dockerfile in this repo.

Logged in [`tasks/todo.md`](../todo.md) during
[task 23](23_bin_target_recalibration.md), which observed but
did not fix it.

### 1.2 The suite does not measure what it claims to

Both `test_coverage_gate.py` and `test_parse_samplesheet.py` invoke
their script as a **child process**:

```
subprocess.run([sys.executable, str(SCRIPT), ...])
```

[`scripts/pytest.sh`](../../scripts/pytest.sh) runs
`coverage run --branch -m pytest`, which traces **only the parent
process**. There is no `conftest.py`, no `.coveragerc`, no
`pytest.ini`/`pyproject.toml`, and no `COVERAGE_PROCESS_START` anywhere
in the repo — so nothing enables subprocess tracing.

Consequence: `bin/coverage_gate.py` and `bin/parse_samplesheet.py`
contribute **zero measured branch coverage** and do not appear in the
`--include='*/bin/*.py'` report at all. C1 and C2 are the two components
[rule 14](../../CONSTITUTION.md) most explicitly covers, and both are at
0 % measured. The other two test modules
(`test_bin_target.py`, `test_plastid_canonicalise.py`) already load
their module in-process and are measured correctly.

**Verify this before changing anything** — it is the premise of the
whole task. Run `bash scripts/pytest.sh scripts/tests/test_plastid_canonicalise.py`
and confirm the report table lists only `plastid_canonicalise.py`.

### 1.3 Why these are one problem

Fixing 1.1 by faking the tool boundary *in-process* fixes 1.2 for free,
because the module is then imported rather than spawned. Fixing 1.1 any
other way (see §6.1) leaves 1.2 untouched.

## 2. Changes to `scripts/tests/test_coverage_gate.py`

Full rewrite of the test module. **[`bin/coverage_gate.py`](../../bin/coverage_gate.py)
is not modified** (except the §4 pragma) — this is a test-side change.

### 2.1 Adopt the existing in-process idiom

The repo already has this pattern in two modules; reuse it rather than
inventing a third. Module load as in
[`scripts/tests/test_bin_target.py`](../../scripts/tests/test_bin_target.py),
CLI entry as in
[`scripts/tests/test_plastid_canonicalise.py`](../../scripts/tests/test_plastid_canonicalise.py):

```
# pseudocode
spec = importlib.util.spec_from_file_location('coverage_gate', CG_PATH)
cg   = importlib.util.module_from_spec(spec)
sys.modules['coverage_gate'] = cg     # lets patch() address it by name
spec.loader.exec_module(cg)

with chdir(tmp), patch('coverage_gate.subprocess.run', side_effect=fake):
    with patch('sys.argv', argv):
        rc = cg.main()
```

Two notes to carry as comments in the file:

- `coverage_gate.py` does `import subprocess`, **not**
  `from subprocess import run`, so `patch('coverage_gate.subprocess.run')`
  patches `run` on the shared module object for the duration of the
  `with` block. Harmless here (single-threaded, nothing else in the
  module shells out) but state it, because it is not obvious.
- `coverage_gate.py` writes `sample_status.json` and `coverage.json` to
  **CWD**. `contextlib.chdir` (3.11+; the image is 3.12) replaces what
  `cwd=tmp` gave the subprocess form.

### 2.2 Fakes must be data-driven, not canned

This is the part that decides whether the rewritten tests still assert
anything. A hard-coded stats table would make the coverage arithmetic a
tautology.

```
# pseudocode
def _stats_table(path):
    """Render a real `seqkit stats -T` table from the real gzip."""
    # gunzip path, count records, sum read lengths,
    # return header line + one data row, tab-separated

def _fake_seqtk(shell_cmd):
    """Emulate `seqtk sample -s S IN F | gzip > OUT`."""
    # regex-parse the command string; Bernoulli-select records with
    # random.Random(seed) — seqtk's real semantics — and write a real gzip

def _fake_run(cmd, **kwargs):
    if cmd[0] == 'seqkit':
        return CompletedProcess(cmd, 0, stdout=_stats_table(Path(cmd[3])))
    if cmd[0] == 'bash':
        _fake_seqtk(cmd[2])
        return CompletedProcess(cmd, 0)
    raise AssertionError(f'unexpected command: {cmd}')
```

Three constraints that are easy to get wrong:

1. **The empty-file case must emit a data row with `sum_len=0`**, as
   real `seqkit stats -T` does. `total_bases()` indexes `out[1]`
   unconditionally; a header-only fake would `IndexError` and make
   `test_empty_input` vacuous rather than passing for the right reason.
2. **Bernoulli selection, not an exact count.** That is what `seqtk
   sample` actually does, it is deterministic given the seed, and at
   n = 100 000 / p ≈ 0.51 the standard deviation is ~0.16 % — two orders
   of magnitude inside the existing ±10 % assertion, so no flake risk.
3. **The fall-through must raise.** If a future change adds a third tool
   call to C2, the fake must fail loudly rather than silently returning a
   `MagicMock` that satisfies every assertion.

### 2.3 Existing cases

The four existing cases (`passthrough`, `low_coverage`, `subsample`,
`empty_input`) keep their bodies, arithmetic comments and assertions
**verbatim** — only the driver call changes. Add `rc == 0` to each.

`test_subsample`'s `assertAlmostEqual(post_subsample_cov, MAX_COV,
delta=MAX_COV * 0.1)` stays meaningful because the fake writes a real
gzip that `total_bases()` then re-reads through the fake stats function.

### 2.4 New cases

| # | Scenario | Purpose |
|---|---|---|
| 5 | `--nominal-size 0` | The false arm of `bases / n if n else 0.0` — the **only** genuinely uncovered branch in C2 |
| 6 | Exactly `max_cov` (51 000 reads → 300×) | Pins `>` vs `>=` in the subsample gate |
| 7 | Exactly `min_cov` (5 100 reads → 30×) | Pins `<` vs `<=` in the low-coverage gate |
| 8 | Stats table with reordered + extra trailing columns | Pins `header.index('sum_len')`. Newer seqkit adds `Q1`/`Q2`; the comment in `total_bases()` invites a hard-coded column 5 |
| 9 | Assert the recorded `bash -c` string shape | `seqtk sample -s 42`, the reads path, the fraction, `\| gzip > <out>`. **This is the CLI contract mocking otherwise erases — assert it rather than lose it** |
| 10 | `side_effect=CalledProcessError` | Pins the module docstring's "non-zero exits are reserved for unexpected tool failures", which Nextflow's `errorStrategy 'ignore'` depends on |

## 3. Changes to `scripts/tests/test_parse_samplesheet.py`

Same conversion, **no mocking needed** — C1 is pure stdlib and its tests
are not failing, merely unmeasured.

Keep the file's existing **pytest style** (plain functions and classes,
`tmp_path`, bare `assert`). Do not convert it to `unittest`; the two
styles coexist in the suite already.

- Replace the `run()` helper with an importlib load +
  `with patch('sys.argv', argv): parse_samplesheet.main()`.
- `parse_samplesheet.py` fails via `_fail()` → `sys.exit(1)`, so the 11
  failure tests become `with pytest.raises(SystemExit) as exc:` and
  assert `exc.value.code != 0`.
- stderr moves from `result.stderr` to pytest's `capsys` fixture. All
  existing stderr-content assertions map across one-for-one.

New cases, for C1's uncovered branches:

| # | Scenario | Branch |
|---|---|---|
| 1 | Sheet with no header row at all | `if not reader.fieldnames` |
| 2 | `reads` cell empty | `if not reads_raw` |
| 3 | Stray pipe in `reads` (`a.fastq.gz\|\|b.fastq.gz`) | empty entry after `split('\|')` |
| 4 | Valid ISO `sample_receipt_date` | The success arm of `date.fromisoformat` — currently only the `ValueError` arm is exercised |

## 4. `# pragma: no cover` on the `__main__` guards

Any importlib-loaded module leaves `if __name__ == '__main__':`
permanently partial;
[`bin/bin_target.py`](../../bin/bin_target.py) and
[`bin/plastid_canonicalise.py`](../../bin/plastid_canonicalise.py) already
carry this artifact today. Add `# pragma: no cover` to that guard in all
four `bin/*.py` so the report reads as a true 100 % rather than "100 %
except a line no test can reach".

One line per file, no behaviour change. Do this **last**, so the
preceding work is measured honestly first.

## 5. Reconcile the spec — three documents disagree

Not optional, and not deferrable: this task's correctness depends on
which of them is authoritative, and the inconsistency will mislead
someone else. Fix it here regardless of implementation detail.

| Source | Where do unit tests live? | Where do they run? |
|---|---|---|
| [spec §5a](../../spec/05-test-data.md#5a-tests) | `scripts/tests/`, "imported directly (no Nextflow harness)" | unstated |
| [spec §5b item 4](../../spec/05-test-data.md#5b-ci) | `scripts/tests/` | "host Python (not per-image containers)" |
| [spec §2.2](../../spec/02-stages.md#22-custom-logic-components) | `tests/unit/` | "the component's own container … matching the runtime environment exactly" |
| [`scripts/pytest.sh`](../../scripts/pytest.sh) + [`scripts/tests/README.md`](../../scripts/tests/README.md) (as built) | `scripts/tests/` | one shared `neoformit/daff-wf5-scripts:test` image |

Note that §5a's "imported directly" already describes the change this
task makes — the subprocess-invoking modules were the deviation.

Reconcile to **what is built**: `scripts/tests/`, one shared test image.
Correct §5b item 4 and §2.2 accordingly. The "host Python" claim in §5b
is now false in any case —
[task 23](23_bin_target_recalibration.md) added `mappy` to
`requirements-test.txt` precisely so the shared image matches the
BIN_TARGET container's minimap2 version.

Also delete or repoint the stale `tests/unit/`
directory, which contains only a `README.md` and is what
`.github/workflows/tests.yml`'s dead "Unit tests" step looks for. Its
line 18 — *"Mock external tool calls (seqkit, seqtk, minimap2, etc.) at
the boundary; test only the Python logic"* — is the policy this task
implements and should be preserved into
[`scripts/tests/README.md`](../../scripts/tests/README.md), not lost with
the directory.

## 6. Decisions taken, recorded so they are not reopened

### 6.1 Why not install seqkit + seqtk into the test image

The obvious alternative — follow the `mappy` precedent and add the real
tools to the `TEST=1` Docker layer — is rejected:

- **The precedent does not transfer.** `mappy` is a pip wheel. `seqkit`
  is an amd64-only Go binary distributed as a GitHub release tarball
  (needs `curl` + `ca-certificates`, and a `TARGETARCH` branch the image
  does not currently need). `seqtk` has no binary release at all: source
  tarball + `make` + zlib, and the Dockerfile *purges* its toolchain
  immediately after the mappy build.
- **It creates a second hand-maintained pin.** Matching seqkit 2.13 /
  seqtk 1.5 to the out-of-band mulled production image, with nothing
  enforcing the match ([rule 10](../../CONSTITUTION.md) versioning is about
  the reference bundle; there is no equivalent guard here).
- **~60–80 MB and a curl+make layer** on an image every developer builds
  locally on first `scripts/pytest.sh` run.
- **It still leaves C2 at 0 % measured coverage** — §1.2 is untouched by
  installing the tools.

It would only be the right answer if the goal were detecting seqkit
output-format drift in unit tests. That goal belongs to the nightly
integration surface ([rule 15](../../CONSTITUTION.md)); see §8.

### 6.2 Why not PATH shims

Writing fake `seqkit`/`seqtk` executables into a temp dir and prepending
to `PATH` preserves the subprocess invocation, but: it encodes the same
assumed output format, only in bash and further from the assertions;
it is no less a mock; and it also leaves coverage at 0 %. Strictly less
value for more machinery.

## 7. Exit criteria

- `bash scripts/pytest.sh` green **across the whole suite** — no seqkit
  failures, no skips.
- All four `bin/*.py` appear in the branch-coverage report at **100 %**
  with an empty `Missing` column.
- The 4 existing C2 cases and 20 existing C1 cases (9 happy-path,
  11 failure) keep their assertions; none is weakened to accommodate the
  conversion.
- No change to any `bin/*.py` behaviour — the §4 pragma comments are the
  only edits outside `scripts/tests/` and `spec/`.
- **No Docker image rebuild required**; `scripts/requirements.txt` and
  `scripts/requirements-test.txt` unchanged.
- `flake8 bin/ scripts/` clean at 79 cols.
- `-stub-run` still green (should be untouched — confirm, don't assume).
- [spec §5b item 4](../../spec/05-test-data.md#5b-ci) and
  [spec §2.2](../../spec/02-stages.md#22-custom-logic-components) reconciled
  per §5; `tests/unit/` resolved; the boundary-mocking policy recorded in
  [`scripts/tests/README.md`](../../scripts/tests/README.md).
- The [`tasks/todo.md`](../todo.md) seqkit entry removed; the §9 entry added.

## 8. Residual risk — state it in the test module docstring

Faking the boundary stops unit tests catching three regression classes.
Record them in the file so the next reader knows what the tests do
**not** cover:

1. A **renamed** `sum_len` column in `seqkit stats -T` output (case 8
   pins column *position* independence, not the name).
2. Real-tool edge behaviour on real input — most importantly whether a
   future seqkit emits a data row for a zero-byte gzip.
3. The unquoted f-string interpolation of `args.reads` / `args.out_fastq`
   into the `bash -c` line, which breaks on any path containing a space
   or shell metacharacter. Case 9 pins the string we *build*, not that a
   shell accepts it.

Per [rule 15](../../CONSTITUTION.md) these belong to the nightly
`-profile integration` surface, which runs the real
`neoformit/daff-wf5-coverage-gate:python3.12_seqkit2.13_seqtk1.5`
container.
[`tests/integration/assertions.sh`](../../tests/integration/assertions.sh)
already asserts `sample_status.json` exists with `.status == "ok"` for
all three samples; any format break makes `coverage_gate.py` raise, the
file goes missing, and the assertion fails. Covered at the right
cadence, with no fourth test surface.

## 9. Follow-up — a gap this work exposes

**The real `seqtk sample | gzip` subsample branch is exercised by no
test surface at all.** `assertions.sh` checks only
`.status == "ok"`, which is identical for passthrough and subsample, and
no integration fixture appears to exceed its `max_cov` (300 / 500 / 300
against 17 kb / 150 kb / 400 kb nominal per
[spec §2.1.1](../../spec/02-stages.md#211-per-target-limits)). If none
does, the real subsample path has never run outside production.

File as a `tasks/todo.md` entry, **do not fix here**. Two cheap options
for whoever takes it:

- assert `coverage.json .subsampled` / `.pre_subsample_cov` in
  `assertions.sh`, so the branch actually taken is recorded and pinned;
  or
- force subsampling in the integration profile via a
  `--coverage_limits.animal_mt.max_cov` override, or a fourth fixture.

This is worth booking *because* the boundary is now explicitly faked:
the unit tests never covered the real call, and making the mock visible
should not let the gap be absorbed silently.

## 10. Not in scope

- **Wiring `scripts/tests/` into CI.** `.github/workflows/tests.yml`'s
  `python-tests` job runs only flake8; its "Unit tests" step points at
  `tests/unit/`, which holds no tests, so it always no-ops. Already
  tracked in [`tasks/todo.md`](../todo.md). A green, measured suite is the
  **precondition** for that task — land this first, then wire it.
- **Any change to `bin/coverage_gate.py` or `bin/parse_samplesheet.py`
  behaviour.** If a new test reveals a genuine defect, report it — do not
  fix it inside this task.
- **Retuning `params.coverage_limits`.** That is
  [task 25](25_coverage_gate_carryover.md)'s subject.
- **Adding a coverage config file** (`.coveragerc` / `pyproject.toml`)
  to enable subprocess tracing. That would fix §1.2 while leaving §1.1
  broken, and adds config to a repo that currently has none —
  [rule 19](../../CONSTITUTION.md).

## Outcomes

- Premise (§1.2) verified before changes: `bash scripts/pytest.sh
  scripts/tests/test_plastid_canonicalise.py` reported only
  `plastid_canonicalise.py` in the coverage table.
- `TestCoverageGate` in `test_coverage_gate.py` rewritten to the
  in-process idiom with a data-driven `_fake_run`/`_stats_table`/
  `_fake_seqtk` (§2.2); all 4 existing cases kept their bodies and
  arithmetic, `rc == 0` added; cases 5-10 added per §2.4. Case 8
  (`sum_len` column reordering) is asserted directly against
  `total_bases()` rather than through `main()` — simpler, and pins the
  same rule (`header.index('sum_len')`).
- `test_parse_samplesheet.py` converted to in-process loading; the
  `run()` helper now drives `parse_samplesheet.main()` under
  `patch.object(sys, 'argv', ...)`, catches the `SystemExit` `_fail()`
  raises, and folds `(code, capsys stderr)` back into the same
  `Result` shape the existing 20 tests already assert against, so no
  existing assertion needed to change. 4 new cases added per §3.
- `# pragma: no cover` added to the `__main__` guard in all four
  `bin/*.py` (§4), last, after the suite was green.
- `tests/unit/` deleted outright (held only a README, no tests); its
  boundary-mocking policy line folded into `scripts/tests/README.md`.
  `spec/05-test-data.md` §5b item 4 and `spec/02-stages.md` §2.2
  reconciled to describe what's built: `scripts/tests/`, one shared
  `neoformit/daff-wf5-scripts:test` image (§5).
- `tasks/todo.md`: the task-23 seqkit entry removed; a new entry added
  recording that the real `seqtk sample | gzip` subsample branch is
  exercised by no test surface (§9) — not fixed here per instruction.
- Exit criteria confirmed: `bash scripts/pytest.sh` — 129 passed, all
  four `bin/*.py` at 100% branch coverage with an empty `Missing`
  column; `flake8 bin/ scripts/` clean at 79 cols;
  `nextflow run . -profile stub -stub-run` still green.
- No deviations from the brief; no `bin/*.py` behaviour changed outside
  the §4 pragma comments.
