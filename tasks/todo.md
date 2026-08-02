# To-do items

Grouped by when an item should be picked up. Everything under
**Carry-forward into future task briefs** is a precondition discovered
while building an earlier stage — fold it into the named task's brief
when that task is drafted, then delete it from here. The remaining
sections are unscheduled backlog.

## Carry-forward into future task briefs

### `COLLATE` (P4)

- (task 24, 2026-08-03) `COLLATE`'s inputs must handle the
  withheld-substitution state: a `plant_pt` sample with
  `plastid_canonicalisation.substitution_applied: false` and
  `substitution_withheld_reason: "no_c3_selection"` emits an **empty**
  `target.fasta` alongside a populated `plastid_isoforms/`. That must
  route to the `no_assembly` bundle path, not to the full bundle with
  empty downstream outputs
  ([task 24 §3.2](completed/24_plastid_substitution_guard.md)).

### Per-sample report (P4)

- (task 24, 2026-08-03) The "Assembly quality assessment" section must
  render the withheld-substitution warning and the `target_source`
  provenance line, per
  [spec §3.6](../spec/03-organelles.md#36-plastid-quadripartite-canonicalisation)
  and [spec §6a.2](../spec/06a-reports.md).

### `ORGANELLE_MAP`

- Its rendering pass should walk both `path1` and `path2` when
  `plastid_isoforms/` is present, per
  [spec §3.6 step 5](../spec/03-organelles.md#36-plastid-quadripartite-canonicalisation).

## Open design questions

- At the end, do we extract all barcodes or just one? Probably better to
  extract all.
- Consider iterative recruitment for tricky/cryptic samples
- Switch flye to --nano-corr if reads are confirmed error-corrected (maybe just
  make this a WF param)
- If we pass a genome_size hint, that could perhaps be a samplesheet column in
  future
- Examine script config (e.g. hard-coded config in Python scripts like
  bin_target.py) that should perhaps be set as NF params.
- (task 23, 2026-08-02) `params.bin_target_thresholds` values are
  prototype defaults calibrated on three fixtures — see
  [spec §9 item 10](../spec/07-open-questions.md#9-fine-tuning-post-prototype-benchmarking).
  The `low_coverage_fraction` of 0.05 has never fired on real data, so
  it is untested against a genuine NUMT/NUPT.

## Test + CI backlog

- Refactor integration tests and assertions.sh into an nf-test suite
- Flake8 all Python scripts with new 79 line length
- Extract the container-coverage and bare-params heredoc checks in
  `.github/workflows/tests.yml` (`lint` job) into standalone scripts
  (e.g. `scripts/check_containers.sh`, `scripts/check_bare_params.sh`)
  so they can be run locally as well as in CI.
- `tests.yml`'s `python-tests` job only runs flake8; its "Unit tests"
  step looks for `tests/unit/test_*.py`, which doesn't exist —
  `scripts/tests/` (the real bin/*.py test suite, run via
  `scripts/pytest.sh` inside the `neoformit/daff-wf5-scripts:test`
  container) is never actually invoked by CI. Wire that up.
- (task 23, 2026-08-02) `scripts/tests/test_coverage_gate.py` fails
  (4 tests) inside the `neoformit/daff-wf5-scripts:test` image because
  `seqkit` isn't installed there — `bin/coverage_gate.py` shells out to
  it. Pre-existing, but it means `scripts/pytest.sh` is not green across
  the whole suite. **Now owned by
  [task 27](27_unit_test_boundary_mocking.md)** — delete this entry when
  that task lands.
- `nextflow -resume` didn't re-run `BIN_TARGET` after an in-place edit
  to `bin/plastid_canonicalise.py` — reused a stale cached result.
  Confirm whether this pipeline's cache mode actually hashes `bin/`
  script contents into the task hash; if not, document that `-resume`
  isn't safe after editing `bin/*.py`.
