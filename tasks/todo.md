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

### Integration fixtures

- (task 25, 2026-08-03) **`INT-PLANT-01-mt` needs replacing with a
  deeper-sequenced plant sample.** Once C2 gates on target-assigned
  bases, this fixture soft-fails at 27.29× against the 30× MIN — its
  recruited pool is 70.5 % plastid and only ~726 kb of it is genuinely
  mitochondrial. The soft-fail is the *correct* verdict, so the fixture
  is not tunable: it cannot support any assertion downstream of the
  coverage gate.

  Consequence: `INT-PLANT-01-mt` no longer exercises METAFLYE,
  BANDAGE_NG, BIN_TARGET or anything after them. The `plant_mt` branch
  of C3 — including the sibling-organelle discrimination
  [task 23](completed/23_bin_target_recalibration.md) built and the
  `emit: all` multi-contig path — is now covered only by unit tests.
  [Task 26](26_binning_marker_genes.md)'s §2 go/no-go gate also loses
  the real-data evidence it was scoped to consume; re-read that task's
  premise before starting it.

  A replacement must clear 30× *mitochondrial* depth after the sibling
  split — note that [task 28](28_plastid_masked_mt_panel.md) does **not**
  rescue this fixture (a masked broad panel puts it at 28.37×, still
  under MIN), so the two are independent. `tests/integration/expected/plant_mt/{assembly,bin}_bounds.json`
  are retained but unasserted, and `tests/integration/assertions.sh`
  restores the `plant_mt` blocks by adding the sample back to
  `ASSEMBLING_SAMPLES`. See
  [spec §5](../spec/05-test-data.md) for the fixture-staging procedure.

## Open design questions

- At the end, do we extract all barcodes or just one? Probably better to
  extract all.
- Consider iterative extension recruitment for tricky/cryptic samples (risks
  walking into nuclears reads via NUMTs - hard to solve)
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
  step looks for `tests/unit/test_*.py` (that directory was removed —
  task 27 — since it held only a stale README, no tests) —
  `scripts/tests/` (the real bin/*.py test suite, run via
  `scripts/pytest.sh` inside the `neoformit/daff-wf5-scripts:test`
  container) is never actually invoked by CI. Wire that up.
- `nextflow -resume` didn't re-run `BIN_TARGET` after an in-place edit
  to `bin/plastid_canonicalise.py` — reused a stale cached result.
  Confirm whether this pipeline's cache mode actually hashes `bin/`
  script contents into the task hash; if not, document that `-resume`
  isn't safe after editing `bin/*.py`.
- (task 27, 2026-08-03) **The real `seqtk sample | gzip` subsample
  branch in `bin/coverage_gate.py` (C2) is exercised by no test surface
  at all.** Unit tests fake the tool boundary (necessarily — see
  [task 27](completed/27_unit_test_boundary_mocking.md) §6);
  `tests/integration/assertions.sh` checks only
  `sample_status.json .status == "ok"`, which is identical for
  passthrough and subsample, and no integration fixture appears to
  exceed its `max_cov`. Two cheap options: assert
  `coverage.json .subsampled` / `.pre_subsample_cov` in
  `assertions.sh` so the branch actually taken is pinned, or force
  subsampling in the integration profile via a
  `--coverage_limits.animal_mt.max_cov` override or a fourth fixture.
