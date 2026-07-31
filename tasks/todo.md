# To-do items

- Refactor integration tests and assertions.sh into an nf-test suite
- At the end, do we extract all barcodes or just one? Probably better to
  extract all.
- Check Flye polishing requirement - maybe not necessary any more
- Think about iterative recruitment for tricky/cryptic samples
- Switch flye to --nano-corr if reads are confirmed error-corrected
- If we pass a genome_size hint, that could possibly be a samplesheet column in
  future
- When the `ORGANELLE_MAP` task is drafted, its rendering pass should walk
  both `path1` and `path2` when `plastid_isoforms/` is present, per
  [spec §3.6 step 5](../spec/03-organelles.md#36-plastid-quadripartite-canonicalisation).
- Examine script config (e.g. hard-coded config in Python scripts like
  bin_target.py) that should perhaps be set as NF params.
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
- (task 20 integration run, 2026-07-31) `INT-ANIMAL-01`'s real
  `-profile integration` run has `bin_metadata.json.circular == false`
  — the task 18 end-overlap circularity check (`check_circularity` in
  `bin/bin_target.py`) isn't detecting circularity on the real
  METAFLYE assembly for this fixture. Pre-existing task 18 behaviour,
  unrelated to task 20 (C4); the `assertions.sh` circularity check
  (task 18 block) currently fails on this. Needs investigation —
  possibly `END_OVERLAP_WINDOW`/`END_OVERLAP_MIN_IDENTITY` tuning, or
  Flye not producing a genuine end-overlap on this fixture.
- (task 20 integration run, 2026-07-31) On the real `-profile
  integration` run, **both** `INT-PLANT-01-pt` and `INT-PLANT-01-mt`
  have `n_target_selected == 0` — `select_primary`/`classify_contigs`
  (task 18, `bin/bin_target.py`) selects no `target_candidate` contigs
  on either fixture (`ref_identity_pct`/`aligned_frac`/ORF thresholds
  not cleared). Pre-existing task 18 behaviour, unrelated to task 20.
  For `plant_pt` this is currently masked: C4's plastid
  canonicalisation runs unconditionally on the raw
  `assembly_graph.gfa` regardless of C3's primaries selection, and on
  the `canonical` branch overwrites `target.fasta` with `path1.fasta`
  — so `target.fasta` ends up populated and correct-length even
  though C3 itself picked nothing. `plant_mt` has no such rescue path,
  so its `target.fasta` stays empty. Needs investigation into task
  18's classification thresholds (`MIN_REF_IDENTITY`,
  `COVERAGE_SPIKE_FACTOR`, `MIN_ORF_AA`) against these real fixtures.
- (task 20, 2026-07-31) `nextflow run . -profile integration -resume`
  did **not** re-execute `BIN_TARGET` after editing
  `bin/plastid_canonicalise.py`'s depth-tag regex (an in-place edit to
  an already-existing, already-run `bin/` script) — it reused the
  stale cached `INT-PLANT-01-pt` result from before the fix. Had to
  `rm -rf work/ .nextflow*` and rerun from scratch to get a trustworthy
  result. Investigate whether Nextflow 25.10.2's default cache mode
  actually hashes `bin/` script contents into the task hash on this
  pipeline (per `main.nf`/`nextflow.config` — no explicit
  `cache`/`resume` settings found); if it doesn't by default here,
  either enable a stricter cache mode or document that `-resume` isn't
  safe after editing `bin/*.py` files.
