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
- `nextflow -resume` didn't re-run `BIN_TARGET` after an in-place edit
  to `bin/plastid_canonicalise.py` — reused a stale cached result.
  Confirm whether this pipeline's cache mode actually hashes `bin/`
  script contents into the task hash; if not, document that `-resume`
  isn't safe after editing `bin/*.py`.
