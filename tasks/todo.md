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

### Benchmark data (distinct from the fixtures below)

- (2026-08-04) **Source a proper benchmark set.** The Tier 2 fixtures
  are correctness tests and cannot measure workflow performance
  ([spec §5.1](../spec/05-test-data.md#51-the-integration-fixtures-are-correctness-tests-not-benchmarks)):
  subsampled to a 60-min CI budget, old public SRA predating the
  R10.4.1 SUP chemistry the pipeline assumes, one sample per target,
  and `plant_pt` / `plant_mt` are the same accession. Every threshold
  in `nextflow.config` is provisional until swept against real data,
  and performance should be materially better than any fixture figure.

  Needs: modern ONT (R10.4.1, SUP-basecalled) at realistic submission
  depth, more than one taxon per target, and `plant_pt` / `plant_mt`
  from *different* samples. This is what
  [spec §9](../spec/07-open-questions.md#9-fine-tuning-post-prototype-benchmarking)
  tunes against — several of its items are currently unanswerable
  without it.

- (2026-08-04) **Record fixture provenance.** Each
  `tests/integration/expected/*/` should carry instrument, flowcell
  chemistry, basecaller + model and run date alongside the organism and
  SRA accession it already records. Without it we cannot judge how far
  fixture behaviour generalises, or notice when it stops being
  representative.

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
  split — note that [task 28](completed/28_plastid_masked_mt_panel.md) does **not**
  rescue this fixture. After task 28's panel change *and* its §10 RECRUIT
  MAPQ repair the sample sits at **28.48×** against the 30× MIN, with
  carry-over down to 0.6887. Better than the 27.29× above, still a
  soft-fail. The two items stay independent. `tests/integration/expected/plant_mt/{assembly,bin}_bounds.json`
  are retained but unasserted, and `tests/integration/assertions.sh`
  restores the `plant_mt` blocks by adding the sample back to
  `ASSEMBLING_SAMPLES`. See
  [spec §5](../spec/05-test-data.md) for the fixture-staging procedure.

### Reference data

- (task 28, 2026-08-03) **Mask the `plant_mt` panel against the full
  RefSeq plastid set, not the 101-genome `plant_pt` panel.** Task 28
  masks NUPTs by aligning RefSeq Viridiplantae mitogenomes against
  `recruit/plant_pt.fa` — 101 genomes, 15 Mb. NUPTs derive from the
  *donor lineage's own* plastome, so a NUPT from a lineage absent from
  those 101 genomes is not found and stays unmasked. This is the most
  likely explanation for the residual signal task 28 measured: under
  C3's own `map-ont` metric, plastid contigs still score 0.31–0.35
  against the masked panel (up from 0.01–0.06 against the Vigna seed).
  The sibling margins stay wide (−0.65, −0.69) so nothing is broken,
  but a fuller mask would push them back toward the −0.9 range.

  `validate/refseq_pt.fa` (2.36 GB, 15 233 genomes) is the obvious
  reference. **It is not free:** an `asm20 -I 8G` alignment of the
  259 MB mitogenome set against it was OOM-killed on a 128 GB
  workstation. Needs a chunked or lower-`-I` strategy, and a re-run of
  task 28 §5.1's measurements to show the margin actually improves
  before the extra build cost is worth taking.

- (task 28, 2026-08-03) **`plant_mt.mmi` is now 536 MB** (was 1.5 MB),
  and the bundle ~1.9 GB (was 865 MB). Every `plant_mt` sample's
  RECRUIT, C2 coverage split and C3 binning loads that index. Confirm
  the peak RSS fits the CI runner budget in
  `conf/integration.config`'s `resourceLimits` and the production
  Azure profile, and consider whether `minimap2 -I` needs pinning so
  index construction and consumption agree.

### Recruitment

- (task 28 §10, 2026-08-04) **`params.recruit_thresholds` ship at 0/0 —
  a deliberate no-op.** The merged-aligned-extent floors that replaced
  RECRUIT's MAPQ filter cannot currently be raised: every non-zero value
  moves an integration fixture outside its expected coverage bounds
  (`animal_mt` and `plant_pt` leave their bands at `min_aligned_frac`
  0.1–0.2; `plant_mt` collapses 28.48× → 7.34×). Revisit once the
  deeper `plant_mt` fixture above lands, then sweep per
  [spec §9 item 1](../spec/07-open-questions.md#9-fine-tuning-post-prototype-benchmarking).
  Purity data for the sweep is in task 28 §10.2.

- (task 28 §10, 2026-08-04) `INT-ANIMAL-01`'s assembly length is not
  stable across recruitment changes — it exceeds `max_cov`, so a
  different recruited pool gives a different subsample and a different
  assembly (348,580 → 215,078 bp on this change, both inside bounds).
  If a tighter `assembly_bounds.json` is ever wanted for `animal_mt`,
  the subsample seed needs pinning first.

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
