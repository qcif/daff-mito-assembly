# To-do items

Grouped by when an item should be picked up. Everything under
**Carry-forward into future task briefs** is a precondition discovered
while building an earlier stage — fold it into the named task's brief
when that task is drafted, then delete it from here. The remaining
sections are unscheduled backlog.

## Carry-forward into future task briefs

### Run-level provenance (task 45_run_report.md)

- (task 43a, 2026-09-02) **`pipeline_commit` is `"unknown"` in current
  real output.** `workflow.commitId` is null when the pipeline runs
  from a working directory rather than a cloned project — verified on
  a real `-profile integration` run of `INT-ANIMAL-01`. The per-sample
  provenance panel
  (`scripts/report/templates/components/provenance.html`)
  renders that honestly rather than hiding the row; task
  45_run_report.md owns run-level provenance and should either resolve
  a real commit hash at invocation time or document why `"unknown"` is
  an acceptable steady state for non-cloned runs.

### `ORGANELLE_MAP`

- Its rendering pass should walk both `path1` and `path2` when
  `plastid_isoforms/` is present, per
  [spec §3.6 step 5](../spec/03-organelles.md#36-plastid-quadripartite-canonicalisation).

### Per-sample report — per-position coverage depth track

- (task 43b, 2026-09-04) spec §6a.2's Assembly tab asks for a
  per-contig coverage plot; no stage in this pipeline computes
  per-position depth, so task 43b_report_stage_tabs.md renders
  **per-contig mean coverage** (Flye's `assembly_info.txt`, one figure
  per contig) as a bar chart instead — that is fully served by data in
  hand. A real depth track needs a new alignment stage (e.g. a
  `samtools depth` pass over the recruited reads against the selected
  target contig) — a pipeline change, not a reporting one. Scope as its
  own task if a real per-position track is wanted.

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

- (task 25, 2026-08-03; superseded 2026-09-02 by task 35) **`INT-PLANT-01-mt`
  needs replacing with a deeper-sequenced plant sample.** Task 35's two-floor
  gate means this fixture no longer soft-fails: at **28.48×** mitochondrial
  depth (after task 28's panel change and RECRUIT MAPQ repair) it clears the
  10× hard floor and assembles as `low_coverage`, sitting between the two
  floors rather than below both. The recruited pool is still 70.5 % plastid
  and only ~726 kb of it genuinely mitochondrial, so this remains correct
  behaviour, not a fixture bug — a sample that only clears the *hard* floor
  is exactly what `low_coverage` exists to describe, and it is worth *far*
  less as a validation fixture than one that clears the *warn* floor: it can
  exercise the `plant_mt` assembly arm but not validate its output quality.

  Consequence: `INT-PLANT-01-mt` now reaches METAFLYE, BANDAGE_NG, BIN_TARGET
  and beyond for the first time (task 35 §8) — the `plant_mt` branch of C3,
  including the sibling-organelle discrimination
  [task 23](completed/23_bin_target_recalibration.md) built and the
  `emit: all` multi-contig path, is exercised by real data again, not only
  unit tests. But it is still left off `tests/integration/assertions.sh`'s
  `ASSEMBLING_SAMPLES`, so none of the downstream biology assertions
  (assembly bounds, binning bounds, barcode/annotation checks) run against
  it yet — wiring those up, replacing this fixture with one that clears the
  *warn* floor, or both, is **task 36**, not resolved here. [Task
  26](26_binning_marker_genes.md)'s §2 go/no-go gate also still lacks the
  real-data evidence it was scoped to consume until task 36 lands; re-read
  that task's premise before starting it.

  A replacement that clears 30× *mitochondrial* depth after the sibling
  split would let the fixture validate rather than merely exercise the
  `plant_mt` arm. `tests/integration/expected/plant_mt/{assembly,bin}_bounds.json`
  are retained but unasserted, and `tests/integration/assertions.sh`
  restores the `plant_mt` blocks by adding the sample back to
  `ASSEMBLING_SAMPLES`. See
  [spec §5](../spec/05-test-data.md) for the fixture-staging procedure.

### Barcode extraction (`EXTRACT_BARCODES` / C5)

- (task 30, 2026-08-05) **Only 1/6 `animal_mt` barcode loci clear the
  60% identity floor on `INT-ANIMAL-01` — this is a reference-panel
  problem, and the fix is a panel change, not a threshold nudge.**
  Diagnosed in [task 30 §10.1](completed/30_unified_locus_pass.md):
  the assembly is 98.281% nucleotide-identical to *Acyrthosiphon
  pisum* over 12,277 bp (i.e. near-conspecific quality — better input
  data cannot help), while protein identity to the panel is 42–65%
  because the `animal_mt` panel's 10 reps/gene span all of Metazoa
  (penguin, flatfish, flatworm, crab, beetles, mosquito, earwig, moth,
  leafhopper, psyllid) and contain **no aphid**. Query coverage stays
  at **85–99.5%** across all six loci — the genes are complete and
  intact; only the nearest reference is distant.

  Three separable pieces of work, in rough priority order:

  1. **Reconsider what the identity floor is for.** It currently
     conflates "is this barcode real and intact" (wanted) with "does
     our bundle happen to cover this taxon" (an artifact). For
     biosecurity the under-referenced taxon is the case that matters
     most and this floor preferentially rejects it. The ORF/
     internal-stop check already tests intactness independently.
     Consider gating on **query coverage** instead of, or alongside,
     identity — it cleanly separated real-but-distant (0.85–0.99) from
     everything else on both fixtures. `--min-identity` is currently
     hard-coded to 60 in `modules/local/extract_barcodes.nf`; it should
     become a `nextflow.config` param (per-target, likely per-gene)
     before any sweep.
  2. **Taxonomically stratify the `animal_mt` panel.** Metazoa needs
     far more than the 10 reps/gene that suffice for angiosperms —
     either many more reps, or deliberate order/family-level coverage
     of taxa likely to be submitted (insects especially), rather than
     task 29's diversity-maximising one-per-genus sampling. Note this
     interacts with task 29's §3 go/no-go caveat about representative
     *count* being confounded with panel *breadth*.
  3. **Only then sweep the threshold**, against the deeper benchmark
     set (see Benchmark data above) — tuning against this single
     fixture would be circular (spec §5.1).

  Not a code defect; `tests/integration/assertions.sh` already treats
  per-locus misses as `WARN` and hard-fails only on zero recovery.

  **(2026-08-15) Re-measure before acting on any of the above.** The
  42–65 % protein identities were measured while `MINIPROT_CDS` was
  aligning under the standard genetic code instead of the target's —
  see task 38_miniprot_genetic_code.md. On a client `animal_mt`
  intercept, correcting the table raised per-gene identity by roughly
  6–7 points across the panel and removed spurious `Frameshift` flags
  entirely. Panel breadth is still a genuine and probably dominant
  factor, but the figures above overstate it by an unknown margin and
  item 1's premise (identity conflating "real and intact" with "covered
  by our bundle") should be re-argued on corrected numbers. Query
  coverage, the separator item 1 proposes gating on, is unaffected by
  the bug.

  **(2026-09-01) Re-measured on `INT-ANIMAL-01` after task 38 landed:
  2/6 now clear the 60% floor, up from 1/6.** Per-locus protein identity
  is now 48.4–72.8% (COX1 0.728, CYTB 0.615, COX2 0.605, COX3 0.531, ND1
  0.525, ATP6 0.484) — up from the pre-fix 42–65%, and COX2/CYTB now
  clear the floor where nothing but one locus did before. `ATP6`'s best
  panel hit still carries a `Frameshift=3` flag even under the now-correct
  table 5 — but the `StopCodon=3` flag that co-occurred with it under the
  wrong table 1 is gone (checked directly: table 1 gave 62 alignments /
  123 `StopCodon` / 125 `Frameshift`; table 5 gives 74 alignments / 26
  `StopCodon` / 148 `Frameshift`). Read that as: task 38 fixed the
  false-stop artifact this task exists to fix, but the residual
  `Frameshift` on divergent low-identity hits like `ATP6` (best hit only
  48.4% identity — no aphid in the panel) is very plausibly a genuine
  indel between distant homologs, not a translation-table bug — i.e. it's
  item 2 below (panel breadth), not something task 38 owns. `COX1` at
  72.8% identity still fails on `internal_stop_codon` despite clearing
  the floor — worth a follow-up look at whether that's a genuine
  premature stop or a fragment-boundary artifact of a partial alignment,
  independent of the identity-floor discussion above.
- (task 30, 2026-08-05) **`codon_blocks()` treats intron ops (`N`/`U`/`V`)
  the same as frameshift ops (`F`/`G`): dropped from translation, not
  spliced into a proper multi-exon reading frame.** This is correct
  for every locus in the current barcode panel (verified
  intron-free for the tested fixtures) but is a simplification for
  genuinely intron-containing genes — `plant_mt`'s `nad1` is the known
  case (trans-spliced, multiple exons). `INT-PLANT-01-mt` never
  reaches `EXTRACT_BARCODES` (soft-fails the coverage gate — see the
  Integration fixtures item above), so this has never been exercised
  against a real intron. Revisit if/when a `plant_mt` fixture that
  reaches this stage lands.

### Annotation (`ANNOTATE` / C8)

- (task 39, 2026-09-01) **miniprot's `animal_mt` CDS calls over-extend
  across frameshifts, and that over-extension is what creates false
  positional collisions with neighbouring genes.** On `INT-ANIMAL-01`
  miniprot's winning `ATP6` spans 1200–1847 while MITOS2 calls `atp6`
  at 1197–1615 — ~230 bp longer — and every `ATP6` hit in `cds.gff`
  carries `Frameshift=3`. The over-extended tail is what `atp8_1`
  collided with (64% of its length), and only the §2.4 fraction
  threshold plus guarding-before-clustering saved `ATP8`. Two
  consequences worth following up: (a) a heavily-frameshifted
  alignment is a lower-confidence feature and should score as such in
  task 40_annotation_confidence_scores.md; (b) if the over-extension
  is systematic it may warrant trimming miniprot spans at frameshift
  boundaries, which would make the rescue guard's job easier and the
  gene map more accurate. Not acted on in task 39 — the rescue works
  as-is, and trimming CDS spans would touch the barcode-coherence
  invariant (task 30) that stage 12 output is used verbatim.

  **(2026-09-01) (a) addressed by task 40.** A heavily-frameshifted
  miniprot span now does read as lower-confidence: it translates into
  a garbled protein past the frameshift point, which shows up as
  reduced `pident`/`qcovhsp` against the panel rather than a silently
  "corrected" call. (b) — trimming the over-extended span itself — is
  untouched and still open.

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
- (task 39, 2026-09-01) **Relative links inside completed task files
  break on move.** Task briefs are written with `](../CONSTITUTION.md)`
  / `](../spec/...)`, correct from `tasks/` but not from
  `tasks/completed/`. Seven files are affected (21, 29, 30, 31, 38, 39,
  41 — 100+ links). The spec already avoids the mirror-image problem by
  never linking *to* task files. Either move to root-relative links, add
  the extra `../` as part of the move step, or add a CI link-checker —
  whichever is picked, do it in one pass rather than per-task, since a
  half-converted set is worse than a consistently-wrong one.
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
