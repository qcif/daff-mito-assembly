# Task 25 — `COVERAGE_GATE` over-estimates coverage when a sibling organelle is recruited

**Phase:** P3 (from [spec §6](../../spec/06-phases.md)).
**Goal:** Establish whether
[`COVERAGE_GATE`](../../spec/02-stages.md#21-coverage-gate-2-stage-6) (C2)
is passing `plant_mt` samples whose true mitochondrial coverage is well
below the 30× MIN, because the recruited read pool it measures is
dominated by plastid. If confirmed, decide and implement the correction.

**Status:** investigation-first. The defect is inferred from strong
circumstantial evidence, not yet directly measured — §2 is the
measurement that must precede any change.

**Related:** found while investigating
[task 23](23_bin_target_recalibration.md). Independent of it — different
stage, different failure mode — but the two share a root cause
(RECRUIT's positive selection does not separate plant organelles).

## 1. The inference

C2 computes `estimated_cov = total_recruited_bases / nominal_organelle_size`,
with the nominal size fixed per `assembly_target`
([§2.1.1](../../spec/02-stages.md#211-per-target-limits)). It has no way to
know which organelle the recruited bases came from.

On `INT-PLANT-01-mt` the gate reported `status: "ok"`. But
[spec §3.7](../../spec/03-organelles.md#37-target-binning-criteria-per-assembly-target)
established that 155 291 bp of the 270 138 bp assembly — **57 %**, and
the two highest-coverage contigs at 137× and 168× — is chloroplast.
The genuine mitochondrial contigs assembled at **5–7×**, far under the
30× MIN the gate is supposed to enforce.

If the assembled coverage reflects the recruited pool, the gate's
estimate was inflated by plastid reads and the sample should have
soft-failed as `low_coverage`. Instead it proceeded to assembly and
produced an under-assembled ~115 kb mitogenome presented as a normal
result.

This is [principle 7](../../CONSTITUTION.md) in reverse: a degraded sample
produced output indistinguishable from a confident one. The whole point
of the coverage gate is to catch this before assembly, and here the
gate's own input is what hid it.

**Why `plant_mt` specifically:** the abundance gradient runs
plastid → mitochondrion in plant tissue, and the two share enough
homology (plus genuine NUPT insertions in the mitogenome) for plastid
reads to be recruited against `plant_mt.mmi`. A `plant_pt` run has the
gradient in its favour and is unlikely to be affected; `animal_mt`
has no sibling organelle. So this is a one-target problem, which
constrains how much machinery is justified in fixing it.

## 2. Measurement (do this first)

Before changing anything, confirm the inference on the existing fixture.
The `INT-PLANT-01-mt` recruited FASTQ and `coverage.json` are already
produced by `-profile integration`:

1. Map the recruited reads (post-gate `gated.fastq.gz`) against both
   `plant_mt.mmi` and `plant_pt.mmi`; classify each read by which panel
   it aligns to better — the same merged-interval comparison
   [task 23 §2.1–2.2](23_bin_target_recalibration.md) introduces for
   contigs.
2. Report the plastid/mito split of recruited bases, and recompute
   `estimated_cov` from the mito-assigned bases alone.
3. Compare against the assembled per-contig coverage (5–7× mt,
   137–168× pt) as a cross-check.

**Decision gate:** if mito-only estimated coverage is above MIN, the
gate's verdict was right for the wrong reason and this task reduces to
a reporting change (§3.1). If it is below MIN, the gate has a real false
`ok` and §3.2 or §3.3 is required.

Keep this as a throwaway diagnostic — do not add a script to `bin/`
before the decision gate resolves.

## 3. Options if confirmed

Presented in ascending cost. **Do not pick before §2 completes.**

### 3.1 Report-only (minimum)

Leave the gate as-is; have C3 report, post-assembly, the fraction of
assembled bases assigned to a sibling organelle, and surface a report
warning where it exceeds a threshold: *"57 % of assembled sequence is
chloroplast — mitochondrial coverage may be substantially lower than the
pre-assembly estimate."*

Task 23 already computes exactly this data for binning, so the marginal
cost is a metadata field plus a report block. Honest and cheap. Does not
prevent the wasted assembly, and the sample still reports `ok`.

### 3.2 Sibling-aware coverage estimate in C2

Have C2 split recruited reads by panel before estimating, using the same
sibling map as C3, and gate on target-assigned bases only.

Correct at the right stage and preserves the gate's purpose. Costs a
second minimap2 pass over the recruited pool inside C2 (moderate — the
pool is already small post-recruitment) and gives C2 a reference-panel
dependency it does not currently have, which is a real change to a
deliberately thin component
([§2.2 C2](../../spec/02-stages.md#22-custom-logic-components): "Python
stdlib only + `seqkit`/`seqtk` in `PATH`"). Would need a container
change.

### 3.3 Sibling-aware recruitment in RECRUIT

Filter at recruitment: drop reads aligning better to a sibling panel
than to the declared one.

Fixes the root cause and improves assembly quality, not just the gate —
METAFLYE would no longer spend its depth budget on chloroplast. But it
sits closest to [principle 5](../../CONSTITUTION.md), which is deliberate
about recruitment being *coarse* positive enrichment with precise
separation deferred to binning. Competitive filtering at RECRUIT is
still positive selection (no host reference involved), so it is not a
principle violation — but it is a design-principle boundary and needs
explicit sign-off per the amendment procedure. It also risks discarding
genuine mt reads spanning NUPT insertions, which are real mitochondrial
sequence that legitimately aligns better to the plastid panel.

**Preliminary lean:** 3.1 unconditionally (it is nearly free and
strictly increases auditability), plus 3.2 if §2 confirms a false `ok`.
Hold 3.3 for [§9 item 1](../../spec/07-open-questions.md#9-fine-tuning-post-prototype-benchmarking)
(RECRUIT filter strictness), which is already scheduled to revisit
recruitment filtering with real data — this is evidence for that item,
not a reason to pre-empt it.

## 4. Knock-on: the `plant_mt` integration fixture

[Task 23 §5.3](23_bin_target_recalibration.md) widens
`expected/plant_mt/bin_bounds.json` rather than pinning it, because
until this task resolves we cannot distinguish "C3 missed contigs" from
"the assembly genuinely lacks them at 5–7× coverage". Once resolved,
tighten those bounds to the evidence.

If §2 shows the mitogenome is simply under-recruited, the honest
conclusion may be that **`INT-PLANT-01-mt` is a poor fixture for
asserting mitogenome completeness** and should either be replaced with a
deeper-sequenced plant sample or have its assertions scoped to what it
can actually demonstrate. Per [rule 19](../../CONSTITUTION.md), a fixture
that cannot support its assertions is a long-term maintenance liability.
Flag this outcome explicitly rather than tuning bounds until they pass.

## 5. Exit criteria

- [x] §2 measurement complete, with the plastid/mito split of recruited
  bases and a recomputed mito-only `estimated_cov` recorded in this
  task's outcomes section.
- [x] A decision recorded against §3, with rationale.
- [x] If 3.1: sibling-organelle fraction in `bin_metadata.json` and a report
  warning above threshold; assertion in `assertions.sh`. *(The warning is
  specified for the report in spec §6a.2 and emitted as a metadata flag;
  the renderer itself is a P4 stub.)*
- [x] If 3.2: C2 gates on target-assigned bases; unit tests to
  [rule 14](../../CONSTITUTION.md) 100 % branch coverage including the
  "sibling panel unavailable" fallback; container change documented;
  `INT-PLANT-01-mt`'s expected gate status updated to whatever the
  corrected estimate warrants — **including `low_coverage` if that is
  the honest answer**, with `assertions.sh` updated to expect it.
  *(`low_coverage` was the honest answer and is what it now expects.)*
- [x] [spec §2.1](../../spec/02-stages.md#21-coverage-gate-2-stage-6) updated
  with the limitation and the resolution, whichever option is taken.
  *(New §2.1.5.)*

## Outcomes

Completed 2026-08-03. §2 confirmed a real false `ok`; §3.1 and §3.2 both
implemented. §3.3 held for
[spec §9 item 1](../../spec/07-open-questions.md#9-fine-tuning-post-prototype-benchmarking)
as the brief proposed.

### §2 Measurement

Throwaway diagnostic (`minimap2 -x map-ont` against `plant_mt.mmi` and
`plant_pt.mmi`, merged query intervals per read, read assigned to the
panel with the higher merged aligned fraction — the §3.7.2 metric C3
applies to contigs). `INT-PLANT-01-mt`'s gate was a passthrough, so the
recruited FASTQ and `gated.fastq.gz` are byte-identical; the recruited
file was used. Nothing was added to `bin/` before the decision gate
resolved.

**Recruited pool, 1366 reads / 37 038 413 bp** (matching the gate's own
`total_recruited_bases` exactly):

| Class | Reads | Bases | Share |
|---|---|---|---|
| plastid-assigned | 980 | 26 119 082 | **70.5 %** |
| mito-assigned | 385 | 10 916 501 | 29.5 % |
| tie | 1 | 2 830 | 0.0 % |
| unaligned | 0 | 0 | 0.0 % |

**Recomputed `estimated_cov`** (400 kb nominal, MIN 30×, MAX 300×):

| Basis | Coverage | Verdict |
|---|---|---|
| all recruited bases (as gated) | 92.60× | `ok` |
| mito-assigned read length | **27.29×** | **`low_coverage`** |
| mito merged-aligned bases only | 1.82× | `low_coverage` |

**§2.3 cross-check.** The assembly's mitochondrial contigs total
100 455 bp at 7–8×, i.e. ≈0.73 Mbp of read bases — which matches the
726 531 bp of genuinely mito-homologous sequence in the pool almost
exactly. The plastid side agrees too: 155 291 bp at 137×/168× ≈ 23.4 Mbp
against 26.1 Mbp plastid-assigned. The assembly is not missing
mitochondrial contigs; the depth was never there.

**Decision gate result:** mito-only coverage is **below** MIN under every
estimator tried, so the gate had a real false `ok` and §3.1 alone was
insufficient.

**Sibling exposure is one-directional, as §1 predicted.** Running the
same diagnostic on `INT-PLANT-01-pt` gives 98.1 % plastid / 1.9 % mito —
265.57× → 260.43×, still comfortably `ok`. `animal_mt` has no sibling
panel and is untouched. The correction changes exactly one fixture's
verdict.

### Decisions taken

Both confirmed with the user before implementation:

- **§3 option: 3.1 + 3.2.** C2 gates on target-assigned bases; C3
  additionally reports the assembled sibling fraction. 3.3
  (sibling-aware recruitment) not attempted — it needs a
  [principle 5](../../CONSTITUTION.md) sign-off and risks discarding genuine
  NUPT-spanning mt reads, and §9 item 1 is already scheduled to revisit
  recruitment filtering with real data.
- **§4 fixture: accept the `low_coverage` verdict.** `INT-PLANT-01-mt`
  now soft-fails and asserts only against the gate. The alternative —
  overriding `min_cov` in the integration profile to keep the sample
  assembling — was rejected as tuning a fixture until it passes. The
  cost is recorded honestly below.

### Implementation

**C2 — [`bin/coverage_gate.py`](../../bin/coverage_gate.py)** (spec §2.1.5):
maps the recruited pool against the declared panel and every sibling,
assigns each read to its best-scoring panel on merged aligned fraction,
and gates on `target_assigned_bases`. Ties fall to the sibling
([rule 18](../../CONSTITUTION.md)). `sample_status.json` gains
`coverage_basis`, `target_assigned_bases`, `sibling_assigned_bases`,
`target_merged_aligned_bases`, `sibling_organelle_fraction`,
`sibling_panels_scored`, and the per-class read counts.

Fallback per [principle 8](../../CONSTITUTION.md): a missing `.mmi` (either
the declared panel or a sibling) or a minimap2 failure warns to stderr,
records `coverage_basis: "total_recruited"` with
`sibling_panels_scored: []`, and estimates from the whole pool. A
degraded bundle never fails the sample.

**C3 — [`bin/bin_target.py`](../../bin/bin_target.py)** (§3.1):
`sibling_organelle_summary()` adds a `sibling_carryover` block to
`bin_metadata.json` — assembled bases, sibling bases, fraction, and a
`sibling_carryover_warning` above
`params.bin_target_thresholds.<target>.sibling_warn_fraction`
(0.30, provisional). The report requirement is recorded in
[spec §6a.2](../../spec/06a-reports.md#6a2-section-outline-our-pipelines-content)
for P4 — COLLATE and the renderer are still stubs, so there is no
report to write the warning into yet.

**Container.** C2 is no longer stdlib + `seqkit`/`seqtk`. Rebuilt via
`mulled-build build 'python=3.12,seqkit=2.13,seqtk=1.5,minimap2=2.31'`,
pushed as
`neoformit/daff-wf5-coverage-gate:python3.12_seqkit2.13_seqtk1.5_minimap2-2.31`
and pinned by digest in [`conf/containers.config`](../../conf/containers.config).
Invocation recorded in [task 2](2_containers.md) per
[spec §1a](../../spec/01-pipeline-flow.md#1a-engineering-constraints).
`COVERAGE_GATE` now takes the reference bundle as a staged `path` input
([rule 13](../../CONSTITUTION.md)) — `main.nf` passes the same
`ch_organelle_refs` value channel BIN_TARGET already used, so there is no
new reference data and no new channel.

### Deviations from the brief

1. **Estimator choice.** §3.2 says "gate on target-assigned bases only"
   without specifying whether that means the *read length* of
   target-assigned reads or only their *merged-aligned* bases. Read
   length is used: merged-aligned bases (1.82×) counts only sequence
   homologous to the single *Vigna* reference and badly under-states a
   divergent mitogenome. Both are below MIN here, so the verdict is
   robust to the choice; the stricter figure is recorded as
   `target_merged_aligned_bases` for the auditor.
2. **The declared panel missing is also a soft fallback.** §3.2's exit
   criterion names only the "sibling panel unavailable" case. C3 treats
   a missing declared panel as fatal because it cannot bin at all; C2
   can still gate, so it degrades to a whole-pool estimate rather than
   failing. Covered by its own test.
3. **`post_subsample_cov` is now scaled, not re-measured.** It was
   `total_bases(out) / nominal_size`, which under the new basis would
   report post-subsample *pool* depth against a target-only pre-value.
   It is now `est × (bases_out / bases_in)` — exact in expectation for a
   uniform `seqtk sample`, and avoids a second mapping pass.
4. **New `expected/*/coverage_bounds.json` fixtures.** §3.2 asked for
   the expected gate status to be updated. Rather than hard-code three
   statuses in `assertions.sh`, per-target bounds files were added
   matching the existing `assembly_bounds.json` / `bin_bounds.json`
   idiom, carrying expected status, expected `coverage_basis`,
   estimated-coverage bounds, and sibling-fraction bounds.
5. **The 4 pre-existing `test_coverage_gate.py` subprocess tests were
   left failing.** They fail on a missing `seqkit` in the shared test
   image — [task 27](27_unit_test_boundary_mocking.md)'s subject, not
   this task's. New coverage was added in-process with the tool boundary
   faked, which is the idiom task 27 adopts, so none of it needs
   rework. Task 27's scope is unchanged; a note was added there that its
   `_fake_run` must now also handle a `minimap2` argv.

### Knock-on cost of the §4 decision — recorded, not mitigated

`INT-PLANT-01-mt` no longer reaches METAFLYE, BANDAGE_NG, BIN_TARGET or
anything downstream. Consequences, all logged in
[`tasks/todo.md`](../todo.md):

- The `plant_mt` branch of C3 — including the sibling-organelle
  discrimination [task 23](23_bin_target_recalibration.md)
  built and the `emit: all` multi-contig path — is now covered by unit
  tests only. Task 23's `contig_1`/`contig_6` plastid regression
  assertion is withdrawn.
- [Task 26](../26_binning_marker_genes.md)'s §2 go/no-go gate loses this
  fixture as evidence; its §5 integration assertion has been struck
  through with a pointer here.
- `expected/plant_mt/{assembly,bin}_bounds.json` are retained but
  unasserted, annotated to say so. Task 23 §5.3 asked that they be
  tightened once this task resolved — they cannot be, because the
  assembly they described should not have been produced.
- A replacement fixture must clear 30× **mitochondrial** depth after the
  sibling split, not 30× recruited.

**Contributing cause worth its own item:** `plant_mt.fa` in the
reference bundle is a **single** *Vigna radiata* mitogenome (398 kb, one
sequence), against 101 plastid genomes in `plant_pt.fa` and 1853 in
`animal_mt.fa` — GetOrganelleDB's `SeedDatabase/embplant_mt.fasta` ships
one seed. A distant single reference recruits weakly on genuine
mitochondrial sequence and strongly through its own NUPT content, which
is much of why this pool came back 70 % plastid. Logged against
[spec §4.1](../../spec/04-reference-data.md#41-setup-task-source-reference-panel-from-getorganelle)'s
"re-evaluate after P1" clause.

### Verification

`-profile integration` (2026-08-03), 40 tasks succeeded / 0 failed. The
published `outdir` was cleared and republished from the work cache
before asserting — `publishDir` never deletes, so the previous run's
`INT-PLANT-01-mt` assembly artefacts would otherwise have been asserted
against. A note to that effect was added to the `assertions.sh` header.

Stage counts confirm the branch: RECRUIT and COVERAGE_GATE ran 3 of 3;
METAFLYE through ORGANELLE_MAP ran **2 of 2**; COLLATE returned to 3 of
3 for the soft-fail bundle.

| Sample | `coverage_basis` | estimated_cov | sibling frac | Status |
|---|---|---|---|---|
| `INT-ANIMAL-01` | `total_recruited` | 325.81× | — | `ok` (subsampled) |
| `INT-PLANT-01-pt` | `target_assigned` | 260.43× | 0.0193 | `ok` |
| `INT-PLANT-01-mt` | `target_assigned` | **27.29×** | **0.7053** | **`low_coverage`** |

The pipeline reproduces the §2 diagnostic exactly —
`target_assigned_bases: 10916501`, `target_merged_aligned_bases: 726531`,
`sibling_assigned_bases: 26121912`. `animal_mt` is byte-identical to its
pre-change estimate, confirming the no-sibling path is untouched.

- `bash tests/integration/assertions.sh` — **all assertions passed**,
  including the six added here (per-sample expected status, expected
  `coverage_basis`, estimated-coverage bounds, sibling-fraction bounds in
  both directions, no-assembly-on-soft-fail, and C3's `sibling_carryover`
  cross-check).
- `-profile stub -stub-run` — green; stub blocks unchanged.
- `scripts/pytest.sh scripts/tests/test_coverage_gate.py -k TestSiblingSplit`
  — 16 passed, **99 % branch coverage** of `bin/coverage_gate.py` (the
  one uncovered statement is the `__main__` guard). C2 was at **0 %
  measured** before this task.
- `scripts/pytest.sh scripts/tests/test_bin_target.py` — 54 passed, 99 %
  branch coverage of `bin/bin_target.py`.
- `flake8 bin/ scripts/` — clean at 79 columns.

**One assertion bug found and fixed while verifying:** jq's `//`
operator treats `false` as empty, so the first draft of the
`sibling_carryover` check reported a legitimately-`false`
`sibling_carryover_warning` as a missing key. It now tests with
`has()`. Worth remembering — the same trap applies to any boolean
`bin_metadata.json` field.

### Follow-up

- `tasks/todo.md` gains two items: replacing `INT-PLANT-01-mt` with a
  deeper-sequenced plant fixture (with the full list of what went dark),
  and the single-sequence `plant_mt` recruit panel.
- [Spec §9 item 1](../../spec/07-open-questions.md#9-fine-tuning-post-prototype-benchmarking)
  extended to cover `plant_mt` and to name competitive sibling-panel
  filtering at RECRUIT as the §3.3 experiment, with this task's 70.5 %
  measurement as its motivating evidence.
- The §3.1 report warning is specified in
  [spec §6a.2](../../spec/06a-reports.md#6a2-section-outline-our-pipelines-content)
  but not rendered — COLLATE and the report renderer are P4 stubs. The
  data it needs is in `bin_metadata.json` and `sample_status.json` now.
