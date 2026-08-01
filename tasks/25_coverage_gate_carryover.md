# Task 25 — `COVERAGE_GATE` over-estimates coverage when a sibling organelle is recruited

**Phase:** P3 (from [spec §6](../spec/06-phases.md)).
**Goal:** Establish whether
[`COVERAGE_GATE`](../spec/02-stages.md#21-coverage-gate-2-stage-6) (C2)
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
([§2.1.1](../spec/02-stages.md#211-per-target-limits)). It has no way to
know which organelle the recruited bases came from.

On `INT-PLANT-01-mt` the gate reported `status: "ok"`. But
[spec §3.7](../spec/03-organelles.md#37-target-binning-criteria-per-assembly-target)
established that 155 291 bp of the 270 138 bp assembly — **57 %**, and
the two highest-coverage contigs at 137× and 168× — is chloroplast.
The genuine mitochondrial contigs assembled at **5–7×**, far under the
30× MIN the gate is supposed to enforce.

If the assembled coverage reflects the recruited pool, the gate's
estimate was inflated by plastid reads and the sample should have
soft-failed as `low_coverage`. Instead it proceeded to assembly and
produced an under-assembled ~115 kb mitogenome presented as a normal
result.

This is [principle 7](../CONSTITUTION.md) in reverse: a degraded sample
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
([§2.2 C2](../spec/02-stages.md#22-custom-logic-components): "Python
stdlib only + `seqkit`/`seqtk` in `PATH`"). Would need a container
change.

### 3.3 Sibling-aware recruitment in RECRUIT

Filter at recruitment: drop reads aligning better to a sibling panel
than to the declared one.

Fixes the root cause and improves assembly quality, not just the gate —
METAFLYE would no longer spend its depth budget on chloroplast. But it
sits closest to [principle 5](../CONSTITUTION.md), which is deliberate
about recruitment being *coarse* positive enrichment with precise
separation deferred to binning. Competitive filtering at RECRUIT is
still positive selection (no host reference involved), so it is not a
principle violation — but it is a design-principle boundary and needs
explicit sign-off per the amendment procedure. It also risks discarding
genuine mt reads spanning NUPT insertions, which are real mitochondrial
sequence that legitimately aligns better to the plastid panel.

**Preliminary lean:** 3.1 unconditionally (it is nearly free and
strictly increases auditability), plus 3.2 if §2 confirms a false `ok`.
Hold 3.3 for [§9 item 1](../spec/07-open-questions.md#9-fine-tuning-post-prototype-benchmarking)
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
can actually demonstrate. Per [rule 19](../CONSTITUTION.md), a fixture
that cannot support its assertions is a long-term maintenance liability.
Flag this outcome explicitly rather than tuning bounds until they pass.

## 5. Exit criteria

- §2 measurement complete, with the plastid/mito split of recruited
  bases and a recomputed mito-only `estimated_cov` recorded in this
  task's outcomes section.
- A decision recorded against §3, with rationale.
- If 3.1: sibling-organelle fraction in `bin_metadata.json` and a report
  warning above threshold; assertion in `assertions.sh`.
- If 3.2: C2 gates on target-assigned bases; unit tests to
  [rule 14](../CONSTITUTION.md) 100 % branch coverage including the
  "sibling panel unavailable" fallback; container change documented;
  `INT-PLANT-01-mt`'s expected gate status updated to whatever the
  corrected estimate warrants — **including `low_coverage` if that is
  the honest answer**, with `assertions.sh` updated to expect it.
- [spec §2.1](../spec/02-stages.md#21-coverage-gate-2-stage-6) updated
  with the limitation and the resolution, whichever option is taken.
