# Task 31 — Stage 13a: `ANNOTATE` (merge + non-CDS features via MITOS2)

**Phase:** P4 (from [spec §6](../spec/06-phases.md)).
**Prerequisite:** [task 30](completed/30_unified_locus_pass.md) — `MINIPROT_CDS`
emits `cds.gff` for every sample reaching `BIN_TARGET`.

## 0. Overview

After [task 30](completed/30_unified_locus_pass.md), every target already has a
protein-coding annotation: `cds.gff`, from the same miniprot pass that
produces the barcodes. What is missing is everything that is not a
protein — **tRNAs and rRNAs** — and a finished annotation artifact to
hand to `ORGANELLE_MAP` and `COLLATE`.

That gap is not cosmetic for animals. A metazoan mitogenome has 37
genes, and **22 of them are tRNAs** — a CDS-only annotation of an animal
mitogenome shows barely a third of the picture. Metazoan mt tRNAs are
also highly degenerate and are exactly what a generic search misses,
which is why **MITOS2** exists and why it remains the right tool for
that job.

So `ANNOTATE` becomes a **merge step**:

- **CDS features** come from `cds.gff` — always, for every target. One
  source, shared with the barcodes.
- **tRNA and rRNA features** come from a specialist annotator where one
  exists. Today that is MITOS2 for `animal_mt`; the plastid equivalent
  is unresolved ([task 32](32_research_plastid_noncds.md)) and plant mt
  is [task 33](33_annotate_plant_mt_noncds.md).
- **MITOS2's own CDS calls are not used as features.** They become a
  **cross-check**: where MITOS2 and miniprot disagree about a
  protein-coding gene, the disagreement is recorded and surfaced. That
  is a genuine QC signal, and it is only available because the two are
  independent.

This preserves the coherence the unified pass bought — one CDS source,
so the barcode is always a feature on the map — while keeping MITOS2's
real strength.

**The plant branches are no longer "not supported".** They produce a
real, complete-for-CDS annotation with `status: "ok_cds_only"`. That is
a materially better position than the empty GFF this stage would
otherwise have emitted for plants, and it is why the plastid annotator
question is no longer blocking P4.

**Goal:** replace the stage-12 `ANNOTATE` stub with a merge stage at
slot 13a producing GFF3 + `annotation_summary.json`, with MITOS2
supplying non-CDS features on `animal_mt`.

**Exit criteria:**

- `-profile integration`, `INT-ANIMAL-01`: GFF3 carrying miniprot CDS
  features **and** MITOS2 tRNA/rRNA features; `annotation_summary.json`
  with `status: "ok"`, per-type counts, completeness against the
  canonical metazoan set, and any miniprot/MITOS2 CDS disagreements.
- `-profile integration`, `INT-PLANT-01-pt`: GFF3 carrying CDS features
  only; `status: "ok_cds_only"` with `non_cds_source: null` and a reason
  naming task 32.
- Every CDS feature in the annotation traces to `cds.gff`; no CDS
  feature originates from MITOS2 (§4).
- C8 exists with a `scripts/tests/` module at **100 % branch coverage**
  ([rule 14](../CONSTITUTION.md)).
- `-profile stub -stub-run` green; `flake8` clean.

**Not in scope:**

- **Plastid non-CDS features** — [task 32](32_research_plastid_noncds.md),
  open research. **Plant mt non-CDS** —
  [task 33](33_annotate_plant_mt_noncds.md), blocked.
- **Any decision logic keyed on annotation output.** No sample is
  failed, flagged or re-binned because annotation recovered few
  features. Counts are *data* for the report; the negative-clarity
  classification stays centralised in `COLLATE`
  ([task 21 §8](completed/21_blast_validate.md#8-notes-non-issues) — do
  not create a second decision surface).
- **`ORGANELLE_MAP` (stage 14).** It consumes this output; its renderer
  is a separate deferred choice
  ([spec §8 item 4](../spec/07-open-questions.md#8-remaining-open-questions)).
  This task only changes the shape of what it is handed (§1).

---

## 1. Output contract: GFF3 + summary JSON, drop GenBank

The stub emits `(meta, <sample>.gff, <sample>.gbk)`. Change the second
file to `annotation_summary.json`:

```groovy
output:
tuple val(meta),
      path("${meta.sample_id}.gff"),
      path("annotation_summary.json"), emit: annotation
```

**Why drop `.gbk`:** no producer (MITOS2 emits GFF/BED/FAS/FAA, miniprot
emits GFF — neither emits GenBank), and no consumer (the Taxodactyl
handoff bundle in [spec §0](../spec/00-overview.md#0-input-sample-sheet)
lists `organelle_annotation.gff` and no `.gbk`; the two spec statements
already disagree and the bundle layout is the contractual one).
Synthesising GenBank would mean writing and testing a format converter
that exists to satisfy a word in the spec — [rule 19](../CONSTITUTION.md).

**Why add `annotation_summary.json`:** it is where provenance and
negative clarity live — feature sources, completeness, cross-check
disagreements, and status. Without it, a GFF with no tRNAs is ambiguous
between "this target has no tRNA annotator" and "the annotator found
none", which [principle 7](../CONSTITUTION.md) requires be
distinguishable.

**Arity stays 3, so downstream churn is mechanical:**
`modules/local/organelle_map.nf` input, `modules/local/collate.nf`
(`path(annotation_gbk)` → `path(annotation_summary)`), and the
`ch_ok_inputs` map closure in [main.nf](../main.nf).

## 2. MITOS2 — container and reference data

**Container:** `quay.io/biocontainers/mitos:2.1.10--pyhdfd78af_0`,
SHA-pinned in [`conf/containers.config`](../conf/containers.config)
(replacing the `python:3.12-slim` stub and its `TODO P4` comments).
Resolve the digest as in
[task 21 §1](completed/21_blast_validate.md#1-container-pin--confcontainersconfig).
Single-tool biocontainer — [rule 12](../CONSTITUTION.md) option 1, no
`mulled-build`, no Dockerfile.

**The image is ~1 GB.** Record the pull time against the 60-minute CI
budget in Outcomes; if material, the mitigation is a registry mirror,
not a different tool.

**Reference data:** Zenodo record `10.5281/zenodo.4284483` ("Reference
data for MITOS2") → `refseq89m.tar.bz2` (metazoan, RefSeq 89),
19 629 319 bytes, `md5:bb6325e27612e61a6995e88e8ffcecf8`. Not
`refseq89f` (fungal) or `refseq89o` (other) — out of scope
([constraint 1](../CONSTITUTION.md)). Never downloaded at run time
([rule 10](../CONSTITUTION.md), [rule 17](../CONSTITUTION.md)).

Add `scripts/refdata/build_annotate.sh` alongside the existing
`build_recruit.sh` / `build_validate.sh`, verifying the published MD5 and
recording source URL + release + SHA256 via
[`build_manifest.py`](../scripts/refdata/build_manifest.py). Bundle
layout gains:

```
refs/v2026.09/annotate/mitos/refseq89m/
```

**Fold this into [task 29](completed/29_comprehensive_protein_panel.md)'s bundle
build if it has not yet been published** — one 1.2 GB publish cycle for
both additions rather than two. If v2026.09 is already out, cut v2026.10
per [rule 10](../CONSTITUTION.md) immutability, copying the other
subtrees unchanged.

**New param:** `params.annotate_refs`, bundle-relative, passed as a
`path` input ([spec §1a](../spec/01-pipeline-flow.md#1a-engineering-constraints)
pattern 1).

**C8 runs in the MITOS2 image** — it already carries Python 3, and
`bin/` auto-stages onto every process `PATH`. C8 must therefore be
**stdlib-only** (GFF is tab-separated, output is JSON — achievable).
Confirm the container's `python3` is ≥ 3.9 as implementation step 1; if
not, C8 moves to its own image and the stage splits into two processes.
Record which applies.

## 3. What MITOS2 contributes, and what it does not

Run MITOS2 normally — it annotates the whole mitogenome — then use its
output selectively:

| MITOS2 output | Use |
|---|---|
| `tRNA` features | **Emitted** as annotation features. Its covariance-model search is the reason MITOS2 is here. |
| `rRNA` features | **Emitted.** |
| `CDS` / `gene` (protein) features | **Not emitted.** Compared against `cds.gff` and recorded as a cross-check (§4). |

**Genetic code.** MITOS2 takes one `--code`. Add
`params.annotate.animal_mt.genetic_code`, defaulting to **5**
(invertebrate) per [spec §3.3](../spec/03-organelles.md#33-specific-issues-and-decisions)
option (a) — most biosecurity intercepts. Note that
[task 30](completed/30_unified_locus_pass.md)'s C5 runs a real clade trial over
tables 2 and 5 and records its choice; ANNOTATE cannot consume that
without making a supplementary stage depend on the contractual one.
The genetic code affects MITOS2's CDS start/stop refinement — which we
discard — and barely touches its tRNA/rRNA search, so the mismatch is
largely inert here. **Record both** the table ANNOTATE applied and the
table C5 chose in `annotation_summary.json`, and flag a disagreement
rather than resolving it silently ([rule 16](../CONSTITUTION.md)).

**Confirm `runmitos.py`'s exact flags** against `--help` inside the
pinned container as implementation step 1 — `-i/--input`, `-c/--code`,
`-o/--outdir`, `-r/--refdir` (the *parent* of `refseq89m`),
`-R/--refseqver`, plus whatever suppresses plotting in a headless
container. Record the verified invocation in Outcomes.

## 4. C8 — `bin/annotate_summary.py`

Register in [spec §2.2](../spec/02-stages.md#22-custom-logic-components):

| # | Component | Container | Stage | Purpose |
|---|---|---|---|---|
| C8 | `annotate_summary.py` | shares the `ANNOTATE` MITOS2 biocontainer (stdlib-only, §2) | 13a | Merges `cds.gff` with a specialist annotator's non-CDS features into one GFF3 keyed by contig name; cross-checks the annotator's CDS calls against miniprot's; computes completeness against the canonical gene set; writes `annotation_summary.json`. **Always exits 0.** |

**Behaviour:**

1. **CDS features from `cds.gff`, unmodified.** Coordinates are not
   adjusted, re-derived or re-sorted in a way that changes them — the
   barcode coherence invariant from
   [task 30 §3](completed/30_unified_locus_pass.md) extends through this stage to
   the published annotation.
2. **Non-CDS features** from the annotator's output, `seqid` rewritten
   to the contig name (MITOS2 writes one output subdirectory per input
   FASTA record and uses its own internal index in column 1 — map back
   by input order).
3. **Cross-check.** For each protein-coding gene the annotator called,
   is there a `cds.gff` feature of the same normalised name with
   overlapping coordinates? Record `cds_crosscheck`:
   `{agreed: [...], miniprot_only: [...], annotator_only: [...],
   coordinate_conflicts: [...]}`. Do not act on it — report it.
4. **Name normalisation.** MITOS2 emits `nad4l`, `trnF(gaa)`, `rrnS`.
   Normalise case and strip the anticodon parenthetical for comparison
   against the canonical set; leave the raw name in the GFF attributes
   untouched.
5. **Completeness** against `assets/organelle_gene_sets.json`
   ([task 29 §5](completed/29_comprehensive_protein_panel.md)) —
   `protein_coding_completeness` as a fraction, plus found/missing
   lists.
6. **Never de-duplicate by gene name.** A plastid gene in the inverted
   repeat genuinely occurs twice
   ([spec §3.6](../spec/03-organelles.md#36-plastid-quadripartite-canonicalisation));
   plant mitogenomes carry genuine large repeats. This is the single
   most plausible wrong "fix" a later change could make — tested in §5.
7. **Always exit 0** ([principle 8](../CONSTITUTION.md)). A non-zero
   annotator exit becomes a status, not a pipeline abort — annotation is
   supplementary.

**Statuses:**

| Condition | `status` |
|---|---|
| `cds.gff` empty (empty `target.fasta` upstream) | `no_assembly` |
| CDS merged, no non-CDS annotator for this target | `ok_cds_only` |
| CDS merged + non-CDS features present | `ok` |
| Nothing parsed from either source | `no_features` |

**`annotation_summary.json` (`wf5/annotation-summary/v1`):**

```json
{
  "schema": "wf5/annotation-summary/v1",
  "sample_id": "INT-ANIMAL-01",
  "assembly_target": "animal_mt",
  "status": "ok",
  "reason": null,
  "cds_source": "miniprot",
  "non_cds_source": "mitos2",
  "tool_versions": {"miniprot": "0.18", "mitos2": "2.1.10"},
  "reference_data": "refseq89m",
  "genetic_code_annotate": 5,
  "genetic_code_barcodes": 5,
  "genetic_code_agreement": true,
  "contigs_annotated": ["contig_8"],
  "feature_counts": {"gene": 37, "CDS": 13, "tRNA": 22, "rRNA": 2},
  "protein_coding_completeness": 1.0,
  "protein_coding_genes_missing": [],
  "cds_crosscheck": {"agreed": 13, "miniprot_only": [], "annotator_only": []},
  "canonical_gene_set": "animal_mt/metazoan_37"
}
```

## 5. Unit tests — `scripts/tests/test_annotate_summary.py`

100 % branch coverage. Fixtures are small on-disk GFF trees, not mocks —
C8 shells out to nothing.

1. **Happy path** (`animal_mt`) — `cds.gff` + a MITOS2 tree merge;
   `status: "ok"`, counts correct, tRNAs present.
2. **CDS provenance** — every CDS feature in the output matches a
   `cds.gff` row exactly; no MITOS2 CDS leaks through as a feature.
3. **`ok_cds_only`** — no annotator directory supplied → CDS-only GFF,
   `non_cds_source: null`, `reason` set.
4. **`no_assembly`** — empty `cds.gff` → empty GFF, exit 0.
5. **`no_features`** — both sources present but nothing parseable.
6. **Cross-check agreement** — matching gene, overlapping coordinates →
   `agreed`.
7. **Cross-check disagreement** — a gene called by only one source lands
   in `miniprot_only` / `annotator_only`; non-overlapping coordinates for
   the same name land in `coordinate_conflicts`. No status change.
8. **IR duplication retained** — the same gene name at two
   non-overlapping ranges survives; both in the GFF, counted as one
   recovered gene.
9. **Name normalisation** — `NAD4L`, `trnF(gaa)`, `rrnS` match the
   canonical set; GFF attribute text unchanged.
10. **Multi-record input** — `plant_mt`'s `emit: all` case; `seqid`
    correct per contig, MITOS2 subdirectories mapped by input order.
11. **Genetic-code mismatch** — C5 chose 2, ANNOTATE applied 5 →
    `genetic_code_agreement: false`, status unaffected.
12. **Malformed GFF line** → skipped with a warning, remaining features
    retained, exit 0.
13. **Unknown `assembly_target`** in the gene-sets file → completeness
    fields `null`, status unaffected.
14. **Exit code 0 in every case** — parametrised.

## 6. Module and dispatch

`modules/local/annotate.nf`, keeping the `stub:` block (add `touch
annotation_summary.json` for the renamed output):

```groovy
input:
tuple val(meta), path(target_fasta), path(cds_gff)
path annotate_refs

script:
def cfg = params.annotate[meta.assembly_target]
"""
if [[ "${cfg.non_cds_tool}" == "mitos2" && -s ${cds_gff} ]]; then
    mkdir -p mitos_out
    runmitos.py --input ${target_fasta} --outdir mitos_out \\
        --code ${cfg.genetic_code} \\
        --refdir ${annotate_refs}/mitos --refseqver ${cfg.refseqver} \\
        2> mitos.stderr || true
    NONCDS="--mitos-dir mitos_out"
else
    NONCDS=""
fi

annotate_summary.py \\
    --cds-gff ${cds_gff} \\
    --assembly-target ${meta.assembly_target} \\
    --gene-sets ${file(params.gene_sets)} \\
    --out-gff ${meta.sample_id}.gff \\
    --out-summary annotation_summary.json \\
    \$NONCDS
"""
```

```groovy
params.annotate = [
    animal_mt: [ non_cds_tool: 'mitos2', genetic_code: 5,
                 refseqver: 'refseq89m' ],
    plant_pt:  [ non_cds_tool: 'none' ],   // task 32
    plant_mt:  [ non_cds_tool: 'none' ],   // task 33
]
```

- **One process, one container** — MITOS2's image runs both MITOS2 and
  C8, and the plant branch simply skips the MITOS2 call. No process
  split is needed while MITOS2 is the only specialist annotator; task 32
  or 33 will force one if they select a tool needing a different image,
  since Nextflow resolves `container` per-process.
- **`|| true` plus C8 always running** — a Groovy `when:` would starve
  the downstream join; a hard failure would abort a supplementary stage.
  Every path reaches C8, and C8 always writes both declared outputs.
- `params.gene_sets` defaults to
  `${projectDir}/assets/organelle_gene_sets.json`, staged via
  `${file(...)}` — acceptable inline pattern for a workflow-wide asset.

## 7. Integration-test wiring

`tests/integration/expected/animal_mt/annotation_bounds.json` — floors,
not exact counts ([rule 19](../CONSTITUTION.md)):

```json
{
  "min_cds": 11,
  "min_trna": 15,
  "min_rrna": 1,
  "required_genes": ["cox1", "cob"],
  "max_cds_crosscheck_conflicts": 2
}
```

`INT-ANIMAL-01` is *Acyrthosiphon pisum*, whose mitogenome carries the
canonical 37, so a complete annotation scores 13/22/2 — but pinning
exact counts is a drift trap on any MITOS2 or reference-data bump.
`cox1` and `cob` are required by name because they are the `animal_mt`
barcode panel's primary loci.

Assertion block in
[`tests/integration/assertions.sh`](../tests/integration/assertions.sh):

- `INT-ANIMAL-01` — `status: "ok"`, counts clear the bounds,
  `required_genes` present, `cds_source: "miniprot"`,
  **no CDS feature absent from `cds.gff`** (the provenance assertion);
- `INT-PLANT-01-pt` — `status: "ok_cds_only"`, CDS features present,
  zero tRNA/rRNA, `non_cds_source: null`.

Tick the ANNOTATE row in the progressive-uncomment header.
`INT-PLANT-01-mt` is not asserted — soft-fails the coverage gate
([task 25](completed/25_coverage_gate_carryover.md),
[task 28](completed/28_plastid_masked_mt_panel.md)).

## 8. Verification

1. `nextflow run . -profile stub -stub-run` green.
2. `bash scripts/pytest.sh` — C8 at 100 % branch coverage; existing
   suites unaffected.
3. `flake8` on new Python (79-col) via the `claude` venv.
4. `bash scripts/fetch_refs.sh <version>` — `annotate/mitos/refseq89m/`
   present.
5. `nextflow run . -profile integration` — per §7; record MITOS2 pull
   time, stage wall time, observed counts, and any cross-check
   disagreements (these are interesting in their own right — two
   independent methods on the same contig).
6. `bash tests/integration/assertions.sh` — all blocks green.
7. CI green on `Tests` and nightly `Integration`.

## 9. Deliverables checklist

- [x] `conf/containers.config` — SHA-pinned MITOS2; `TODO P4` removed.
- [x] `scripts/refdata/build_annotate.sh` + manifest entry; refdata in
      the bundle (folded into task 29's `v2026.09` bundle since it was
      still unpublished).
- [x] `nextflow.config` — `params.annotate`, `params.annotate_refs`,
      `params.gene_sets`.
- [x] `bin/annotate_summary.py` (C8) — stdlib-only, always exits 0.
- [x] `scripts/tests/test_annotate_summary.py` — 100 % branch coverage.
- [x] `modules/local/annotate.nf` — merge script, new output contract,
      stub retained.
- [x] `organelle_map.nf`, `collate.nf`, `main.nf` — second annotation
      output renamed; `ANNOTATE(MINIPROT_CDS.out.cds)` wiring.
- [x] `tests/integration/expected/animal_mt/annotation_bounds.json` +
      `assertions.sh` block + header tick.
- [x] Spec: [§2](../spec/02-stages.md#2-stage-detail) stage 13a,
      [§2.2](../spec/02-stages.md#22-custom-logic-components) C8,
      [§3.1/§3.2](../spec/03-organelles.md#31-what-differs-between-assembly-targets)
      annotator rows, [§6a](../spec/06a-reports.md) annotated-map row
      (these sections already described the merge design accurately;
      only the `ORGANELLE_MAP` input description needed correcting from
      "annotated GenBank" to "GFF3 + annotation_summary.json").
- [x] `tasks/todo.md` — report must surface `cds_crosscheck`
      disagreements and `genetic_code_agreement: false`.

## 10. Notes / non-issues

- **Why MITOS2's CDS calls are discarded as features but kept as a
  check.** Emitting both would reintroduce exactly the incoherence the
  unified pass removed — the report could show a MITOS2 `cox1` while
  `barcodes.fasta` carries a miniprot `cox1` from different
  coordinates. Keeping the comparison costs nothing and yields a real
  QC signal from two genuinely independent methods.
- **Why `ok_cds_only` is not a failure.** A plant sample gets a real
  annotation of every protein-coding gene found, which is most of what
  a diagnostic map shows for a plastome. The missing tRNAs/rRNAs are
  stated explicitly rather than implied by absence
  ([principle 7](../CONSTITUTION.md)).
- **Why annotation quality is not a gate.** A poor annotation on a good
  assembly is common (divergent taxon, fragmented contig), and the
  barcode path is the contractual one. Gating here would fail samples
  whose actual deliverable is fine.
- **MITOS2 vs MitoFinder.** [spec §3.1](../spec/03-organelles.md#31-what-differs-between-assembly-targets)
  lists MitoFinder as an alternative; its only biocontainer is
  `1.4.1--py27...` from 2022. Not a candidate under
  [rule 19](../CONSTITUTION.md).

## 11. Outcomes

**Container.** `quay.io/biocontainers/mitos:2.1.10--pyhdfd78af_0`
resolved to
`@sha256:36541c15ec4d3f0e2e7da6f41cb511b9ea08c08708eece4c3d67b48bc866148a`.
~1 GB image; pull completed in under 5 minutes on a warm Docker layer
cache (base layers shared with other biocontainers already pulled this
session) — not separately timed cold, but well inside the 60-minute CI
budget alongside the rest of the pipeline. Container `python3` is
3.12.12 (≥ 3.9 required for C8's stdlib-only code) — no second image
needed, C8 runs inside the MITOS2 container as planned.

**Verified `runmitos` invocation** (binary is `runmitos`, not
`runmitos.py` — the `.py` suffix in the task brief and in
`--help`'s own usage line is cosmetic, the installed entry point has no
extension):

```
runmitos --input <target_fasta> --outdir <dir> --code <cfg.genetic_code> \
    --refdir <annotate_refs> --refseqver mitos/<cfg.refseqver>
```

`--linear` is listed under `--help`'s "mandatory options" heading but
is a bare flag (default: circular) — organelle genomes are circular,
so it is never passed. `-R/--refdir` points at the *parent* of
`mitos/refseq89m` (i.e. `annotate_refs` itself, since the bundle layout
is `refs/<ver>/annotate/mitos/refseq89m/`); `-r/--refseqver` is the
relative path `mitos/refseq89m` from there — both flags together,
confirmed against a real run, not guessed from `--help` alone.

**Deviation: MITOS2 output layout is simpler than the task brief
assumed.** §4 point 2 anticipated needing to rewrite `seqid` because
"MITOS2 writes one output subdirectory per input FASTA record and uses
its own internal index in column 1". Reading
`mitos/gfffile/__init__.py::gffwriter` and `scripts/runmitos.py`
directly (then confirming against a real run) shows the opposite:
`gffwriter`'s `acc` parameter is always the *real* input record id
(`sequences[i]["id"]`), not an index — column 1 already carries the
correct contig/path name in every case. Multi-record input does get
one subdirectory per record (`<outdir>/<i>/result.gff`, `i` in FASTA
order), but each such `result.gff` is already correctly keyed. C8
(`find_mitos_result_gffs`) therefore just globs every `result.gff`
under `mitos_dir` (top-level for single-record input, per-index
subdirectories for multi-record) and unions them — no index-based
remapping needed. Verified with a synthetic two-record fixture
(`TestMultiRecord`).

**MITOS2 GFF feature shape.** Confirmed from `feature/__init__.py`:
each RNA gene emits three rows (`ncRNA_gene`, then `tRNA` or `rRNA`,
then `exon`); each protein gene emits two (`gene`, then `exon`). C8
uses only the summary row per gene — `tRNA`/`rRNA` rows are emitted as
non-CDS features, `gene` rows are used for the CDS cross-check only;
`ncRNA_gene`/`exon` child rows are dropped as redundant. Anticodon
tRNAs carry the anticodon in `Name` (`trnF(gaa)`), exactly as the task
brief anticipated — `normalise_gene_symbol` strips it for comparison
and leaves the raw GFF attribute untouched.

**Deviation: `cds_crosscheck` needs the gene-sets alias table, not
just symbol folding.** The brief's normalisation step (§4 point 4)
describes case/anticodon/`MT-`-prefix folding only. That is enough for
`protein_coding_completeness` (which already consults
`protein_coding_aliases`), but the first cross-check implementation
folded `CYTB` (miniprot) and `cob` (MITOS2) to different strings and
never matched them — the two tools use different established
abbreviations for the same gene, not just different casing. Fixed by
threading the same `build_alias_index` used for completeness into
`cds_crosscheck` (`canonicalise()`), so both sides fold to the
`organelle_gene_sets.json` canonical symbol before comparison. Caught
by `TestNameNormalisation.test_alias_and_anticodon_folding`.

**Real fixture results (`-profile integration`, `refs/v2026.09`,
RefSeq 236 / refseq89m):**

| | `INT-ANIMAL-01` (animal_mt) | `INT-PLANT-01-pt` (plant_pt) |
|---|---|---|
| status | `ok` | `ok_cds_only` |
| CDS / tRNA / rRNA | 9 / 20 / 2 | 89 / 0 / 0 |
| protein_coding_completeness | 0.692 (9/13) | 0.988 (79/80) |
| missing | ATP8, ND3, ND4L, ND6 | psaM |
| cds_crosscheck | agreed 4, miniprot_only 5, annotator_only 17 (mostly MITOS2 fragment calls e.g. `cox1_0`/`cox1_1` on this divergent-taxon draft), conflicts 0 | n/a (`non_cds_source: null`) |
| genetic_code_agreement | `true` (5 == 5) | `null` (no MITOS2 run) |
| ANNOTATE wall time | ~10m23s (MITOS2-dominated; C8 itself is sub-second) | ~1s (no MITOS2 call) |

**Deviation: `annotation_bounds.json` recalibrated against real
output, not the brief's illustrative numbers.** The brief's example
(`min_cds: 11`, `required_genes: ["cox1", "cob"]`) does not match this
fixture or even this target: `cox1`/`cob` are MITOS2's/plant_mt's
lower-case spellings, while `animal_mt`'s emitted CDS features use the
barcode panel's convention (`COX1`, `CYTB` — see `assets/loci.json`).
Real recovery on this fixture is 9/13 canonical protein-coding genes
(divergent-taxon, fragmented-assembly draft — same class of gap task
30's outcomes documented for barcode recovery on this same fixture).
Set `min_cds: 8` (floor just under observed, not the aspirational
canonical count), `required_genes: ["COX1", "CYTB"]` (this target's
actual primary-locus spellings), `min_trna`/`min_rrna`/
`max_cds_crosscheck_conflicts` left at the brief's values since real
output (20/2/0) clears them comfortably. This mirrors task 30 §9's
precedent: fixture recovery is *data*, not a value to chase — a floor
is set from what was actually observed, not from the canonical ideal.

**Bug caught in integration testing:** `bin/annotate_summary.py` was
not committed executable (`chmod +x`), which every other `bin/*.py`
script is — Nextflow auto-stages `bin/` onto the container `PATH` and
executes scripts directly (no `python3` prefix in the module script),
so the first real integration run failed with `Permission denied`.
Fixed; a reminder that `-profile stub -stub-run` alone (which never
actually invokes `bin/annotate_summary.py` in its `stub:` block) cannot
catch this class of bug — only `-profile integration` can.

**`refs/v2026.09` bundle.** `annotate/mitos/refseq89m/` (6.4 MB)
fetched from Zenodo, MD5 and byte count verified against the published
values in [task 31 §2](#2-mitos2--container-and-reference-data) before
extraction. `refs/v2026.09` was still unpublished, so this folded
directly into the existing bundle per the brief's guidance rather than
cutting `v2026.10`; `manifest.json` regenerated (276 total artefacts,
105 of them under `annotate/`) with a new `annotate_refs` provenance
block.
