# Task 29 — Comprehensive organellar protein panel (reference data)

**Phase:** P3/P4 boundary (from [spec §6](../spec/06-phases.md)).
**Blocks:** [task 30](30_unified_locus_pass.md),
[task 31](31_annotate_merge.md).

## 0. Overview

The reference bundle currently ships a **barcode** protein panel: five
or six genes per assembly target
([`assets/loci.json`](../assets/loci.json), realised as per-gene `.faa`
files under `refs/v2026.08/proteins/`). That is enough to extract
barcodes and nothing else.

The pipeline is about to ask more of it. Under the architecture agreed
for stages 12–13, a **single miniprot pass** against a comprehensive
organellar protein set supplies *both* the CDS annotation shown in the
report *and* the barcode sequences shipped to Taxodactyl — the barcodes
being a subset of the annotation rather than the product of a separate,
independently-run alignment. The motivation is coherence: if the report
draws `rbcL` at one set of coordinates and the barcode FASTA carries a
sequence extracted at another, those are two claims that can silently
disagree, and nothing in the pipeline would notice. Making one the
subset of the other removes the disagreement by construction
([rule 18](../CONSTITUTION.md)).

That architecture needs a protein set covering the **full protein-coding
complement** of each organelle, not just the barcode loci. Building it
is this task. It is deliberately separated from the stages that consume
it ([task 30](30_unified_locus_pass.md),
[task 31](31_annotate_merge.md)) because it lands in the versioned
reference bundle — a build, publish and checksum cycle with its own
verification — and because it carries a **go/no-go measurement** (§3)
that must pass before either consumer is worth building.

**Why now, and why this is cheap.** `MINIPROT_EXTRACT` is still a P0
stub: [`modules/local/miniprot_extract.nf`](../modules/local/miniprot_extract.nf)
`touch`es its outputs and `bin/validate_barcodes.py` (C5) does not
exist. So this is not a migration — there is no working barcode path to
put at risk, and no rework to pay for. The unified design costs nothing
extra if adopted before stage 13 is implemented, and would cost a
rewrite of the contractual output path if adopted after.

**Goal:** extend the reference bundle with a comprehensive
protein-coding gene set per assembly target, verify that the broader
panel does not degrade barcode-locus alignment, and publish the result
as a new versioned bundle.

**Exit criteria:**

- `refs/v2026.09/proteins/<target>/<GENE>.faa` covers the full
  protein-coding complement per target (§1), with representative counts
  per §2 and provenance in `manifest.json`.
- §3's panel-competition measurement run on all three integration
  fixtures and **passed** — no barcode locus lost, none degraded.
- Bundle built, checksummed, uploaded, and profile/CI versions bumped.
- `assets/loci.json` **unchanged** — see §4.
- `assets/organelle_gene_sets.json` created, carrying the canonical gene
  set per target (§5).

**Not in scope:**

- Any pipeline stage. Nothing in `modules/`, `main.nf` or `bin/`
  changes here. Consumers are tasks 30 and 31.
- tRNA and rRNA references. The unified pass is protein-coding only;
  non-CDS features come from specialist annotators
  ([task 31](31_annotate_merge.md)).
- Changing which loci are barcodes. That is `loci.json`, it is
  contractually coupled to Taxodactyl
  ([principle 9](../CONSTITUTION.md)), and it is untouched (§4).

---

## 1. What the panel becomes

Per-gene `.faa` files under `proteins/<target>/`, same layout as today —
only the gene list grows:

| Target | Today (barcode panel) | Comprehensive complement | Notes |
|---|---|---|---|
| `animal_mt` | 6 | **13** | The canonical metazoan set: `atp6 atp8 cob cox1 cox2 cox3 nad1 nad2 nad3 nad4 nad4l nad5 nad6`. Barely an expansion. |
| `plant_pt`  | 6 | **~79** | Plastid PCGs — `psa*`, `psb*`, `pet*`, `atp*`, `ndh*`, `rpo*`, `rps*`, `rpl*`, `rbcL`, `matK`, `ccsA`, `cemA`, `clpP`, `accD`, `ycf1/2/3/4`. |
| `plant_mt`  | 5 | **~40** | Plant mt PCGs — `atp*`, `ccm*`, `cob`, `cox*`, `nad*`, `mtt*`, `rps*`, `rpl*`, `matR`. |

**Source: records already in the bundle.** `validate/refseq_pt.fa` and
`validate/refseq_mt_viridiplantae.fa` were downloaded for
[§4.2](../spec/04-reference-data.md#42-refseq-organelle-dbs-for-blast_validate),
and the metazoan equivalent for `animal_mt`. Take CDS translations from
the associated RefSeq annotations rather than issuing new NCBI queries —
fewer moving parts, and the provenance is already recorded.

Extend [`scripts/refdata/build_proteins.py`](../scripts/refdata/build_proteins.py)
with a full-complement mode. Keep the existing per-gene `.faa` output
shape: [task 30](30_unified_locus_pass.md) concatenates them at runtime,
and per-gene files keep the barcode subset trivially selectable and the
bundle diff-able.

**`plant_mt` — NUPT discipline.** Build the `plant_mt` set from
**mitochondrial CDS records only**. Plastid-derived proteins leaking in
would let a NUPT insertion be annotated as a functional mitochondrial
gene — the same class of error
[task 28](completed/28_plastid_masked_mt_panel.md) masked out of the
recruitment panel. Do **not** mask the protein sequences themselves (a
mitochondrial protein is a mitochondrial protein regardless of where its
homologue also occurs); the discipline is entirely in which source
records are used. Record the resulting gene list in `manifest.json` so
the set is auditable.

## 2. Representatives per gene

5–10 divergent representatives per gene, matching the existing barcode
panel's convention and the reasoning in
[spec §9 item 11](../spec/07-open-questions.md#9-fine-tuning-post-prototype-benchmarking).

Rationale for not going broader: miniprot runtime scales with query
count, and the *marginal* sensitivity gain from the 50th representative
of a gene is small compared with the first five spanning the clade. Note
that §9 item 11 remains open — it is a sweep to run against the P3
datasets, and this task should not pre-empt it. Pick representatives to
maximise phylogenetic spread within the target's clade, not to maximise
count.

**Size check.** `plant_pt` at 79 genes × 10 representatives ≈ 790
proteins ≈ a few MB — negligible against the ~1.9 GB bundle. Record
actual sizes in Outcomes.

## 3. Go/no-go: does the broader panel change barcode calls?

**This is the gate. Run it before building anything downstream.**

Barcodes are the pipeline's contractual output
([brief.md §5](../brief.md)). Broadening the query set changes what
competes at each locus: with ~790 proteins in play instead of ~60,
miniprot may select a different best alignment at the `rbcL` locus than
it would have with an `rbcL`-only query. That is *probably* an
improvement — more evidence, better model selection — but "probably" is
not good enough for the contractual output, so measure it.

**Method.** No pipeline code is needed. The three integration fixtures
have already produced binned `target.fasta` files under
`tests/integration/output/<sample>/bin_target/`. For each, run miniprot
twice — once against the current barcode panel, once against the
comprehensive panel — and compare at the barcode loci only:

```
miniprot --gff -t 2 target.fasta <barcode panel>.faa       > narrow.gff
miniprot --gff -t 2 target.fasta <comprehensive panel>.faa > broad.gff
```

**Compare, per barcode locus in `loci.json` for that target:**

1. Is the locus still called at all?
2. Are the coordinates the same?
3. Is the translated sequence identical, and if not, is the alignment
   identity and query coverage higher or lower?

**Pass criteria — all three must hold:**

- **No barcode locus lost.** A locus called under the narrow panel and
  absent under the broad panel is a hard fail.
- **No locus degraded.** Coordinate shifts are acceptable *only* where
  identity or query coverage improves; a shift with worse metrics is a
  fail.
- **`INT-PLANT-01-mt` is exempt** — it soft-fails the coverage gate and
  has no `target.fasta`
  ([task 25](completed/25_coverage_gate_carryover.md),
  [task 28](completed/28_plastid_masked_mt_panel.md)). Note the gap;
  do not add a fixture or relax the gate for it.

**On failure:** do not proceed to tasks 30–31 as designed. The fallback
is two miniprot invocations sharing one process and one container —
comprehensive for annotation, barcode-only for extraction — which keeps
the container and reference-data work from this task but loses the
coherence guarantee. Record the measurement either way; it is the
evidence for the architecture.

**Record the full comparison table in Outcomes**, per locus per fixture.
This is the only place that evidence will exist.

## 4. `assets/loci.json` does not change

Under the unified design, the barcode panel becomes a **selector** over
the comprehensive set: a locus is a barcode where its gene symbol
appears in `loci.json` for that origin. The file keeps its schema
(`wf5/loci-panel/v1`), its content, and its contractual coupling to
Taxodactyl's accepted-loci set
([principle 9](../CONSTITUTION.md)) — arguably expressed more cleanly
than before, since the panel is now purely a statement about *which
loci are barcodes* rather than doubling as a statement about which
proteins get aligned.

Symbol matching must be **case-insensitive**: `loci.json` uses `COX1`
for `animal_mt` and `cox1` for `plant_mt`, and the comment in the file
already notes NCBI symbols are case-insensitive. Task 30 owns the
matching logic and its tests; flagged here because the two files must
agree on gene naming and this task sets the `.faa` filenames.

**One caveat to carry into task 30:** `loci.json`'s `$description`
notes the panel is protein-coding only, because rRNA barcodes (12S, 16S)
and tRNA barcodes (trnL) need a nucleotide extraction path miniprot
cannot provide. That limitation is unchanged and unaffected here.

## 5. `assets/organelle_gene_sets.json` — canonical sets

Create alongside the panel, as versioned config
([principle 9](../CONSTITUTION.md) / [rule 18](../CONSTITUTION.md)):

```json
{
  "schema": "wf5/gene-sets/v1",
  "version": "2026.08",
  "sets": {
    "animal_mt": {
      "name": "metazoan_37",
      "description": "Canonical metazoan mitogenome: 13 PCG, 22 tRNA, 2 rRNA",
      "protein_coding": ["atp6", "atp8", "cob", "cox1", "..."],
      "rrna": ["rrnL", "rrnS"],
      "trna_count": 22
    },
    "plant_pt": { "...": "~79 PCG, ~30 tRNA, 4 rRNA" },
    "plant_mt": { "...": "~40 PCG" }
  }
}
```

This is the **completeness yardstick** — what an annotation is measured
against to report "recovered 71 of 79 canonical plastid genes". It is
deliberately a separate file from `loci.json`: conflating a
completeness yardstick with the Taxodactyl-coupled barcode panel would
make the two drift together for no reason.

Distinct from the `.faa` panel too: the gene set is a *list of expected
names*, the panel is *sequences to align*. They will usually agree, but
the panel may omit a gene with no usable representative, and the
yardstick should still count it missing.

Note for [task 26](shelved/26_binning_marker_genes.md) (shelved): it
wanted a marker-gene set for `BIN_TARGET` and was blocked partly on the
absence of exactly this reference data. This removes that obstacle. It
does **not** un-shelve the task — its §2 gate also needs P3 data and a
sample with known NUMT content — but its premise should be re-read.

## 6. Bundle build and publish

New version `v2026.09` per [rule 10](../CONSTITUTION.md) immutability —
`refs/v2026.08/` is published and must not be mutated. Copy the
`recruit/`, `validate/` and `annotate/` subtrees unchanged; do not
re-derive them, which would churn every downstream fixture bound for a
protein-panel addition.

- `scripts/refdata/build_proteins.py` — full-complement mode;
  `manifest.json` entries with source URL, RefSeq release and SHA256.
- **Fold in [task 31](31_annotate_merge.md)'s MITOS2 reference data**
  (`annotate/mitos/refseq89m/`, 19.6 MB from Zenodo) if that task is
  being worked in the same cycle — one 1.2 GB publish for both
  additions rather than two. Coordinate before building; otherwise
  task 31 cuts its own version.
- `scripts/refs-v2026.09.sha256`. **Note:** `scripts/` currently carries
  only `refs-v2026.07.sha256` even though the bundle in use is v2026.08,
  so [`scripts/fetch_refs.sh`](../scripts/fetch_refs.sh) cannot verify
  the bundle CI actually fetches. Fix that gap here — write both the
  missing v2026.08 file and the new v2026.09 one.
- Upload `v2026.09/refs.tar.gz` to the `refdata-wf5` public container
  per [spec §4a](../spec/04a-azure-blob-storage.md); confirm HTTP 200 and
  byte size.
- Bump the version in [`conf/integration.config`](../conf/integration.config),
  [`conf/azure.config`](../conf/azure.config) and
  [`.github/workflows/integration.yml`](../.github/workflows/integration.yml).

## 7. Verification

1. Gene counts per target match §1; spot-check five genes per target for
   correct symbol, sensible length, and no stop codons mid-sequence.
2. `plant_mt` set contains no plastid-only gene symbols (`psbA`, `rbcL`,
   `ndh*`) — the NUPT discipline check from §1.
3. §3's measurement run and passed; comparison table recorded.
4. `bash scripts/fetch_refs.sh v2026.09` on a clean checkout — fetches,
   checksums, `proteins/` populated.
5. `nextflow run . -profile integration` still green — nothing consumes
   the new panel yet, so this is purely a regression check that the
   bundle swap broke nothing.
6. `flake8` on `scripts/refdata/build_proteins.py` via the `claude` venv.

## 8. Deliverables checklist

- [ ] `scripts/refdata/build_proteins.py` — full-complement mode.
- [ ] `refs/v2026.09/proteins/<target>/*.faa` + `manifest.json` entries.
- [ ] §3 measurement run, passed, table recorded in Outcomes.
- [ ] `assets/organelle_gene_sets.json` — all three targets.
- [ ] `scripts/refs-v2026.08.sha256` (backfill) and
      `scripts/refs-v2026.09.sha256`.
- [ ] Bundle uploaded; `conf/integration.config`, `conf/azure.config`,
      `.github/workflows/integration.yml` bumped.
- [ ] Spec: [§4.3](../spec/04-reference-data.md#43-protein-panel-and-barcode-selector)
      rewritten (panel is now comprehensive; `loci.json` is a selector),
      [§4.4](../spec/04-reference-data.md#44-consolidated-build-script)
      bundle layout updated.
- [ ] `assets/loci.json` **unchanged** (assert this in review).

## 9. Notes / non-issues

- **Why not build this inside task 30.** It lands in the versioned
  bundle, which is a build/publish/checksum cycle with its own
  verification, and §3's gate must pass before either consumer is worth
  building. Bundling them would mean discovering a §3 failure after
  writing the code that assumes it passed.
- **Why per-gene `.faa` rather than one concatenated file per target.**
  The barcode subset stays trivially selectable, the bundle stays
  diff-able, and the existing layout and build script are reused.
  Concatenation at runtime is one `cat`.
- **Runtime cost is not a concern.** ~790 proteins against a 150 kb
  target is seconds. miniprot was already being run; this changes what
  it is run *against*, not how often.
- **This does not settle the plastid annotator question.** The unified
  pass gives `plant_pt` a CDS annotation; tRNA, rRNA and intron
  structure remain open in
  [task 32](32_research_plastid_noncds.md).

## 10. Outcomes

_(fill on completion: gene counts and panel sizes per target, §3's
per-locus comparison table for each fixture, bundle URL + byte size.)_
