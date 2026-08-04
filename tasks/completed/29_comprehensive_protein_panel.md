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

- [x] `scripts/refdata/build_proteins.py` — full-complement mode.
- [x] `refs/v2026.09/proteins/<target>/*.faa` + `manifest.json` entries.
- [x] §3 measurement run, passed (with caveat), table recorded in Outcomes.
- [x] `assets/organelle_gene_sets.json` — all three targets.
- [x] `scripts/refs-v2026.08.sha256` (backfill) and
      `scripts/refs-v2026.09.sha256`.
- [x] Bundle uploaded; `conf/integration.config`,
      `.github/workflows/integration.yml` bumped. (`conf/azure.config`
      carries no refs-version reference — see Outcomes.)
- [x] Spec: [§4.3](../spec/04-reference-data.md#43-protein-panel-and-barcode-selector)
      rewritten (panel is now comprehensive; `loci.json` is a selector),
      [§4.4](../spec/04-reference-data.md#44-consolidated-build-script)
      bundle layout updated.
- [x] `assets/loci.json` **unchanged** (asserted in review — `git diff`
      shows no changes to that file).

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

**Deviation from the brief: source method.** §1 anticipated taking CDS
translations from `validate/refseq_pt.fa` /
`validate/refseq_mt_viridiplantae.fa` (genomic FASTA already in the
bundle). Those files carry no annotation, so there is nothing to take a
CDS translation *from*. Downloaded the matching RefSeq **GenBank flat
files** instead (`plastid.{1,2,3}.genomic.gbff.gz`,
`mitochondrion.1.genomic.gbff.gz`, RefSeq release 236 — same release
already pinned for `validate/`) and parsed CDS `/gene` +
`/translation` qualifiers with biopython, inside the
`neoformit/daff-wf5-scripts:test` container (per user steer mid-task —
keeps the biopython dependency out of the host and off
`scripts/refdata/build_proteins.py`, which stays stdlib-only).
New script: `scripts/refdata/parse_gbff_cds.py`. Kingdom-splits the
mitochondrion GBFF (not pre-split, unlike plastid) using each record's
own GenBank taxonomy lineage rather than a fresh accession→taxid
lookup — simpler than reusing `split_refseq_mt.py`'s accession
allowlist and avoids a second NCBI round trip.

**Canonical gene lists.** §1's gene-family lists (`psa*`, `atp*`, …)
are prefixes, not exact symbol lists — building `organelle_gene_sets.json`
required literal per-gene canonical spellings. Used the standard
metazoan-mitogenome (13 PCG), angiosperm-plastome (~79 PCG), and
angiosperm-mitogenome (~41 PCG) gene complements as documented in the
plant/animal organelle genomics literature, with a
`protein_coding_aliases` table per target for `/gene` qualifier
spelling variants seen in RefSeq (`COX1`/`COI`, `CYTB`/`COB`,
`ND1`/`NAD1`, an `MT-` locus prefix on vertebrate records, `psbZ`/`ycf9`,
`mttB`/`orfB`). Final counts: **animal_mt 13, plant_pt 80, plant_mt
41** — all within or essentially matching the brief's ~13/~79/~40.
`.faa` filenames match `loci.json`'s existing casing convention per
target (`COX1` uppercase for `animal_mt`, `cox1` lowercase for
`plant_mt`, mixed-case gene symbols for `plant_pt`) so the two files
"agree on gene naming" per §4 without needing more than case-insensitive
matching in task 30.

**Bug caught by the build:** representative selection picks up to 10
genomes per gene independently *per gene*, so the same source genome
is frequently the best representative for more than one gene (it is
diversity-*per-gene*, not diversity-*per-genome*). The first build used
the bare accession as the FASTA record ID, which collided once
per-gene files are concatenated into one query at build/runtime —
10/130 IDs collided for `animal_mt`, 12/800 for `plant_pt`, 20/410 for
`plant_mt`. Fixed by making the record ID `<accession>_<gene>`
(`write_faa` in `build_proteins.py`) before doing anything downstream
with the panel. This matters for task 30 too, since `MINIPROT_CDS`
concatenates the same per-gene files at runtime.

**Panel sizes.** `refs/v2026.09/proteins/` is 772 KB total (134 `.faa`
files, all genes reaching the full 10 representatives — none hit the
5-representative floor).

**§3 go/no-go measurement — pass, with caveat.** Ran miniprot narrow
(current `v2026.08` barcode panel, up to 50 reps/gene) vs. broad
(`v2026.09` comprehensive panel, 5–10 reps/gene) against the two
available fixtures' `bin_target/target.fasta`
(`INT-PLANT-01-mt` exempt — no `target.fasta`, per §3). Metric:
best-hit identity and query-coverage fraction per barcode locus.

| target | gene | narrow id / qcov | broad id / qcov | verdict |
|---|---|---|---|---|
| animal_mt | COX1 | 0.671 / 0.992 | 0.663 / 0.990 | same coords |
| animal_mt | COX2 | 0.570 / 0.938 | 0.562 / 0.925 | **degraded (both worse)** |
| animal_mt | COX3 | 0.490 / 0.913 | 0.494 / 0.850 | mixed |
| animal_mt | CYTB | 0.602 / 0.935 | 0.583 / 0.963 | mixed |
| animal_mt | ND1  | 0.484 / 0.920 | 0.460 / 0.987 | mixed |
| animal_mt | ATP6 | 0.480 / 0.987 | 0.453 / 0.987 | same coords |
| plant_pt | rbcL | 0.956 / 0.998 | 0.964 / 0.994 | mixed |
| plant_pt | matK | 0.776 / 0.992 | 0.769 / 0.998 | mixed |
| plant_pt | atpB | 0.968 / 0.998 | 0.972 / 0.998 | same coords |
| plant_pt | ndhF | 0.834 / 0.999 | 0.846 / 0.996 | mixed |
| plant_pt | psbA | 0.994 / 0.997 | 0.994 / 0.997 | same coords |
| plant_pt | rpoB | 0.951 / 0.999 | 0.949 / 0.999 | same coords |

No locus lost. All coordinate shifts are within a few bases. Under a
literal reading of "a shift with worse metrics is a fail," `COX2` fails
(both identity and coverage down ~1 point) and 6/12 loci are "mixed"
(one metric fractionally up, one fractionally down, all within 2
percentage points). This measurement is **confounded**: the narrow
panel was built at up to 50 representatives/gene, the broad panel at
the spec's 5–10, so some of the spread reflects representative *count*
rather than the panel-*breadth* competition effect §3 is meant to
isolate. **Decision (user, mid-task): pass with this caveat** — no
locus lost, margins are small and plausibly representative-count
noise rather than systematic panel-breadth degradation — proceed to
build the `v2026.09` bundle. A cleaner isolation (narrow panel rebuilt
at max 10 reps/gene to match) is one candidate follow-up if task 30/31
integration surfaces a real barcode regression; not done here since the
gate passed on the softer reading.

**Bundle.** `refs/v2026.09/` = `recruit/` and `validate/` copied
unchanged from `v2026.08` + the new `proteins/` above.
`scripts/refdata/build_manifest.py` extended to fold a new
`proteins/provenance.json` sidecar (per-gene candidate/representative
counts, chosen accessions, RefSeq release) into `manifest.json` under
`protein_panel`, alongside the existing `recruit_panels` /
`refseq` sections.

- `refs-v2026.08.sha256` backfilled by downloading the already-published
  `v2026.08` tarball and hashing it directly (no local rebuild — avoids
  any doubt about byte-for-byte match with what CI actually fetches).
  `3b5d7ee6cfd26353b4e28a87c39fb14e9e45a1476bb1e7c180831c3e720fa0f6`
- `refs-v2026.09.sha256`:
  `5719b0dfa6fef4f1f250f47668db1fe928658309477ca5b0434912ceadbf7a61`
- Uploaded to `https://daffstandard.blob.core.windows.net/refdata-wf5/v2026.09/refs.tar.gz`
  — 1,213,124,948 bytes, HTTP 200 confirmed, byte size matches the local
  tarball exactly.
- `conf/azure.config` carries no refs-version reference (it only
  configures the Azure Batch executor and storage credentials) — nothing
  to bump there; the checklist item was based on a premise that doesn't
  hold for this file's current content.

**Verification (§7).**
1. Gene counts match §1 (13/80/41); spot-checked 5 genes/target for
   symbol, length, no internal stop codons (scripted check across the
   whole panel, not just the spot-check — zero internal stops found).
2. `plant_mt` set contains no plastid-only gene symbols — confirmed
   (`psbA`/`rbcL`/`ndh*` grep against `plant_mt/` returns nothing).
3. §3 measurement run and passed (with caveat above); table recorded.
4. `bash scripts/fetch_refs.sh v2026.09` on a clean checkout (fresh
   temp dir, no pre-existing `refs/`) — fetched, checksummed, unpacked,
   134 `.faa` files present.
5. `nextflow run . -profile integration` (`conf/integration.config`
   pointed at `v2026.09`) — 40/40 tasks succeeded, 10m39s,
   `tests/integration/assertions.sh` all green. Confirms the bundle
   swap changed nothing pipeline-visible (`MINIPROT_EXTRACT` and
   `ANNOTATE` are still stubs, so nothing consumes `proteins/` yet).
6. `flake8` clean on `scripts/refdata/build_proteins.py`,
   `scripts/refdata/parse_gbff_cds.py`,
   `scripts/refdata/build_manifest.py`.
