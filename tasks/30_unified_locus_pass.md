# Task 30 — Stages 12–13: unified locus pass (`MINIPROT_CDS` + `EXTRACT_BARCODES`)

**Phase:** P3/P4 (from [spec §6](../spec/06-phases.md)).
**Prerequisite:** [task 29](completed/29_comprehensive_protein_panel.md) — the
comprehensive protein panel is in the bundle and its §3 go/no-go has
passed.

## 0. Overview

This task implements barcode extraction — the pipeline's **contractual
output** ([brief.md §5](../brief.md)) — for the first time. Stage 13 has
been a `touch`-only stub since P0.

It implements it in the unified shape: **one miniprot pass** against the
comprehensive organellar protein panel produces every protein-coding
feature on the binned contigs, and the barcodes are the **subset** of
those features whose gene symbol appears in
[`assets/loci.json`](../assets/loci.json). The same alignment that draws
`rbcL` on the report's gene map is the one that yields the `rbcL`
sequence in `barcodes.fasta`. They cannot disagree, because they are the
same object.

The alternative — which this replaces before it is ever built — would
run miniprot twice over two different protein sets and hope the two
agreed. They usually would. "Usually" is not a property you want in a
biosecurity audit trail ([rule 18](../CONSTITUTION.md)).

**Two processes, because two containers.** miniprot is an off-the-shelf
biocontainer; the barcode validation is C5, custom Python in its own
image ([spec §2.2](../spec/02-stages.md#22-custom-logic-components)).
Splitting them also makes the subset relationship explicit in the DAG
rather than buried inside one script.

**Goal:** implement stages 12 and 13 as `MINIPROT_CDS` (all CDS
features) and `EXTRACT_BARCODES` (validated barcode subset), replacing
the `MINIPROT_EXTRACT` stub.

**Exit criteria:**

- `-profile integration` produces, for every sample reaching
  `BIN_TARGET`:
    - `results/<sample>/cds/<sample>.cds.gff` — all protein-coding
      features found on `target.fasta`;
    - `results/<sample>/barcodes/barcodes.fasta` — validated barcode
      subset, `results/<sample>/barcodes/<sample>.validation.tsv` —
      per-locus pass/fail with reasons.
- Every record in `barcodes.fasta` traces to a feature in
  `cds.gff` by identical coordinates — asserted, not assumed (§6).
- `INT-ANIMAL-01` recovers `COX1` and `CYTB`; `INT-PLANT-01-pt` recovers
  `rbcL` and `matK`, matching `expected/*/expected_loci.txt`.
- C5 exists with a `scripts/tests/` module at **100 % branch coverage**
  ([rule 14](../CONSTITUTION.md)).
- `-profile stub -stub-run` green; `flake8` clean.

**Not in scope:**

- **tRNA and rRNA features.** miniprot is protein-to-genome only.
  `loci.json`'s own `$description` already records that rRNA barcodes
  (12S, 16S) and tRNA barcodes (trnL) need a separate nucleotide
  extraction path — unchanged, and out of scope here as it always was.
  Non-CDS *annotation* features come from
  [task 31](31_annotate_merge.md).
- **The annotation output contract.** `cds.gff` is an input to
  `ANNOTATE`, not the final annotation. Task 31 owns
  `annotation_summary.json` and the merge.
- **Changing which loci are barcodes.** `loci.json` is unchanged and
  contractually coupled to Taxodactyl
  ([principle 9](../CONSTITUTION.md)).

---

## 1. Stage numbering

The unified design inverts a dependency: annotation now consumes the
locus pass rather than running beside it. The flow diagram in
[spec §1](../spec/01-pipeline-flow.md) must match execution order.

| Slot | Was | Becomes |
|---|---|---|
| 12 | `ANNOTATE` | **`MINIPROT_CDS`** |
| 13 | `MINIPROT_EXTRACT` | **`EXTRACT_BARCODES`** |
| 13a | — | **`ANNOTATE`** (new sub-slot, [task 31](31_annotate_merge.md)) |
| 14–16 | `ORGANELLE_MAP`, `COLLATE`, `RUN_REPORT` | unchanged |

Stage 13 keeps its meaning — *the barcode stage* — which is how the spec
refers to it in a dozen places. `ANNOTATE` takes a sub-slot rather than
forcing a renumber of 14–16; the spec already tolerates non-contiguous
numbering, keeping stage 8 as a reserved slot for exactly this
stability reason ([spec §2 stage 8](../spec/02-stages.md#2-stage-detail)).

Update [spec §1](../spec/01-pipeline-flow.md) flow diagram,
[§2](../spec/02-stages.md#2-stage-detail) stage table, and
[§3.2](../spec/03-organelles.md#32-per-stage-parameter-selection). In
§2 stage 13's notes, the sentence *"Runs independently of ANNOTATE — the
barcode panel is the contractual output; ANNOTATE is supplementary"*
must be rewritten: the barcode panel is still the contractual output and
ANNOTATE is still supplementary, but the dependency now runs
ANNOTATE → EXTRACT_BARCODES' shared source, and the coherence that buys
is the point.

## 2. `MINIPROT_CDS` (stage 12)

Replace `modules/local/annotate.nf`'s slot with
`modules/local/miniprot_cds.nf`.

```groovy
input:
tuple val(meta), path(target_fasta), path(blast_tsv)
path protein_panel          // refs/<ver>/proteins/ — value channel

output:
tuple val(meta), path(target_fasta),
      path("${meta.sample_id}.cds.gff"), emit: cds

script:
"""
if [[ -s ${target_fasta} ]]; then
    cat ${protein_panel}/${meta.assembly_target}/*.faa > panel.faa
    miniprot --gff -t ${task.cpus} ${target_fasta} panel.faa \\
        > ${meta.sample_id}.cds.gff
else
    : > ${meta.sample_id}.cds.gff
fi
"""
```

- **Container:** `quay.io/biocontainers/miniprot:0.18--h577a1d6_0`,
  SHA-pinned in [`conf/containers.config`](../conf/containers.config)
  (replacing the `TODO P3` stub). Resolve the digest as in
  [task 21 §1](completed/21_blast_validate.md#1-container-pin--confcontainersconfig).
- **Empty-query guard.** Two upstream states legitimately produce an
  empty `target.fasta` — C3 selected nothing, or a `plant_pt` sample
  withheld the C4 substitution
  ([task 24 §3.2](completed/24_plastid_substitution_guard.md)). Same
  short-circuit as [task 21 §2](completed/21_blast_validate.md), and for
  the same reason: an empty output file, not a Nextflow error.
- **Panel staged as a `path` input**, not `${params.x}` —
  [spec §1a](../spec/01-pipeline-flow.md#1a-engineering-constraints)
  pattern 1, materialised as a value channel at the top of `main.nf`
  alongside `ch_organelle_refs`.
- **`target_fasta` is re-emitted** so `EXTRACT_BARCODES` can pull
  sequence without a second join.
- **No `--outn` / hit-count limit by default.** The panel is curated and
  target-specific; capping hits per query would silently drop a genuine
  second copy — and a plastid gene inside the inverted repeat *is* a
  genuine second copy ([spec §3.6](../spec/03-organelles.md#36-plastid-quadripartite-canonicalisation)).
  Do not de-duplicate by gene name anywhere in this pass.

## 3. `EXTRACT_BARCODES` (stage 13) — C5 `bin/validate_barcodes.py`

C5 is already specified in
[spec §2.2](../spec/02-stages.md#22-custom-logic-components); this task
implements it, with one change to its inputs — it reads `cds.gff` rather
than running miniprot itself.

```
validate_barcodes.py \
    --cds-gff <sample>.cds.gff \
    --target-fasta target.fasta \
    --assembly-target animal_mt \
    --locus-panel assets/loci.json \
    --genetic-codes 2,5 \
    --min-identity 60 \
    --out-fasta barcodes.fasta \
    --out-coords <sample>.coords.gff \
    --out-tsv <sample>.validation.tsv
```

**Behaviour:**

1. **Subset.** Select features whose gene symbol appears in
   `loci.json[<origin>]`. Matching is **case-insensitive** —
   `loci.json` uses `COX1` for `animal_mt` and `cox1` for `plant_mt`,
   and NCBI symbols are case-insensitive
   ([task 29 §4](completed/29_comprehensive_protein_panel.md)).
2. **Extract** the nucleotide sequence at each selected feature's
   coordinates from `target.fasta`, respecting strand.
3. **Validate** per [spec §2.2 C5](../spec/02-stages.md#22-custom-logic-components):
   length, ORF under the target-appropriate NCBI table, internal stops,
   protein-identity floor.
4. **`animal_mt` clade trial** ([spec §3.3](../spec/03-organelles.md#33-specific-issues-and-decisions)):
   try tables 2 (vertebrate) and 5 (invertebrate), pick the one yielding
   a valid ORF, record the chosen table. `params.genetic_code_tables`
   already carries `animal_mt: [2, 5]`, `plant_pt: [11]`,
   `plant_mt: [1]` — a single-element list means no trial.
5. **Emit** `barcodes.fasta` (validated only), `coords.gff` (the
   selected features, coordinates **byte-identical** to their `cds.gff`
   source rows), and `validation.tsv` recording every panel locus
   including those that dropped out and why.
6. **Always exit 0** ([principle 8](../CONSTITUTION.md)) — a sample with
   no recoverable barcode is a `no_barcode` result
   ([principle 7](../CONSTITUTION.md)), not a pipeline error.

**Container:** `wf5/barcode-validate:<tag>` per spec §2.2 (Python +
`biopython`). Add to [`.github/workflows/images.yml`](../.github/workflows/images.yml)
if not already present; source stages at runtime from `bin/`.

**The coherence invariant, stated once so it can be tested:** every
record in `barcodes.fasta` corresponds to exactly one feature in
`cds.gff`, at identical coordinates. C5 must never re-align, re-trim, or
adjust coordinates — if a barcode needs different bounds than the
annotation shows, that is a bug in one of them, not something to paper
over. Unit test §5 case 2 and integration assertion §6 both pin this.

## 4. Wiring

```groovy
MINIPROT_CDS(BLAST_VALIDATE.out.validated, ch_protein_panel)
EXTRACT_BARCODES(MINIPROT_CDS.out.cds)
// ANNOTATE(MINIPROT_CDS.out.cds) — task 31
```

`COLLATE`'s join takes `EXTRACT_BARCODES.out.barcodes` in place of
`MINIPROT_EXTRACT.out.barcodes`; the tuple shape
`(meta, barcodes.fasta, coords.gff, validation.tsv)` is unchanged, so
`collate.nf` needs no input change.

## 5. Unit tests — `scripts/tests/test_validate_barcodes.py`

100 % branch coverage on `bin/validate_barcodes.py`. External tools
mocked at the subprocess boundary per
[task 27](completed/27_unit_test_boundary_mocking.md) — though C5 shells
out to nothing, so fixtures are small on-disk GFF + FASTA pairs.

1. **Happy path** — a `cds.gff` with 20 features, 6 of them panel loci;
   exactly those 6 are extracted and validated.
2. **Coherence invariant** — every emitted record's coordinates match
   its `cds.gff` source row exactly. The regression test for §3's
   invariant.
3. **Case-insensitive symbol match** — `COX1` in the panel matches
   `cox1` in the GFF and vice versa.
4. **Strand handling** — a minus-strand feature yields the reverse
   complement.
5. **`animal_mt` clade trial** — a sequence valid under table 5 but not
   table 2 selects 5 and records it; one valid under both selects
   deterministically (first in list) and records that.
6. **Neither table valid** → locus fails with an explicit reason in
   `validation.tsv`, absent from `barcodes.fasta`, exit 0.
7. **Internal stop codon** → fails with that reason.
8. **Below identity floor** → fails with that reason.
9. **Panel locus absent from `cds.gff`** → present in `validation.tsv`
   as `not_found`, absent from FASTA. Negative clarity: a locus that was
   looked for and not found is different from one never in the panel.
10. **All loci fail** → empty `barcodes.fasta`, populated
    `validation.tsv`, exit 0.
11. **Empty `cds.gff`** (empty `target.fasta` upstream) → all loci
    `not_found`, exit 0.
12. **Duplicate gene at two loci** (the plastid IR case) — both are
    emitted, with distinct sequence IDs carrying coordinates. Not
    de-duplicated.
13. **Malformed GFF line** → skipped with a warning, remaining features
    processed, exit 0.
14. **Exit code 0 in every case** — parametrised.

## 6. Integration-test wiring

Existing `tests/integration/expected/<target>/expected_loci.txt` files
already carry the expected locus lists — use them rather than inventing
new fixtures.

Assertion block in
[`tests/integration/assertions.sh`](../tests/integration/assertions.sh),
over `ASSEMBLING_SAMPLES`:

- `barcodes.fasta` non-empty and contains every locus in that target's
  `expected_loci.txt`;
- `validation.tsv` has one row per panel locus (found or not);
- **the coherence assertion** — for each record in `barcodes.fasta`,
  a row exists in `cds.gff` at the same seqid/start/end. This is the
  whole point of the architecture; assert it in CI, not just in unit
  tests.

`INT-PLANT-01-mt` is not asserted — it soft-fails the coverage gate and
never reaches `BIN_TARGET`
([task 25](completed/25_coverage_gate_carryover.md),
[task 28](completed/28_plastid_masked_mt_panel.md)). Standing consequence
already documented in
[task 21 §4](completed/21_blast_validate.md#4-integration-test-wiring).

Tick the extraction row in the progressive-uncomment header.

## 7. Verification

1. `bash scripts/pytest.sh` — C5 green at 100 % branch coverage.
2. `flake8 bin/validate_barcodes.py scripts/tests/test_validate_barcodes.py`
   (79-col) via the `claude` venv.
3. `nextflow run . -profile stub -stub-run` — green; stub process count
   updated for the 12/13 split.
4. `nextflow run . -profile integration` — expected loci recovered for
   both assembling fixtures; record per-locus identity and the chosen
   `animal_mt` genetic-code table in Outcomes.
5. `bash tests/integration/assertions.sh` — new blocks green,
   pre-existing blocks unaffected.
6. CI green on `Tests` and nightly `Integration`.

## 8. Deliverables checklist

- [ ] `conf/containers.config` — SHA-pinned miniprot; `wf5/barcode-validate`
      pinned; `TODO P3` comments removed.
- [ ] `modules/local/miniprot_cds.nf`, `modules/local/extract_barcodes.nf`;
      `miniprot_extract.nf` removed.
- [ ] `bin/validate_barcodes.py` (C5) — always exits 0.
- [ ] `scripts/tests/test_validate_barcodes.py` — 100 % branch coverage.
- [ ] `main.nf` — `ch_protein_panel` value channel, new wiring, `COLLATE`
      join updated.
- [ ] `tests/integration/assertions.sh` — extraction + coherence blocks;
      header tick.
- [ ] Spec: [§1](../spec/01-pipeline-flow.md) flow diagram,
      [§2](../spec/02-stages.md#2-stage-detail) stages 12/13,
      [§2.2](../spec/02-stages.md#22-custom-logic-components) C5 inputs,
      [§3.2](../spec/03-organelles.md#32-per-stage-parameter-selection).

## 9. Notes / non-issues

- **Why the contractual output may depend on the annotation pass.** It
  does not, in the sense that matters: there is no *extra* failure mode.
  The barcode still comes from a miniprot alignment, as it always would
  have; the only change is that the same alignment is also published as
  an annotation feature. What has been removed is the possibility of two
  alignments disagreeing.
- **Why not one process.** Two containers — an off-the-shelf
  biocontainer and a custom-logic image — and
  [rule 14](../CONSTITUTION.md) wants custom logic isolated and
  unit-tested. The split also makes the subset relationship visible in
  the DAG.
- **Do not de-duplicate by gene name.** Stated in §2, §3 and tested in
  §5 case 12. A plastid gene in the inverted repeat genuinely occurs
  twice, and collapsing it destroys the evidence the IR is present.
- **`params.locus_panel`** now points at `assets/loci.json` as a pure
  selector; the protein sequences come from `ch_protein_panel`. Two
  params, two roles — do not merge them.

## 10. Outcomes

_(fill on completion: miniprot + barcode-validate digests, per-locus
identity per fixture, chosen `animal_mt` genetic-code table, stage wall
times.)_
