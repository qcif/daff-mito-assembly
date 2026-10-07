# Task 55 — `EXTRACT_BARCODES`: measure whether barcode regions extend past the called CDS, and if so add configurable per-locus 5′/3′ flanks

**Phase:** P3 (spec/06-phases.md) — barcode extraction, C5.
**Goal:** Establish, with measurements, whether any panel locus's
recognised barcode region runs past the CDS boundaries that stage 12
calls. If it does, let the workflow config declare a fixed number of
flanking nucleotides to include at the 5′ and/or 3′ end of each
barcode. The flanks are recorded explicitly at every output surface,
and the default is zero.

**Depends on:** nothing outstanding. Edits `bin/validate_barcodes.py`,
`bin/collate.py`, `bin/report/report.py` and
`tests/integration/assertions.sh`. Tasks 53 and 54 also touch the
report, so re-read those files before starting if either has landed.

**Structure:** Part A is a measurement with a go/no-go gate (§2.4).
Part B (§3 onward) is implemented **only if Part A says go**. If it
says no-go, the task closes after §2.5.

---

## 1. Overview

Stage 13 (`EXTRACT_BARCODES`, `bin/validate_barcodes.py`) ships each
barcode locus as **exactly** the CDS span that miniprot called in stage
12, nothing more. This is deliberate. Spec §2 stage 13 says
"coordinates are never re-derived": every emitted barcode equals a
`cds.gff` row (`pass`) or is a verbatim sub-span of one (`partial`,
task 51). `tests/integration/assertions.sh` enforces this. So the
answer to "do we currently collect the UTRs?" is **no**. No
untranslated or flanking sequence is ever emitted, and no stage of the
pipeline annotates UTRs.

**Detecting UTRs properly is not feasible here, and mostly not
meaningful:**

- **Animal mitochondria have essentially no UTRs.** The mitogenome is
  transcribed as long polycistronic precursors, which are cut apart at
  the tRNAs that sit between the genes. Mature mRNAs typically start at
  or within a few nucleotides of the start codon. Many genes end in an
  incomplete stop codon (`T` or `TA`) that is completed by
  polyadenylation after transcription.
- **Plant organellar genes do have UTRs,** but only transcript evidence
  (RNA-seq, 5′/3′ RACE) can define their boundaries. This is a DNA-only
  pipeline. No annotator in the pipeline (miniprot, MITOS2) calls UTRs,
  and spec §8 item 3 already ships `plant_pt` and `plant_mt` as
  CDS-only.

So the fallback in the request, a fixed per-barcode 5′/3′ inclusion
declared in workflow config, is the realistic mechanism. **It is only
worth building if the problem is real.** The evidence so far is
mixed:

- Every COX1 called on real data so far is full-length: 1,533 nt
  (511 codons) on `CLIENT-BC05`, `CLIENT-BC06`, `R2-INSECT-barcode05`
  and `R2-INSECT-barcode06`, and 1,536 nt on `R2-INSECT-barcode13`
  (`tests/manual/output/*/…/validation.tsv`). On `INT-ANIMAL-01`,
  miniprot's COX1 row covers residues `1 510` of its panel protein.
  The standard COI barcode (the Folmer LCO1490/HCO2198 amplicon, about
  658 bp) lies near the 5′ end of the gene and **inside** these spans.
- **Hypothesis (to be checked in Part A):** a reference COI barcode
  appears to run "outside the CDS" because of one of three edge
  effects, not because a UTR is missing:
  1. The forward primer binds at the very 5′ end of COX1. Depending on
     which start codon is annotated, it can overlap the upstream
     tRNA-Tyr by a few nucleotides. Insect COX1 start codons are
     notoriously non-canonical.
  2. miniprot sometimes fails to align a divergent panel protein's N-
     or C-terminus, so the called CDS stops short of the real gene end.
     Task 30 measured query coverage of 85–99.5 % on `INT-ANIMAL-01`,
     which leaves up to about 15 % of a protein unaligned at an end.
  3. An incomplete `T`/`TA` stop codon is not part of the aligned span.
- **One canonical plant barcode really is outside a CDS:** the
  *psbA–trnH* intergenic spacer. The panel's `psbA` entry is the gene,
  not the spacer. A spacer barcode is a non-coding region and belongs
  with rRNA/`trnL` in the "separate nucleotide-based extraction path"
  that `assets/loci.json`'s `$description` already defers. It is **out
  of scope** here (§7). A fixed flank on `psbA` is not a substitute.

Cause 2 matters for the design. Terminal truncation varies from sample
to sample, depending on how divergent the taxon is from the panel. A
fixed flank would over-pad one sample and under-pad the next while
labelling both as identical "flanks". Part A must separate truncation
(a CDS-calling problem with a different fix, §2.4) from genuine
non-coding overhang (which a fixed flank is the right tool for).

Constitution constraints that shape Part B:

- **Rule 18 / principle 7:** a flanked barcode must never be mistaken
  for a bare CDS. Flank lengths are recorded at every surface,
  including the FASTA header, `validation.tsv`, `coords.gff`,
  `metadata.json` and the report.
- **Rule 16:** the flank configuration in force is provenance and goes
  into `metadata.json`.
- **Principle 9:** flank lengths live in versioned config, never in
  code. The request says workflow config, so they go in
  `nextflow.config` params (§3.1).
- **Principle 6 and spec §8 item 3 (barcode coherence):** the locus
  shipped to Taxodactyl must stay auditable back to the annotated
  feature. Task 51 relaxed the invariant from *equal to* to *within*.
  This task relaxes it again, to *the source span (or a verbatim
  sub-span of it) plus declared flanks* (§3.3).

---

## 2. Part A — measurement (go/no-go)

### 2.0 Pin down the original observation

Before measuring, record in this task's Outcomes where "COX1 lies
outside the CDS" was seen: which sample, which reference record (BOLD
process ID or GenBank accession), and how it was compared. If it came
from a Taxodactyl hit, record the hit's alignment coordinates on the
query and the subject. That is the case Part A must explain.

### 2.1 Terminal truncation from the miniprot alignment

For every panel locus's winning `cds.gff` row (the row C5 shipped), on
every available sample, read the `Target=<protein> <qstart> <qend>`
attribute and the panel protein's length. Compute:

```
missing_5p_codons = qstart - 1
missing_3p_codons = protein_length - qend
```

Samples: the three integration fixtures (`tests/integration/output/`)
and every sample under `tests/manual/output/` (`CLIENT-BC05/06`, the
`round_2_insect` set). Tabulate per locus, per sample. This uses no
new tools. It is a read-only pass over existing outputs and the panel
FASTA in the reference bundle.

### 2.2 Reference-barcode overhang

For each `animal_mt` sample, and for `plant_pt` where a reference is
available, take a reference **barcode-region** sequence for the
nearest available taxon. Examples: a BOLD COI-5P record for COX1, or a
GenBank rbcLa/matK barcode record. Ideally use the record behind §2.0's
observation, or one from Taxodactyl's own reference set. Align it to
the sample's `target.fasta` (minimap2 or blastn). Then measure, at each
end, how many nucleotides of the reference alignment extend past the
**emitted** barcode span, oriented 5′/3′ relative to the gene.

Note whether each reference record includes the primer sequences. A
primer-inclusive record will overhang by up to about one primer length
(~25 nt) for reason 1 in §1, and that matters for the decision in
§2.4.

### 2.3 Classify each overhang

For each (locus, sample, end) with a non-zero overhang, assign one
cause:

| Cause | Signature |
|---|---|
| Terminal truncation | §2.1 shows missing codons at that end that account for the overhang |
| Primer / start-codon edge | Overhang ≤ ~30 nt at the 5′ end, no §2.1 truncation, reference is primer-inclusive or the gene has an ambiguous start |
| Incomplete stop | ≤ 2 nt at the 3′ end |
| Genuine non-coding barcode region | Anything else that is consistent across samples |

### 2.4 Decision gate

- **Go (Part B):** at least one panel locus shows an overhang that
  matters to Taxodactyl, is **consistent across samples**, and is not
  explained by terminal truncation. "Matters" means it costs query
  coverage or loses primer-adjacent sequence the reference set
  includes. The per-locus flank values are then chosen from the §2.2
  measurements: the largest consistent overhang, rounded up modestly.
  Their provenance (samples, references, date) is written beside them
  in `nextflow.config` (§3.1).
- **No-go:** every overhang is either negligible or explained by
  terminal truncation. Close the task after §2.5. Do **not** build
  flanks to mask truncation.
- **Truncation is the dominant cause** (either outcome): raise it as a
  separate `tasks/todo.md` item under *Barcode extraction*. The
  candidate remedy is to extend a called CDS outward to the panel
  protein's terminus when the alignment stops short. That re-derives
  coordinates, which spec §2 stage 13 forbids today. It is a design
  change with its own coherence implications (MITOS2's calls on
  `animal_mt` are a natural cross-check), and it is not in scope here.

**Confirm the go/no-go call with the user before starting Part B.**
Also confirm with the Taxodactyl side that a barcode carrying a few
dozen non-coding nucleotides at either end is acceptable input. This
matters for its query-coverage statistics and any CDS-anchored
alignment it performs.

### 2.5 Record the findings

Whatever the outcome:

- Write the §2.1–2.3 tables into this task's Outcomes section.
- Add a short paragraph to spec/02-stages.md stage 13 stating that
  barcodes are CDS-bounded (or CDS plus declared flanks, if Part B
  lands), and why UTR detection is not attempted (the §1 reasoning,
  condensed).

---

## 3. Part B — design (only on go)

### 3.1 Configuration

A new param in `nextflow.config`, next to `barcode_partial_min_nt`
under "Barcode extraction". It is keyed by assembly target, then by
gene symbol exactly as written in `assets/loci.json`, and gives a
5′/3′ pair in nucleotides. Example shape (illustrative values only):

```
barcode_flank_nt = [
    animal_mt: [ COX1: [5p: 30, 3p: 0] ],
    plant_pt:  [:],
    plant_mt:  [:],
]
```

- A locus absent from the map has flanks 0/0, which is exactly
  today's behaviour.
- 5′/3′ are **relative to the gene's orientation**, not the contig.
  On a minus-strand gene, the 5′ flank is taken from higher contig
  coordinates.
- A comment beside every non-zero value records its provenance from
  §2.4, in the style of `barcode_partial_min_nt`'s comment. Mark the
  values PROVISIONAL pending spec §9's benchmarking sweep. Add the
  param to spec/07-open-questions.md §9's tuning table.
- `EXTRACT_BARCODES` passes the current target's map to
  `validate_barcodes.py` as a CLI argument, for example a repeated
  `--flank GENE=5P,3P`. Do not pass it as a bare file path (rule 13).
- **Fail loudly on bad config (rule 18).** These are all hard errors
  at argument-parse time, not silent skips:
  - a gene symbol that is not in the target's `loci.json` list (a typo
    would otherwise silently ship no flank)
  - a negative value
  - a non-integer value

### 3.2 Where flanks are applied

Validation is **unchanged**. Length, identity, ORF and internal-stop
checks, the clade trial and the task 51 partial rescue all run on the
CDS exactly as now. Flanks are never translated. Flanking sequence is
non-coding, or belongs to a neighbouring tRNA, and would produce
spurious stops. Flanks are added only when building the **emitted**
sequence, after the locus's status is decided:

```
emitted_span = barcode span decided by validation (pass: whole row;
               partial: the stop-free sub-span)
for each end in (5p, 3p):
    if status == partial and this end was trimmed by the rescue:
        applied[end] = 0          # never pad back over the stop
    else:
        applied[end] = min(configured[end], room to contig edge)
extend emitted_span outward by applied[] (strand-aware)
```

Rules:

- **`partial` barcodes:** a flank is applied only to an end that
  coincides with the source row's own end. An end the rescue trimmed
  gets 0. Padding it would re-admit the CDS sequence holding the stop
  the rescue just removed.
- **Contig edges:** clip the flank at the contig boundary and record
  the applied value, not the configured value. Do not wrap across the
  origin of a circular contig. Wrapping needs circularity information
  C5 does not currently receive, and clipping is honest if it is
  recorded. Log a warning when clipping happens.
- **Spliced (multi-exon) features:** flanks extend only the outermost
  ends: the 5′ end of the first exon and the 3′ end of the last, in
  gene orientation. Internal intron handling is unchanged (see the
  `codon_blocks()` item in `tasks/todo.md`).
- **Overlap with neighbouring features** (an adjacent tRNA, or the
  ATP8/ATP6 overlap) is expected and allowed. A flank is genomic
  sequence, not a claim about feature ownership.

### 3.3 The coherence invariant — *within* becomes *within, plus declared flanks*

Today, for every emitted barcode, the CI-asserted relation is
`source_start ≤ start ≤ end ≤ source_end`. After this task, it becomes
(in contig coordinates, after strand-orienting the flanks):

```
start == cds_start - applied_flank_left
end   == cds_end   + applied_flank_right
source_start ≤ cds_start ≤ cds_end ≤ source_end
```

Here `cds_start`/`cds_end` is the validated coding span: the whole row
for `pass`, the sub-span for `partial`. The coding part of every
barcode therefore stays a verbatim span or sub-span of exactly one
`cds.gff` row. The relationship stays auditable from `validation.tsv`
and `coords.gff` alone. Update the invariant wording in:

- spec/02-stages.md, stage 13 and the C5 row
- spec/07-open-questions.md §8 item 3's rationale ("the locus drawn on
  the gene map is by construction the locus shipped"). It now reads
  "the coding part of the locus shipped".

### 3.4 Output contract

When no flank is applied, all outputs stay **byte-identical** to
today's. This is the regression guard for the default.

When a flank is applied (non-zero after clipping):

- **`barcodes.fasta`**: the record ID keeps the `GENE_seqid_start_end`
  convention, using the **emitted** (flanked) span. A header tag
  records the flanks and the coding span, alongside task 51's
  `partial=`/`source=` tags. For example:
  `>COX1_contig_10_2863_4414 flank=5p:30,3p:0 cds=2863-4384 source=2863-4384`.
- **`validation.tsv`**:
  - `start`/`end`/`length_nt` describe the emitted (flanked) sequence.
  - New columns `cds_start`, `cds_end`, `flank_5p_nt`, `flank_3p_nt`
    hold the validated coding span and the *applied* flank lengths.
    They are blank for loci that did not ship.
  - New columns are **appended**, not inserted. `assertions.sh` reads
    `source_start`/`source_end` positionally (columns 11–12) and those
    positions must not move.
- **`coords.gff`**: the emitted row for a flanked barcode is written by
  `validate_barcodes`. Its feature type must **not** be `CDS`, because
  a flanked span is not a coding sequence; use a non-coding type such
  as `region`. It carries `flank_5p=`/`flank_3p=` attributes and the
  usual `Parent=` link to the source row. The verbatim source row is
  still written as now.
- **`metadata.json`**: the per-barcode entries pick up the new
  columns (check `bin/collate.py` and
  `assets/sample_metadata.schema.json`). The run's
  `barcode_flank_nt` configuration for the sample's target is recorded
  under provenance (rule 16).
- **Report (Barcodes tab)**: a flanked barcode shows its coding length
  and flank lengths. For example "1,522 nt CDS + 30 nt 5′ flank",
  presented like task 51's "emitted / source" length display. A reader
  must never take a flanked length for a CDS length.

---

## 4. Work items (Part B)

1. `nextflow.config`: add `barcode_flank_nt` (§3.1), with values and
   provenance from §2.4.
2. `modules/local/extract_barcodes.nf`: pass the target's flank map
   to the script. Update the stage header comment.
3. `bin/validate_barcodes.py`: parse and validate the flank argument,
   apply flanks per §3.2, and write outputs per §3.4. Reuse the
   existing span and strand helpers (`extract_cds_seq`,
   `map_segment_to_genome`, `build_partial_coords_line`) rather than
   adding a second coordinate path (rule 19).
4. `bin/collate.py`, `assets/sample_metadata.schema.json`: carry the
   new fields and the provenance entry.
5. `bin/report/report.py`,
   `scripts/report/templates/components/barcodes.html`: flank display.
6. `tests/integration/assertions.sh`: replace the *within* assertion
   with §3.3's relation.
7. Spec: spec/02-stages.md stage 13 and C5; spec/07-open-questions.md
   §8 item 3 and §9; spec/06a-reports.md (Barcodes tab);
   `bin/README.md` (C5 entry).
8. Run flake8 on every touched `bin/*.py` and `scripts/tests/*.py`.
   Render a fresh `INT-ANIMAL-01` report and give its path to the user
   (CLAUDE.md).

---

## 5. Tests (Part B)

### 5.1 Unit — `scripts/tests/test_validate_barcodes.py` (extend)

100 % branch coverage on the new code (rule 14). Use small synthetic
contigs, so that no reference data is involved (rule 19). Cases:

- Zero/absent flank config → output byte-identical to the current
  fixtures. Regression guard.
- Plus-strand gene, 5′ and 3′ flanks → correct emitted coordinates,
  sequence, header tags, TSV columns and GFF row type.
- Minus-strand gene → 5′ flank taken from the higher coordinates and
  the sequence reverse-complemented correctly.
- Flank exceeding the contig start and the contig end → clipped,
  applied value recorded, warning logged.
- `partial` with the 3′ end trimmed by the rescue → 5′ flank applied,
  3′ flank 0. The mirror case also.
- Spliced two-exon feature → only the outer ends extended.
- Validation is unaffected: a locus whose flank region contains stop
  codons still `pass`es on its CDS.
- Bad config → parse-time error for an unknown gene symbol, a negative
  value, and a non-integer value.

### 5.2 Unit — `test_collate.py`, `test_report.py` (extend)

- `metadata.json` carries the new fields and the flank provenance.
- The Barcodes tab renders CDS and flank lengths for a flanked locus.
- The Barcodes tab is unchanged for an unflanked locus.

### 5.3 Integration — `tests/integration/assertions.sh`

- Replace the source-span *within* check with §3.3's relation. Read
  the new columns by header name, not position.
- If §2.4 sets a non-zero flank for an `animal_mt` locus:
  - assert that `INT-ANIMAL-01`'s barcode for that locus carries the
    `flank=` header tag and that `flank_5p_nt`/`flank_3p_nt` match the
    configured values (or less, if clipped);
  - keep the assertion keyed on the configured value, not a literal
    number, so a later re-tune does not break CI (rule 19).

### 5.4 Fixtures

No new fixture data. `INT-ANIMAL-01` already exercises a `pass` and a
`partial` (COX1) barcode on a minus-strand gene, which together cover
the strand and trimmed-end rules. `conf/stub.config` needs no change,
because the stub writes empty outputs.

---

## 6. Acceptance criteria

Part A (always):

- [ ] §2.0's original observation recorded and explained.
- [ ] §2.1 truncation table and §2.2 overhang table in Outcomes,
      covering every available sample.
- [ ] Each overhang classified (§2.3). Go/no-go decided and confirmed
      with the user. Taxodactyl's tolerance of flanked input confirmed
      if go.
- [ ] Spec stage 13 states that barcodes are CDS-bounded (plus
      declared flanks, if applicable) and why UTRs are not detected.
- [ ] If truncation is material, a `tasks/todo.md` item records it with
      its measurements.

Part B (on go):

- [ ] Default (no flanks configured) output is byte-identical to the
      pre-task output on all three integration fixtures.
- [ ] Configured flanks appear in `barcodes.fasta`, `validation.tsv`,
      `coords.gff`, `metadata.json` and `report.html`, with applied
      (post-clip) values.
- [ ] Flanks are never applied over a rescue-trimmed end, never wrap a
      contig, and never alter validation outcomes.
- [ ] Bad flank config fails at parse time.
- [ ] §3.3's coherence relation is asserted in CI.
- [ ] 100 % branch coverage on new C5 code; flake8 clean; stub run and
      integration run green.

---

## 7. Out of scope

- **Detecting UTRs** from sequence alone (§1).
- **Non-coding barcode regions** such as *psbA–trnH*, `trnL`, 12S and
  16S. These need the nucleotide-based extraction path that
  `assets/loci.json`'s `$description` already defers.
- **Extending truncated CDS calls** to the panel protein's terminus
  (§2.4). This is a separate design change, raised as a todo item if
  the measurements warrant it.
- **Wrapping flanks across the origin** of circular contigs.
- **Changing the `loci.json` schema.** Flanks live in
  `nextflow.config`, per the request. If they later need to move with
  the panel instead (principle 9), that is a schema bump and a
  separate task.

---

## Outcomes

_To be completed on execution._
