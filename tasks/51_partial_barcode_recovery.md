# Task 51 — `EXTRACT_BARCODES`: recover a trimmed, `partial` barcode when an otherwise-good locus carries an internal stop

**Phase:** P3 (spec 06-phases.md) — barcode extraction, C5.
**Goal:** When a barcode locus is strongly supported but fails only on
an internal stop codon, ship the longest stop-free stretch of its CDS
as an explicitly-labelled `partial` barcode instead of dropping the
locus. Ship it only if that stretch is at least a configurable minimum
length (default 100 nt).

**Depends on:** nothing outstanding. Touches files that task
47_bandage_interactive_svg.md also edits (`bin/report/report.py`,
`bin/collate.py`, `tests/integration/assertions.sh`), so re-read those
files before starting if task 47 has landed in the meantime.

**Folds in:** the `tasks/todo.md` item under *Barcode extraction* that
ends "`COX1` at 72.8% identity still fails on `internal_stop_codon`
despite clearing the floor — worth a follow-up look…". Once this task
lands, delete that sentence from todo.md. The rest of that item (the
identity-floor and panel-breadth discussion) stays, because it is a
separate issue.

---

## 1. Overview

For each locus in the barcode panel (`assets/loci.json`), stage 13
(`EXTRACT_BARCODES`, `bin/validate_barcodes.py`) takes the best miniprot
protein-to-genome alignment on the assembled organelle contig. It then
checks that the aligned CDS is a plausible open reading frame before
shipping it to Taxodactyl. One of those checks rejects the whole locus
if its translation contains an **internal stop codon** under every
configured genetic code. The check exists for a good reason (brief.md
§7): a nuclear mitochondrial pseudogene (NUMT) typically shows up as a
broken ORF, and a stop codon is the clearest sign of one.

On the `INT-ANIMAL-01` integration fixture (pea aphid mitogenome), 5 of
6 `animal_mt` barcodes pass and **COX1 is dropped**. COX1 is the
standard animal barcode and the one Taxodactyl users care about most.
But this COX1 is not a pseudogene. It is a near-full-length,
high-identity gene with **one** stop codon about 35 residues from its C
terminus, which is almost certainly an ONT assembly error (§2). Well
over 90% of the gene, including the whole Folmer barcode region, is an
intact reading frame. The same feature is already in the annotation
GFF and on the report's gene map. We draw the gene but withhold its
sequence.

This task keeps the stop-codon check but changes what happens when it
fires. Rejecting the whole locus is the wrong response when the stop
sits in a short, isolated part of an otherwise strong alignment. We
trim to the longest stretch that *is* a clean ORF and ship it, clearly
labelled `partial`, so the user gets a usable but visibly truncated
barcode instead of nothing. If no clean stretch is long enough to be
useful (under `params.barcode_partial_min_nt`, default 100 nt), the
locus still fails exactly as it does today.

This follows CONSTITUTION principle 7 ("a partial recovery is not a
negative… must be distinguishable from a full-depth `ok` at every
surface, *and* must not be suppressed or reported as a failure merely
for being incomplete"). It also follows principle 18: the trim is
recorded, auditable and labelled, never silent.

---

## 2. What was verified before this brief was written

These results come from `tests/integration/output/INT-ANIMAL-01/` (run of
2026-09-24), with `bin/validate_barcodes.py`'s own functions run
inside `neoformit/daff-wf5-scripts:test`. They are settled, so do not
re-litigate them.

`INT-ANIMAL-01.validation.tsv`:

```
COX1  fail  internal_stop_codon  contig_10  2863  4384  -  0.9353  (blank)  1522
```

The winning hit is `MP001229`, from reference `NC_035301.1_COX1`, covering
protein 1–510 of 510. It has `Identity=0.9353`, `Positive=0.9667`,
`Frameshift=4`, `StopCodon=1`, and CIGAR `62M1G152M2G44M2G183M2G29M1I35M`.
That is the highest-identity barcode hit in the sample; the five loci
that pass range from 0.77 to 0.88.

Running `codon_blocks()` and translating under table 5 (the selected
table in `genetic_code.json`) gives:

| Block | nt | aa | Stops |
|------:|----:|----:|-------|
| 0 | 186 | 62 | — |
| 1 | 456 | 152 | — |
| 2 | 132 | 44 | — |
| 3 | 549 | 183 | — |
| 4 | 192 | 64 | **aa 29** (`TAA`, aligned to ref `K`; miniprot `cs` `*taaK`) |

- The stop sits right next to the only `I` op in the CIGAR, in the last
  block, just after the fourth frameshift break. This looks like a
  local indel/homopolymer error in the ONT assembly that miniprot
  resolved imperfectly, not a genuine truncation.
- The first four blocks and the first 29 codons of block 4 are clean.
  Counting the 7 frameshift-skipped bases, the stop-free 5′ stretch is
  sequence offsets 0–1416, i.e. **1,417 nt ≈ 93% of the 1,522 nt
  span**. It runs from genome coordinates 2968 to 4384 on the minus
  strand. Recompute this in the implementation; do not hard-code it
  into tests. The Folmer (LCO1490/HCO2198) region, roughly the first
  ~700 nt of COX1, is entirely inside it.
- The annotation already carries this feature. `organelle_annotation.gff`
  has the same `MP001229` row (and a `NC_063091.1_cox1` hit with the
  same `StopCodon=1`), so the gene map draws COX1 at 2863–4384 while
  `barcodes.fasta` omits it.
- Falling back to another hit in the cluster is a worse fix, so it is
  **rejected**. The cluster holds 1,091 COX1 hits: 1,027 fail
  `internal_stop_codon`, 21 fail on identity, and 43 pass. The passing
  hits top out at 0.738 identity. The best of them (`NC_086791.1`)
  "passes" with a spurious `328N` intron. The rest pass only because
  their alignments end before the stop codon. Switching to them would
  give a lower-quality alignment with arbitrary boundaries, and it
  would break the "barcode = annotated feature" story.
- This is fixture-specific, not systemic. Client samples `CLIENT-BC05`
  and `CLIENT-BC06` (`tests/manual/output/`) pass COX1 at 1,533 nt with
  no stop. Keep the change narrow: it must not alter any currently
  passing locus.
- Barcodes that already pass also carry indel errors: COX2, COX3, CYTB,
  ND1 and ATP6 all have `Frameshift=1..6`. Those frameshifted bases ship
  as-is today, which is intended (validate_barcodes.py module
  docstring). A trimmed partial therefore contains the same kind of
  imperfection as a normal pass, and that is acceptable. Do **not** try
  to "correct" frameshifts or remove the stop codon's bases in place.
  Only trim.

---

## 3. Design

### 3.1 When a locus becomes `partial`

A candidate becomes `partial` only if **all** of these hold:

1. It already clears the length and identity-floor checks. Order in
   `validate_orf` is unchanged: `invalid_length` and
   `identity_below_floor` still fail first and are never rescued.
2. It fails the internal-stop check under **every** configured table.
   A candidate that passes under any table is a plain `pass`, exactly as
   today.
3. The longest stop-free segment (§3.2), under the table that gives
   the longest one, is ≥ `params.barcode_partial_min_nt`.

Otherwise it fails with `internal_stop_codon`, exactly as today.

Consider, but do not require without evidence, a second guard: skip
the rescue when the hit has many internal stops. Many stops point to a
real pseudogene rather than a single assembly error. If you add one,
put it in config next to the minimum length (e.g.
`barcode_partial_max_stops`). Otherwise leave a one-line note in the
Outcomes section of this task explaining why it was not needed. The
length floor alone already rejects hits whose stops are spread through
the gene.

### 3.2 Finding the longest stop-free segment

Work on the same codon blocks the pass/fail check already uses: use
`codon_blocks()` when a CIGAR is present, otherwise the
whole-sequence-truncated fallback. There should be one definition of a
codon and one translation path (CONSTITUTION rule 19).

- Frameshift/intron breaks between blocks do **not** end a segment,
  because frameshifts are tolerated in a normal pass. Only a stop codon
  ends one.
- The stop codon's own three bases are excluded from both neighbouring
  segments.
- The normal terminal-stop allowance still applies. A stop as the final
  codon of the final block is not internal, and it may be kept at the
  end of a segment that reaches it.
- To choose between tables, use the longest segment. Break ties by
  configured table order, as the existing clade trial does. Record the
  chosen table in `genetic_code`, as for a pass.
- For ties between segments of equal length within one table, prefer
  the 5′-most, so the result is deterministic.

The block logic needs a map from sequence offsets back to genome
coordinates. Every block boundary has to carry its start offset in the
spliced sequence, and the dropped frameshift bases must be counted.
`codon_blocks()` currently discards the cursor. Extend it, or add a
sibling helper, so it returns the offsets. **Do not** change what the
function returns to existing callers:
`bin/translate_annotation_cds.py` imports from this module.

Pseudocode:

```
segments = []
for block (with seq_offset) in blocks:
    for codon_index, aa in translate(block):
        if aa == '*' and not (last block and last codon):
            close current segment at seq_offset + codon_index*3
            open new segment at seq_offset + codon_index*3 + 3
close final segment at end of last block
best = max(segments, key=length)   # nt, counted on the spliced sequence
```

Take the emitted sequence as the contiguous slice
`seq[best.start:best.end]` of the spliced CDS. That keeps any
frameshift bases inside the segment. It ships the assembly verbatim,
just shorter, which is consistent with §2's last bullet.

Map offsets to genomic coordinates through the feature's CDS exon list
(`feature['cds']`) and strand, with the same ordering
`extract_cds_seq` uses. On `-` strand, offset 0 is the highest
coordinate. A multi-exon partial should not occur in the current
panels (see the todo.md note on `codon_blocks()` and introns). Don't
build a general spliced-subrange emitter for it. If the segment would
span an exon junction, fail the locus with `internal_stop_codon` and
leave a comment saying why. Do not guess.

### 3.3 The coherence invariant — relaxed from *equal to* to *within*

Today the invariant is that every `barcodes.fasta` record matches a
`cds.gff` row **exactly**. It is asserted in
`tests/integration/assertions.sh` and stated in spec 02-stages.md stage
13, spec 02-stages.md C5 and spec 05-test-data.md. A partial breaks that
by construction, so the invariant becomes:

> Every `barcodes.fasta` record is **either** identical to a `cds.gff`
> source row (`pass`) **or** a contiguous sub-span of exactly one
> `cds.gff` source row (`partial`). Both cases are auditable back to
> that row from `validation.tsv` and `coords.gff` alone.

The guarantee this invariant was protecting still holds: there is no
second alignment and no re-derived gene model. The barcode is still
cut from the annotated feature, only shorter. Update the three spec
passages above to say so. Keep the new wording as strict as possible,
because the spec explicitly flags this invariant as "the one thing that
would silently rot".

### 3.4 Output contract

**`validation.tsv`.** Add a new `status` value `partial`, alongside
`pass`, `fail` and `not_found`. Its `reason` is `internal_stop_codon`,
so the TSV still records why the locus is incomplete. Add columns
(append them; don't reorder the existing ones):

| Column | `pass` | `partial` | `fail` / `not_found` |
|---|---|---|---|
| `start`, `end` | source row | **emitted** sub-span | source row / blank (unchanged) |
| `length_nt` | emitted length | **emitted** length | source length / blank (unchanged) |
| `source_start`, `source_end` | = `start`, `end` | the `cds.gff` row's span | same as `start`/`end` / blank |
| `source_length_nt` | = `length_nt` | full spliced length | same as `length_nt` / blank |
| `n_internal_stops` | 0 | count under the chosen table | count under the first table, for `internal_stop_codon` fails; blank otherwise |

`start`/`end`/`length_nt` always describe **what ships**. That keeps
the barcode ID (`{gene}_{seqid}_{start}_{end}`), `report.py`'s
`_barcode_id()` lookup and the existing columns correct with no
special cases. The `source_*` columns are what an auditor uses to find
the row in `cds.gff`.

**`barcodes.fasta`.** Keep the ID format unchanged. For a partial, add
a FASTA description after the ID, separated by whitespace, e.g.
`>COX1_contig_10_2968_4384 partial=internal_stop_trimmed source=2863-4384`.
This follows standard FASTA semantics: the ID token is unchanged, so
tools that key on the ID are unaffected. Leave `pass` records exactly
as they are today, with no description.

**`coords.gff`.** For a partial, write the source `cds.gff` row(s)
verbatim, as for a pass, so the source is always reproduced exactly.
Then add one derived line for the emitted span with source column
`validate_barcodes` and attributes that link it back to the source
row (`Parent=<source mRNA ID>`) and mark it (`partial=true`,
`trim_reason=internal_stop_codon`). Either this or an equivalent
lossless encoding is fine. Document the choice in the module
docstring.

**Config.** Add `params.barcode_partial_min_nt = 100` in
`nextflow.config`, with a comment marking it as provisional pending the
spec 07-open-questions.md §9 benchmarking sweep. Pass it to
`validate_barcodes.py` as `--partial-min-nt` from
`modules/local/extract_barcodes.nf`. Don't hard-code it in the script
or module (CONSTITUTION rule 18: versioned config over hardcoded
value). Setting it to `0` must **disable** partial recovery entirely
and restore today's behaviour. Test that explicitly.

**Why 100 and not Taxodactyl's floor.** Taxodactyl accepts partial
sequences down to 20 nt (confirmed with the user, 2026-09-29). That is
the *consumer's* hard floor, not a target. This workflow deliberately
sets a stricter default so it doesn't encourage downstream analysis of
fragments too short to identify anything reliably. Put that reasoning
in the config comment so nobody later "aligns" the default down to 20.
Validate the param: `0` (disabled) or any value ≥ 20 is accepted. A
value from 1 to 19 fails fast at `validate_barcodes.py` argument
parsing with a message naming Taxodactyl's 20 nt floor, because
shipping a partial that the consumer will reject is worse than not
shipping one (rule 18). Keep `20` as a named constant in the script.

### 3.5 Downstream surfaces — partial must be visible everywhere

Principle 7 requires a partial to be distinguishable from a full pass
at every surface, and not reported as a failure.

- **`bin/collate.py`.** `barcodes_section` keeps `n_passed`, still
  counting `pass` only, and adds `n_partial`. Sample status
  `no_barcode` means **zero** recovered barcodes, and a sample with only
  partials is **not** `no_barcode`. Check `n_loci_passed()` and every
  caller that decides `no_barcode`, and count `pass + partial` there.
  `barcodes.fasta` is already copied whole, so nothing changes there.
- **`bin/report/report.py`** (per-sample, Barcodes tab and Overview key
  findings):
  - Partials appear alongside passes with a sequence and a download,
    carrying a visible `partial` badge and warning styling. Show
    emitted vs. source length ("1,417 / 1,522 nt, trimmed at internal
    stop").
  - Keep them out of the "dropped" list.
  - `_barcode_count_finding` and the `barcode_pass_fraction` key
    finding must not count a partial as a full pass. Pick one wording,
    e.g. "5 of 6 loci recovered in full, 1 partial", and use it
    consistently.
  - `BARCODE_DROPOUT_REASONS['internal_stop_codon']` still describes
    true failures. Add separate explanatory text for the partial case.
- **Run report (`bin/run_report.py`, `bin/report/run_report.py`).** The
  panel-recovery column must show partials distinctly (e.g. `5 + 1
  partial / 6`). Don't fold them into the passed count.
- **Annotation.** No change. The annotation GFF keeps the full miniprot
  feature, which is correct, because the gene *is* there. The report's
  explanatory text should mention that a partial barcode is a trimmed
  sub-span of the feature drawn on the map (spec 06a-reports.md tab 4
  already says the two are "the same object"; adjust that sentence).
- **`metadata.json`.** Nothing new is needed beyond the columns that
  flow through `barcodes.loci`. Confirm that the new columns show up in
  metadata and don't break `bin/run_report.py`'s reader.

---

## 4. Work items

1. `bin/validate_barcodes.py`:
   - partial-recovery logic (§3.1, §3.2)
   - offset-carrying block helper
   - new TSV columns, FASTA description and coords.gff line (§3.4)
   - `--partial-min-nt` CLI argument
   - update the module docstring. The "Never re-aligns or re-trims"
     paragraph needs rewording: it still never re-aligns, and it now
     trims only in the one documented, labelled case.
2. `modules/local/extract_barcodes.nf` and `nextflow.config`: new param
   and its wiring. The stub is unchanged.
3. `bin/collate.py`, `bin/report/report.py`, the Barcodes-tab template
   under `bin/report/`, `bin/run_report.py`, `bin/report/run_report.py`:
   §3.5.
4. Spec updates:
   - spec 02-stages.md stage 13 and C5 rows (the invariant wording in
     §3.3, plus the `partial` status)
   - spec 05-test-data.md acceptance criteria (the invariant)
   - spec 06a-reports.md tab 4 (partial listed alongside passes, not in
     the dropout list; "same object" sentence) and the run-report
     recovery column
   - `bin/README.md` if C5's description changes
5. `tasks/todo.md`: remove the COX1 sentence named under **Folds in**.

---

## 5. Tests

### 5.1 Unit — `scripts/tests/test_validate_barcodes.py` (extend)

Run with `scripts/pytest.sh`. Build the fixtures synthetically:
hand-written GFF with `##PAF` lines and short synthetic contigs. Don't
copy the 12 MB fixture `cds.gff`. Cover:

- **No regression:** a clean ORF is still `pass`, with `source_*` equal
  to `start`/`end`, no FASTA description, and byte-identical FASTA
  record.
- **Rescue at the 3′ end:** a stop near the 3′ end, `+` strand, gives
  `partial` with the 5′ segment, correct coordinates, `n_internal_stops
  == 1`, and emitted sequence equal to the exact slice.
- **`-` strand mapping:** the same case on `-` strand gives coordinates
  that are correct against the reverse-complemented slice. This is
  where off-by-ones will hide.
- **Rescue at the 5′ end:** a stop near the 5′ end selects the 3′
  segment.
- **Frameshifts:** a CIGAR with `G`/`F` breaks before and after the
  stop gives offsets that account for the skipped bases (mirror the §2
  COX1 CIGAR shape).
- **Length floor:**
  - longest segment `< min` gives `fail internal_stop_codon`, with
    `n_internal_stops` populated
  - a segment exactly `== min` gives `partial`
- **Disabled:** `--partial-min-nt 0` reproduces today's behaviour
  (fail).
- **Consumer floor:** `--partial-min-nt` values 1–19 are rejected at
  argument parsing with a message naming the 20 nt floor. 20 is
  accepted.
- **Earlier failures are never rescued:** `identity_below_floor` and
  `invalid_length` candidates still fail with their own reasons.
- **Table choice:** a candidate that passes under table 2 but has stops
  under table 5 is a plain `pass` under table 2, not a partial under
  table 5. A candidate with stops under both tables picks the table
  with the longer segment.
- **Terminal stop:** a terminal stop is never counted as internal and
  never ends a segment early.
- **Exon junction:** a segment that would span an exon junction fails
  (§3.2).
- **Regression guard:** the existing callers of `codon_blocks()` /
  `extract_cds_seq`, including `translate_annotation_cds`'s tests,
  still pass unchanged.

Keep branch coverage of `validate_barcodes.py` at or above its current
level.

### 5.2 Unit — `test_collate.py`, `test_report.py`, `test_run_report.py` (extend)

- A sample whose only recovered locus is partial is **not**
  `no_barcode`.
- `n_partial` is present and correct.
- The Barcodes tab renders a partial with its badge, both lengths and a
  sequence, and does not list it under dropped loci.
- The key finding and the run-report recovery column both distinguish
  partial from pass.

### 5.3 Integration — `tests/integration/assertions.sh`

- **Replace** the assertion that `INT-ANIMAL-01/report.html` carries
  the COX1 `internal_stop_codon` dropout reason (task 43b §5.8). That
  assertion now describes the bug. Instead assert that INT-ANIMAL-01's
  `validation.tsv` has COX1 as `partial` with `length_nt ≥
  barcode_partial_min_nt`, and that `barcodes.fasta` carries a `>COX1_`
  record with `partial=` in its description.
  - Don't pin the exact 2968–4384 coordinates. The assembly is not
    stable across recruitment changes (todo.md, *Recruitment*).
  - Assert that the emitted span lies inside the source span.
- **Update the coherence check** ("every `barcodes.fasta` record matches
  a `cds.gff` row"):
  - `pass` records still need an exact match.
  - `partial` records need their `source_start`/`source_end` to match a
    `cds.gff` row exactly **and** their emitted span to lie within it.
  - Any other mismatch remains a hard `FAIL`.
- Locus counting (`found`/`candidates`) should count a partial as
  recovered, but print it distinctly (e.g. `OK: INT-ANIMAL-01 6/6
  barcode-panel loci recovered (1 partial)`).

### 5.4 Fixtures

No new data fixture is needed: `INT-ANIMAL-01` already exercises the
case. Regenerate any expected-output files that record per-locus
barcode status or counts, and say which ones in Outcomes. The plant
fixtures should be unaffected. Confirm that their `validation.tsv`
status columns are unchanged apart from the new columns.

---

## 6. Acceptance criteria

- [ ] `INT-ANIMAL-01` ships 6 barcodes: 5 `pass` + COX1 `partial`,
      ≈1.4 kb, stop-free under the recorded table.
- [ ] Every `pass` record in `barcodes.fasta` is byte-identical to the
      pre-task output. Verify by diffing the three integration fixtures'
      `barcodes.fasta` with partial records excluded.
- [ ] `params.barcode_partial_min_nt = 0` reproduces pre-task
      `validation.tsv` status/reason columns exactly.
- [ ] A partial is visibly distinct from a pass in:
  - [ ] `validation.tsv`
  - [ ] `barcodes.fasta`
  - [ ] `metadata.json`
  - [ ] per-sample report (Barcodes tab and key findings)
  - [ ] run report
- [ ] A partial is never listed as a dropout and never causes
      `no_barcode`.
- [ ] Spec and `assertions.sh` state the relaxed coherence invariant
      (§3.3), and CI enforces it.
- [ ] flake8 is clean on touched `bin/*.py` and `scripts/tests/*.py`.
      `scripts/pytest.sh` passes.
- [ ] A fresh `INT-ANIMAL-01` report has been rendered, and its path is
      given to the user for inspection.

---

## 7. Out of scope

- Correcting frameshifts or stop codons in the emitted sequence, e.g.
  by consensus re-polishing or reference-guided base fixing. We ship
  the assembly verbatim, or a verbatim sub-span of it.
- Rescuing loci that fail `identity_below_floor` or `invalid_length`.
  The identity-floor question stays in `tasks/todo.md`.
- Trimming over-extended miniprot spans at frameshift boundaries for
  annotation (separate todo.md item under *Annotation*).
- Multi-exon / trans-spliced partials (`plant_mt` `nad1`). These are
  failed explicitly (§3.2) until a `plant_mt` fixture reaches this
  stage.
- Taxodactyl hand-off. **Resolved (2026-09-29):** Taxodactyl accepts
  partial sequences down to 20 nt, so partials ship in the same
  `barcodes.fasta` with the ID format unchanged (§3.4). The 20 nt figure
  is enforced only as a lower bound on the param (§3.4, *Why 100*), not
  as the default.

---

## Outcomes

(fill in on completion)
