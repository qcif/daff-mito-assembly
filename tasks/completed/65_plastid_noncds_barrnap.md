# Task 65 — `plant_pt` tRNA and rRNA annotation with barrnap

**Phase:** P4b (plant annotation arms). Implements the decision in
task 32_research_plastid_noncds.md §7.4: **option A, barrnap alone,
`plant_pt` only** (decided by the user, 2026-10-13).

**Status:** open. **No preconditions.** Independent of task
63_cds_locus_clustering.md and task 64_protein_panel_label_check.md.
Should land before task 66_plastid_map_structure.md, so that the map's
visual check sees real tRNA/rRNA features with normalised names.

## Overview

A plastid (chloroplast) genome carries about 80 protein-coding genes,
about 30 tRNA genes and 4 rRNA genes. The rRNA genes sit in the
inverted repeat, so there are 8 rRNA loci in all. The pipeline already
annotates the protein-coding genes with miniprot (stage 12). miniprot
aligns proteins, so it cannot find RNA genes. A `plant_pt` sample
today ships a CDS-only annotation, labelled `status: "ok_cds_only"`.

Task 32 measured the options. **barrnap 1.10** (Rfam covariance models
run with Infernal) found on the plastid fixture:

- **29 of 37 tRNAs, with zero false calls.** The 8 it misses are all
  tRNAs that contain an intron (`trnK-UUU`, `trnG-UCC`, `trnL-UAA`,
  `trnV-UAC`, `trnI-GAU` ×2, `trnA-UGC` ×2).
- **6 of 8 rRNA loci**, both inverted-repeat copies of `rrn16`,
  `rrn23` and `rrn5`. It misses `rrn4.5` in both copies, because
  Rfam's bacterial model treats 4.5S as part of 23S, so `rrn23` comes
  out ~94 bp too long at its 3′ end.

One bioconda container runs it in about 6 seconds. The alternative
(tRNAscan-SE plus barrnap) gave the same correct calls at 14× the run
time. Chloë would be better but is unlicensed (task 32 §1, §7.3).

This task adds barrnap to stage 13a for `plant_pt`. tRNA and rRNA
features join the published GFF under plastid gene names (`trnM-CAU`,
`rrn16`). The status becomes `ok`. The gaps barrnap is known to leave
are **written into `annotation_summary.json` and shown in the
report**, so a reader is never left to infer them from what's missing
(principle 7). `plant_mt` stays CDS-only. It needs a NUPT cross-check
first (shelved task 37).

## 1. Constraints

- **Principle 8.** A barrnap failure is a status, never a run abort or
  a lost sample. **Do not rely on `errorStrategy = 'ignore'`**: an
  ignored task emits nothing, and the downstream `join` silently drops
  the sample (`tasks/todo.md`, Intake robustness, same failure mode as
  `NANOPLOT_RAW`). The barrnap step must exit 0 and pass its real exit
  code downstream as data, the way `ANNOTATE` already does for MITOS2.
- **Task 41 semantics.** A configured annotator that exits non-zero,
  or exits 0 with zero features, gives `status: "annotator_failed"`
  with a `reason` that names the tool and the exit code. Never
  `ok_cds_only`.
- **Rule 11/12.** barrnap runs in its own pinned bioconda image, not a
  mulled bundle. The `ANNOTATE` process stays in the MITOS2 image. No
  bespoke Dockerfile.
- **Rule 13.** barrnap's GFF reaches `ANNOTATE` as a staged `path`.
  Samples with no barrnap step still need a stageable placeholder, not
  a string. Use a **file**, never an empty directory (`tasks/todo.md`,
  Azure: empty directories vanish on blob storage).
- **CDS features unchanged** (task 31 §4 point 1). barrnap adds
  features and never touches CDS rows.
- **Raw names kept.** Normalised names go in `Name=`. barrnap's own
  `Name=` and `product=` stay in the attributes (task 31 §4 point 4
  does the same for MITOS2).

## 2. Container

- `quay.io/biocontainers/barrnap:1.10.6--pl5321hdfd78af_0`. Pin it by
  `sha256` digest in `conf/containers.config`, like every other image,
  with the tag in a trailing comment.
- Known quirks (task 32 §7.1), each to be re-checked and recorded:
  - `barrnap --version` prints **1.10.5** in this 1.10.6 build. Record
    both the image tag and the self-reported version in provenance,
    and say which is which.
  - `--all` crashes (broken `transterm` data path). Use
    `--rrna --trna`.
  - There are no `mito`/`euk` kingdoms in 1.10. Use `--kingdom bac`:
    plastid rRNAs are bacterial-type.
  - Peak RSS is **3.2 GB** (Rfam CM database load). Set the process
    `memory` directive explicitly (e.g. 6 GB) rather than inheriting a
    label that may be smaller. Check it against `resourceLimits` in
    `conf/integration.config` and the Azure profile.

## 3. Pipeline changes

### 3.1 Config

`nextflow.config`, `params.annotate`:

```
plant_pt: [ non_cds_tool: 'barrnap', kingdom: 'bac' ],   // task 65
plant_mt: [ non_cds_tool: 'none' ],                      // shelved task 37
```

Fix the stale `// task 33` comment on `plant_mt` to name shelved task
37, and do the same for the `"see task 32 (plant_pt) / task 33
(plant_mt)"` reason string in `bin/annotate_summary.py`.

### 3.2 New process `BARRNAP` (stage 13a, before `ANNOTATE`)

- Input: `(meta, target_fasta)` from the same `ch_cds` the other
  stage 13a inputs use. Runs only for samples whose target config says
  `non_cds_tool == 'barrnap'` (branch the channel; don't gate inside
  the script).
- An empty `target_fasta` (BIN_TARGET withheld a substitution, task 24
  §3.2) gives an empty GFF and exit code 0, the same as
  `MINIPROT_CDS` already does.
- Outputs, all files:
  - `barrnap.gff`;
  - `barrnap.log`: stderr, plus a trailing `--- barrnap exit code: N`
    line (the `mitos.log` pattern), published for audit;
  - `barrnap_exit.txt`: the exit code alone, read by `ANNOTATE`.
- `set +e` around the barrnap call, then capture the code and exit 0.
- Stub: empty `barrnap.gff`, `barrnap_exit.txt` containing `0`, and a
  log line.

### 3.3 `ANNOTATE` wiring

- Add one `path(non_cds_gff)` and one `path(non_cds_exit)` to the
  input tuple. For non-barrnap samples, join a shared placeholder file
  (e.g. `assets/NO_FILE`, the nf-core idiom) so the tuple shape is
  uniform and nothing is passed as a string.
- In the script block, when `non_cds_tool == 'barrnap'`, pass
  `--barrnap-gff`, `--annotator-exit` (from the exit file) and the
  barrnap version to `annotate_summary.py`, mirroring the MITOS2
  branch. Keep the `script:` block under ~10 lines of new bash
  (rule 14). If it grows, the logic belongs in C8.

### 3.4 C8 (`bin/annotate_summary.py`)

Stays stdlib-only (it runs in the MITOS2 image).

1. **Parse** barrnap's GFF3. Keep `tRNA` and `rRNA` rows. barrnap
   writes the input FASTA's own record IDs in column 1, so unlike
   MITOS2 no `seqid` remap is needed. Assert this in a test.
2. **Normalise names** to plastid nomenclature. Observed barrnap
   output on the fixture (2026-10-13):

   | barrnap `Name=` | `product=` | Normalised `Name=` |
   |---|---|---|
   | `tRNA-Leu` | `transfer RNA (caa)` | `trnL-CAA` |
   | `tRNA-Met` | `transfer RNA (cat)` | `trnM-CAU` (see point 3) |
   | `16S_rRNA` | `16S ribosomal RNA` | `rrn16` |
   | `23S_rRNA` | `23S ribosomal RNA` | `rrn23` |
   | `5S_rRNA` | `5S ribosomal RNA` | `rrn5` |

   - tRNA: one-letter amino-acid code from the three-letter name, then
     `-`, then the anticodon from `product=`, upper-cased with T → U.
     The 20 standard amino acids plus `SeC` go in a small constant
     table.
   - rRNA: a constant table keyed on barrnap's `Name=`.
   - **Unrecognised names** (no table entry, no anticodon) keep
     barrnap's raw name, are counted in the summary
     (`non_cds_unnormalised`), and print a warning to stderr. Never
     drop the feature.
3. **The CAU anticodon is ambiguous.** In plastids `trnfM-CAU`,
   `trnM-CAU` and `trnI-CAU` (lysidine-modified, reads AUA as Ile) all
   carry the CAU anticodon, and barrnap names all of them Met. Emit
   `trnM-CAU` and add an attribute saying so (e.g.
   `Note=CAU isotype ambiguous: fMet/Met/Ile`). Do not guess `trnI`
   from position. Rule 18: say it's ambiguous rather than state a
   false precision.
4. **rRNA completeness:** compare normalised rRNA names against the
   `rrna` list for `plant_pt` in `assets/organelle_gene_sets.json`
   (`rrn16`, `rrn23`, `rrn4.5`, `rrn5`). Record found and missing in
   the summary. `rrn4.5` is expected to be missing. tRNA count against
   `trna_count_approx` is descriptive only. Do not compute a tRNA
   "completeness" fraction from an approximate denominator.
5. **Status:** follow the existing table. `ok` when features merged;
   `annotator_failed` on a non-zero exit or zero tRNA/rRNA features
   from a successful run (task 41); `no_assembly` when `cds.gff` is
   empty (unchanged).
6. **Known gaps:** a new summary field (e.g.
   `non_cds_known_gaps`), filled from a constant keyed by
   `non_cds_tool`. For `barrnap` it lists what task 32 measured:
   intron-containing tRNAs not called (name the 8), `rrn4.5` not
   called (its sequence sits inside the over-extended `rrn23`), and
   the CAU isotype ambiguity. Add the CDS-side limits from task 32
   §7.2 that apply to every `plant_pt` sample: the `rps12` 5′ exon is
   not called, and the micro-exons of `petB`, `petD` and `rpl16` are
   missed. Each entry gets a short code plus one sentence, so the
   report can render it and an auditor can read it. Empty list or
   `null` for `mitos2` (no gaps measured yet). State which in the
   schema.
7. **Provenance:** `non_cds_source: "barrnap"`;
   `tool_versions.barrnap` holds the image tag and self-reported
   version (§2); `reference_data` names the Rfam models barrnap ships
   (whatever barrnap itself reports, not a guess).
   `cds_crosscheck` and `cds_rescued` are `null`: barrnap calls no
   CDS.
8. **Schema:** the new fields are a contract change. Bump `$schema`
   in C8, and in C11 (`bin/annotation_scores.py`, which rewrites it),
   following the project's existing v2 → v3 convention, and update
   every reader (`bin/collate.py`, the report, tests).

### 3.5 Report and map

- Per-sample report, Assembly tab, *Genome annotation*: when
  `non_cds_known_gaps` is non-empty, render it as a short "Known
  limitations of this annotation" list beside the feature summary.
  It is informational (`info` style), never a warning or danger badge.
  These are known gaps of the method, not problems with the sample
  (spec §6a, task 46 §3.2's icon discipline).
- The text that currently says plant arms are CDS-only must now say
  so only for `plant_mt`.
- Run report: no change expected. Confirm the annotation status
  column shows `ok` for `plant_pt`.
- `ORGANELLE_MAP`: no code change. `parse_non_cds_features` already
  draws any `tRNA`/`rRNA` row with `Name=` (task 32 §8.1). Confirm the
  normalised names appear in tooltips.

## 4. Pseudocode (C8 additions)

```
AA3_TO_1 = {"Ala": "A", ..., "Met": "M", ..., "SeC": "U"}
RRNA_NAMES = {"16S_rRNA": "rrn16", "23S_rRNA": "rrn23", "5S_rRNA": "rrn5"}

for row in barrnap_gff where type in (tRNA, rRNA):
    if type == tRNA and anticodon := re("\((\w{3})\)", product):
        name = f"trn{AA3_TO_1[aa3]}-{anticodon.upper().replace('T','U')}"
        if anticodon == CAU: note = "CAU isotype ambiguous: fMet/Met/Ile"
    elif type == rRNA and raw in RRNA_NAMES:
        name = RRNA_NAMES[raw]
    else:
        name = raw; unnormalised += 1; warn
    emit row with Name=name, keeping raw Name/product as attributes
```

## 5. Testing

### 5.1 Unit tests (`scripts/tests/test_annotate_summary.py`)

Hold C8 at 100 % branch coverage. Fixtures are small barrnap GFF
strings in `tmp_path`. Copy the header and attribute layout verbatim
from the real barrnap output above, so the tests pin the real format.

- Name normalisation: every row in the §3.4 table; DNA → RNA
  anticodon; lower → upper case; the CAU note; an unknown amino acid;
  a tRNA with no anticodon in `product`; an unknown rRNA name. The
  last three keep the raw name and are counted.
- Raw `Name`/`product` preserved alongside the normalised `Name`.
- Status: `ok` with features; `annotator_failed` on non-zero exit
  (reason names barrnap and the code); `annotator_failed` on exit 0
  with zero features; `no_assembly` on empty `cds.gff` regardless of
  barrnap.
- rRNA found/missing against the gene set (`rrn4.5` missing).
- `non_cds_known_gaps` present for barrnap, empty or `null` for
  mitos2, absent or `null` for `none`. Matches the schema.
- `seqid` passes through unchanged.
- CDS rows byte-identical with and without barrnap input.
- `animal_mt` (MITOS2) and `plant_mt` (`none`) outputs unchanged.
  Existing tests must pass without edits beyond the schema-version
  bump.
- C11 test: schema bump only.

### 5.2 Stub run

The new process and the placeholder join must pass `-stub-run`, and
CI's container-coverage and bare-params lints. The `INT-PLANT-01-pt`
stub sample must reach `COLLATE`.

### 5.3 Integration

Update the `INT-PLANT-01-pt` block in `tests/integration/assertions.sh`:

- status `ok`, `non_cds_source == "barrnap"` (today it asserts
  `ok_cds_only` and `null`, and has a task 41 tripwire that becomes
  live);
- extend `tests/integration/expected/plant_pt/annotation_bounds.json`
  (created by task 63; create it in that shape if task 63 hasn't
  landed) with `min_trna`, `min_rrna`, `required_rrna`
  (`rrn16`, `rrn23`, `rrn5`), and `rrn16`/`rrn23`/`rrn5` in the IR-pair
  list (at least 2 each). Put the measured figures in its `_comment`,
  as the `animal_mt` file does. The rRNA pairs are biology. The tRNA
  floor is a correctness floor, not a tuned threshold (spec §5.1);
- `non_cds_known_gaps` non-empty;
- the GFF has no `tRNA` row whose `Name=` still starts with `tRNA-`
  (normalisation actually ran);
- `INT-PLANT-01-mt` stays `ok_cds_only` (if the fixture reaches
  annotation, which it does under the two-floor gate).

Record wall time and peak RSS of `BARRNAP` on the fixture in Outcomes,
against task 32 §7.1's 6 s / 3.2 GB.

### 5.4 Report check

Render `INT-PLANT-01-pt`'s report. Check that the known-limitations
list, the `ok` status and tRNA/rRNA tooltips on the map all show.
Give the report path to the user (CLAUDE.md).

## 6. Spec updates

- spec/07-open-questions.md §8 item 3: `plant_pt` resolved. Record
  barrnap with its evidence (task 32 §7.1 figures), the known gaps,
  and the rejected alternatives with one-line reasons (GeSeq: web
  only; Chloë: no licence, angiosperm-optimised, kept as a possible
  upgrade; PGA: stale; tRNAscan-SE + barrnap: same recovery, 14× the
  cost). Point to task 32_research_plastid_noncds.md by plain name.
  `plant_mt` stays open.
- spec/07-open-questions.md §9 item 12: note that the §3 thresholds
  from task 32 were not met (78.4 % strict tRNA recovery) and the arm
  shipped anyway, with the gap recorded as data.
- spec/03-organelles.md §3.1 table, `ANNOTATE` row: barrnap for
  `plant_pt`.
- spec/02-stages.md: stage 13a row (new `BARRNAP` step, container,
  memory) and C8 row (barrnap parsing, normalisation, known gaps).
- spec/01-pipeline-flow.md stage diagram: the 13a box text.
- spec/06-phases.md P4b: `plant_pt` done, `plant_mt` open.
- spec/06a-reports.md Assembly tab: the known-limitations list, and
  CDS-only wording narrowed to `plant_mt`.
- spec/containers.md: add the barrnap image.
- spec/brief.md §3.6 names GeSeq. Task 32 flagged that the brief
  should be reconciled. **Ask the user before editing the brief**: it
  is the client-facing document.

## 7. Out of scope

- `plant_mt` non-CDS (shelved task 37). Recording that barrnap would
  likely serve it, subject to a NUPT cross-check, is enough.
- Intron-containing tRNAs and `rrn4.5`. These are the recorded gaps.
  Closing them is Chloë's job if its licence clears (task 32 §7.3).
- tRNAscan-SE as a second opinion.
- Gene labels on the map, and the segment ring (task
  66_plastid_map_structure.md).

## 8. Checklist

- [ ] §2 image pinned by digest; quirks re-verified
- [ ] §3.1–3.3 config, `BARRNAP` process, `ANNOTATE` wiring
- [ ] §3.4 C8 parsing, normalisation, status, known gaps, provenance,
      schema bump
- [ ] §3.5 report renders known gaps; map tooltips checked
- [ ] §5.1 unit tests; 100 % branch on C8; flake8 clean
- [ ] §5.2 stub run and CI lints pass
- [ ] §5.3 integration assertions updated and passing; RSS/time recorded
- [ ] §5.4 report rendered; path given to the user
- [ ] §6 spec updated (brief only with the user's go-ahead)
- [ ] task 32: §7.4 decision recorded (done when this brief was
      written); move task 32 to `tasks/completed/` once the Chloë
      email (§7.3) is sent or explicitly dropped. That is the user's
      call
- [ ] Outcomes section written
