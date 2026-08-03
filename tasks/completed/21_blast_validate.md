# Task 21 — P3 stage 11: `BLAST_VALIDATE`

**Phase:** P3 (from [spec §6](../spec/06-phases.md)).
**Goal:** Replace the P0 stub for stage 11 with a real `blastn`
invocation that sanity-checks each `BIN_TARGET`-selected contig against
the target-appropriate RefSeq organelle BLAST DB. Output is a
per-contig top-hits TSV consumed by `COLLATE` (P4) for inclusion in
the per-sample report and by `RUN_REPORT` for cross-sample summaries.
This is a thin biocontainer wrapper — no bespoke Python, no
[C1–C7 custom logic](../spec/02-stages.md#22-custom-logic-components).

Per [spec §2 stage 11](../spec/02-stages.md#2-stage-detail) the
kingdom-appropriate DB is selected per-sample from the
[reference bundle](../spec/04-reference-data.md#42-refseq-organelle-dbs-for-blast_validate):

| `meta.assembly_target` | BLAST DB (relative to `params.blast_db`) |
|---|---|
| `plant_pt`   | `refseq_pt`             |
| `plant_mt`   | `refseq_mt_viridiplantae` |
| `animal_mt`  | `refseq_mt_metazoa`     |

**Prerequisite:** [task 18 — BIN_TARGET](completed/18_bin_target.md) is real, so
`BIN_TARGET.out.binned` carries a real `target.fasta` (possibly empty —
see the empty-query note in §2) into this stage. Channel wiring already
exists in [main.nf](../main.nf) —
`BLAST_VALIDATE(BIN_TARGET.out.binned)` — and requires no topology
change.

**Note on `plant_pt` input:** on `plant_pt` samples where
[task 20 (C4 plastid canonicalisation implementation)](completed/20_plastid_canonicalise.md)
has landed, `target.fasta` is the canonicalised `path1` sequence
emitted by C4 (LSC + IR + SSC + rc(IR), ~150 kb). This is transparent
to BLAST_VALIDATE — the process still consumes `target.fasta`
unchanged and requires no plant_pt-specific branch. All
plastid-quadripartite awareness stays inside `bin/bin_target.py` +
`bin/plastid_canonicalise.py` per
[spec §3.6 step 5](../spec/03-organelles.md#36-plastid-quadripartite-canonicalisation).

**Exit criteria:**

- `-profile integration` produces, for every sample that reaches
  `BIN_TARGET`:
    - `results/<sample_id>/blast_validate/<sample_id>.blast.tsv` — a
      real BLAST outfmt-6 TSV. For each `target.fasta` record, ≥ 1
      hit row with columns
      `qaccver saccver pident length qcovs evalue bitscore stitle`.
      Empty file if `target.fasta` is empty (documented "no hits"
      case per [spec principle 7](../CONSTITUTION.md)).
- BLAST_VALIDATE runs in the pinned BLAST+ biocontainer
  `quay.io/biocontainers/blast:2.17.0--h66d330f_0` (SHA-pinned in
  [`conf/containers.config`](../conf/containers.config)) — no host
  tools.
- Integration assertion block added to
  [`tests/integration/assertions.sh`](../tests/integration/assertions.sh)
  verifying per-sample top-hit consistency over `ASSEMBLING_SAMPLES`:
  the top `saccver`'s `stitle` contains a substring appropriate for
  the sample (`Acyrthosiphon` for `INT-ANIMAL-01`, plastid keywords
  for `INT-PLANT-01-pt`). See §4 for why `INT-PLANT-01-mt` is not
  asserted.
- `-profile stub -stub-run` still green (stub block unchanged).
- `bin/` unchanged, `scripts/tests/` unchanged — off-the-shelf tool
  wrapper per [spec §5a rule of thumb](../spec/05-test-data.md#5a-tests).

**Not in scope:**

- **Any decision logic based on BLAST output** (e.g. failing a sample
  when identity is low). BLAST_VALIDATE emits *data* for the report;
  interpretation lives in `COLLATE` / `RUN_REPORT` (P4) where the
  three-way negative-clarity classification
  ([spec principle 7](../CONSTITUTION.md)) is centralised. Adding
  any threshold here would be a second, competing decision surface.
- **Reference bundle build.** The RefSeq organelle BLAST DBs are
  built by [`scripts/refdata/build_validate.sh`](../scripts/refdata/build_validate.sh)
  per [spec §4.2](../spec/04-reference-data.md#42-refseq-organelle-dbs-for-blast_validate).
  They already exist: `refs/v2026.08/validate/` (RefSeq release 236)
  carries all three stems and is inside the `refs-v2026.08.tar.gz`
  bundle that CI fetches. Nothing to build here — see §4.
- Alignment against a nucleotide-collection (nt) DB. Spec pins the
  target-specific curated DBs — noise reduction and provenance
  clarity trump broad recall here.
- Custom output format. `outfmt 6` with the columns above is what
  the report renderer expects; SAM/XML would add parsing surface
  without value.

**Cross-cutting rules (from [spec §1a](../spec/01-pipeline-flow.md#1a-engineering-constraints)):**

- Container pinned by SHA, never `latest`.
- No host tools; `blastn` comes from the container.
- BLAST DB path staged via `${file(params.blast_db)}` in the script
  block — acceptable inline `file(params.*)` pattern per
  [spec §1a](../spec/01-pipeline-flow.md#1a-engineering-constraints)
  (workflow-wide reference, not per-sample fanout).

---

## 1. Container pin — `conf/containers.config`

Replace the `python:3.12-slim` stub + `TODO P3` comment for
`BLAST_VALIDATE` in [`conf/containers.config`](../conf/containers.config)
with the pinned BLAST+ biocontainer already staged in the TODO
comment:

```groovy
withName: 'BLAST_VALIDATE' {
    container = 'quay.io/biocontainers/blast:2.17.0--h66d330f_0@sha256:<pinned>'
}
```

**Resolving the SHA:** as in
[task 16 §2](completed/16_metaflye.md#2-container-pin-confcontainersconfig):

```bash
docker pull quay.io/biocontainers/blast:2.17.0--h66d330f_0
docker inspect --format='{{index .RepoDigests 0}}' \
    quay.io/biocontainers/blast:2.17.0--h66d330f_0
```

Record the resolved digest in Outcomes on completion. Single-tool
biocontainer — no `mulled-build`.

## 2. `modules/local/blast_validate.nf` (replace stub)

Replace the current stub `script:` block in
[`modules/local/blast_validate.nf`](../modules/local/blast_validate.nf)
with a real `blastn` invocation. Retain the stub block unchanged for
`-profile stub`.

Pseudocode of the intended shape:

```groovy
script:
def db_name = [
    animal_mt: 'refseq_mt_metazoa',
    plant_pt:  'refseq_pt',
    plant_mt:  'refseq_mt_viridiplantae',
][meta.assembly_target]
"""
if [[ -s ${target_fasta} ]]; then
    blastn \\
        -query ${target_fasta} \\
        -db ${file(params.blast_db)}/${db_name} \\
        -outfmt '6 qaccver saccver pident length qcovs evalue bitscore stitle' \\
        -max_target_seqs 5 \\
        -num_threads ${task.cpus} \\
        -out ${meta.sample_id}.blast.tsv
else
    # Empty target.fasta from BIN_TARGET → empty BLAST output; not an error.
    : > ${meta.sample_id}.blast.tsv
fi
"""
```

Notes:

- **Output shape unchanged.** Module already emits
  `tuple(meta, target_fasta, blast_tsv)` — matches what `ANNOTATE`,
  `MINIPROT_EXTRACT` and the `COLLATE` join consume in
  [main.nf](../main.nf). No `output:` change. The module's `input:`
  already matches `BIN_TARGET.out.binned`'s
  `tuple(meta, target.fasta, secondaries.tsv)`; `secondaries.tsv` is
  staged but unused by `blastn`.
- **Guard on empty `target.fasta`.** Two distinct upstream states
  produce an empty query, and both are legitimate soft-fails rather
  than errors:
    1. BIN_TARGET selected nothing
       ([task 18 §11](completed/18_bin_target.md#11-notes-non-issues));
    2. a `plant_pt` sample where C4 withheld the substitution —
       `substitution_applied: false` with
       `substitution_withheld_reason: "no_c3_selection"` — which
       emits an empty `target.fasta` alongside a populated
       `plastid_isoforms/`
       ([task 24 §3.2](completed/24_plastid_substitution_guard.md)).

  `blastn` on an empty query exits non-zero with a confusing
  message. The empty-query short-circuit keeps BLAST_VALIDATE's
  failure mode aligned with the rest of the soft-fail chain — an
  empty TSV, not a Nextflow error.
- **`-max_target_seqs 5`.** Small: the TSV is a sanity check, not a
  full similarity search. 5 hits per contig is enough to spot a
  mis-binned contig without blowing up the report's table.
- **`stitle` in output.** Included so the assertion (§4) and the
  report can show human-readable species names without a second
  `blastdbcmd` lookup.
- **`process_medium` label.** BLAST on an organelle-scale query
  against a target-specific DB fits in the P3 medium budget; escalate
  only if `plant_mt` fixtures exceed it.
- **Runtime watch on `plant_pt`.** `refseq_pt` is the largest DB in
  the bundle (582 MB `.nsq`, 15 233 genomes) and the integration
  profile caps at 2 vCPU / 14 GB. A ~150 kb canonicalised `path1`
  query under default megablast should be minutes, not tens of
  minutes, but this is the first stage to touch a DB that size —
  record the observed wall time in Outcomes and check it against the
  60-min CI budget.
- **BLAST DB path.** `${file(params.blast_db)}/${db_name}` — the
  bundle-root is the value channel; the per-target subname is derived
  from `meta.assembly_target`. This matches the reference bundle
  layout defined in
  [spec §4.4](../spec/04-reference-data.md#44-consolidated-build-script).

## 3. Fast CI

No change. Stub block still `touch`es the TSV; container-coverage
lint stays satisfied.

## 4. Integration-test wiring

Add a BLAST_VALIDATE assertion block after the BIN_TARGET block
(added in [task 18](completed/18_bin_target.md)) in
[`tests/integration/assertions.sh`](../tests/integration/assertions.sh):

```bash
# BLAST_VALIDATE is real (task 21). Asserted over ASSEMBLING_SAMPLES,
# not SAMPLES: a sample that soft-fails the coverage gate never
# reaches BIN_TARGET and so produces no blast.tsv at all.
declare -A EXPECTED_TITLE_SUBSTR=(
    [INT-ANIMAL-01]="Acyrthosiphon"
    [INT-PLANT-01-pt]="plastid|chloroplast"
)
for sample in "${ASSEMBLING_SAMPLES[@]}"; do
    tsv="$OUTDIR/$sample/blast_validate/${sample}.blast.tsv"
    tgt="$OUTDIR/$sample/bin_target/target.fasta"
    if [[ ! -e "$tsv" ]]; then
        echo "FAIL: $sample blast.tsv missing"
        FAILED=1
        continue
    fi
    # An empty target.fasta (no C3 selection, or a withheld C4
    # substitution) legitimately yields an empty TSV — only demand
    # hits when there was actually a query.
    if [[ ! -s "$tgt" ]]; then
        echo "SKIP: $sample target.fasta empty — no BLAST query expected"
        continue
    fi
    if [[ ! -s "$tsv" ]]; then
        echo "FAIL: $sample blast.tsv empty but target.fasta is non-empty"
        FAILED=1
        continue
    fi
    # Column 8 is stitle. Grep is line-wise; case-insensitive substr match.
    if ! head -n 5 "$tsv" | cut -f8 \
            | grep -qEi "${EXPECTED_TITLE_SUBSTR[$sample]}"; then
        echo "FAIL: $sample top-5 hits do not mention '${EXPECTED_TITLE_SUBSTR[$sample]}'"
        FAILED=1
    else
        echo "OK:   $sample BLAST top hits look on-target"
    fi
done
```

**Why `INT-PLANT-01-mt` is absent.** Since
[task 25](completed/25_coverage_gate_carryover.md) the fixture
soft-fails C2 (28.48× mitochondrial depth against a 30× MIN after
[task 28](completed/28_plastid_masked_mt_panel.md)) and never reaches
`BIN_TARGET`, so there is no `target.fasta` to validate. The
consequence to accept knowingly: the `plant_mt` →
`refseq_mt_viridiplantae` branch of the §2 DB map ships **unexercised
by integration**. It is a one-line map entry with no target-specific
logic behind it, and the fixture replacement tracked in
[`tasks/todo.md`](todo.md) will cover it. Do not add a fourth fixture
or relax the gate for this task's benefit.

**Progressive-assertion note.** No hard identity or length threshold
is enforced here — spec §2 stage 11 is a "sanity check", and the
per-sample report will render the raw numbers. A substring assertion
on `stitle` catches the failure mode that matters (contig mis-binned
to the wrong kingdom / organelle) without introducing brittleness to
per-fixture identity drift.

Also update the progressive-uncomment header at the top of
`assertions.sh` to tick off the `VALIDATE (real BLAST) — Taxonomic-
identification consistency checks` row.

**No new expected-output fixture files** — the assertion is structural
(non-empty TSV) plus a substring check on the sample-specific title.
Adding a full expected-hits table would be a per-fixture drift trap
that violates [constitution rule 19](../CONSTITUTION.md).

**BLAST DB availability in CI — already satisfied; nothing to do.**
The integration profile points `params.blast_db` at
`${projectDir}/refs/v2026.08/validate`
([`conf/integration.config`](../conf/integration.config)), and all
three stems (`refseq_pt`, `refseq_mt_viridiplantae`,
`refseq_mt_metazoa`, RefSeq release 236) are present in the local
`refs/v2026.08/validate/` tree and inside `refs-v2026.08.tar.gz`.
The integration workflow stages them via
`bash scripts/fetch_refs.sh v2026.08`
([`.github/workflows/integration.yml`](../.github/workflows/integration.yml)),
which pulls the bundle from Azure blob and verifies it against
[`scripts/refs-v2026.08.sha256`](../scripts/refs-v2026.08.sha256).

The bundle is uploaded and public: `v2026.08/refs.tar.gz` in the
`refdata-wf5` container returns `HTTP 200` at 1 213 060 376 bytes,
matching the locally built tarball. A clean runner will fetch and
checksum it without intervention — no refs work is in this task's
path.

## 5. Unit tests

**None required.** Off-the-shelf tool wrapper with no
[C-component custom logic](../spec/02-stages.md#22-custom-logic-components).
Per [spec §5a](../spec/05-test-data.md#5a-tests):

- Channel wiring covered by `-stub-run` (fast CI on every push).
- Real-tool behaviour covered by nightly `-profile integration` +
  `assertions.sh`.

If a `bin/parse_blast_hits.py` helper is added later (e.g. to
aggregate best-hit identity per contig for the report), that would
be a C-component with its own unit tests — not this task.

## 6. Verification

1. `nextflow run . -profile stub -stub-run` — 17/17 processes hit
   their stub, no regressions in `tests/output/STUB-01/`.
2. `nextflow run . -profile integration` locally (requires the
   `refs/v2026.08/validate/` DB tree present per §4):
   - `INT-ANIMAL-01/blast_validate/INT-ANIMAL-01.blast.tsv` — ≥ 1
     row; top hit `stitle` mentions `Acyrthosiphon pisum` or a
     close aphid relative.
   - `INT-PLANT-01-pt/…/.blast.tsv` — ≥ 1 row; top hits mention
     `plastid` or `chloroplast`.
   - `INT-PLANT-01-mt` — no `blast_validate/` output expected; the
     sample soft-fails C2 (§4). If a stale TSV is present from an
     earlier run, clear the outdir — `publishDir` never deletes (see
     the header note in `assertions.sh`).
3. `bash tests/integration/assertions.sh` — new BLAST_VALIDATE block
   reports OK for all three samples; pre-existing blocks still green.
4. `flake8` — no Python touched; nothing to run.
5. CI green on `Tests` (fast) + nightly `Integration` on the PR.

## 7. Deliverables checklist

- [x] [`conf/containers.config`](../conf/containers.config) — SHA-pinned
      BLAST+ biocontainer for `BLAST_VALIDATE`, `TODO P3` comment
      removed.
- [x] [`modules/local/blast_validate.nf`](../modules/local/blast_validate.nf) —
      real `blastn` `script:` block with empty-query guard;
      per-target DB name derived from `meta.assembly_target`; stub
      retained.
- [x] [`tests/integration/assertions.sh`](../tests/integration/assertions.sh) —
      BLAST_VALIDATE substring-check block added over
      `ASSEMBLING_SAMPLES`, with the empty-`target.fasta` skip;
      progressive-uncomment header updated.
- [x] Fast CI (`Tests`) + `Integration` workflows both green locally;
      to be confirmed on the PR.

## 8. Notes / non-issues

- **No channel-topology change.** [main.nf](../main.nf) already wires
  `BLAST_VALIDATE(BIN_TARGET.out.binned)`; the emit tuple shape is
  identical to the stub's; downstream
  `ANNOTATE(BLAST_VALIDATE.out.validated)` and
  `MINIPROT_EXTRACT(BLAST_VALIDATE.out.validated)` still consume
  `(meta, target_fasta, blast_tsv)`.
- **Why no decision here.** Centralising the low-confidence /
  no-recovery classification in `COLLATE` keeps the three-signal
  negative-clarity contract
  ([spec principle 7](../CONSTITUTION.md)) auditable in one place;
  scattering pass/fail thresholds across stages fragments that.
- **Guard is a `bash if`, not a Groovy `when:`.** A Groovy `when:`
  would silently skip the process and starve the downstream join
  channel; the `if` inside the script emits an empty TSV that
  matches the process's declared output. Downstream consumers see
  a valid (empty) file and can handle it.
- **BLAST database name construction.** Written as a Groovy map
  inside the script block rather than a `params.blast_db_names`
  table because the mapping is fixed by
  [spec §4.2](../spec/04-reference-data.md#42-refseq-organelle-dbs-for-blast_validate)
  and never varies across profiles. Promoting to `params` would
  add configuration surface for no downstream flexibility gain.
- **`-max_target_seqs 5` and BLAST's caveat.** Recent BLAST releases
  warn that `-max_target_seqs` is applied post-search; this is fine
  here because the search space is a curated single-organelle DB and
  we want a small shortlist for the report, not the "true top N".

## 9. Outcomes

- Container:
  `quay.io/biocontainers/blast:2.17.0--h66d330f_0@sha256:78e8ebf231ce498b0015100a891693ef2acbbf7df60b33299517ca2b4891e567`.
- Observed top hits per sample (`-profile integration`, local run):
    - `INT-ANIMAL-01`: top hit `NC_011594.1`, "Acyrthosiphon pisum
      mitochondrion, complete genome" at 98.28 % identity.
    - `INT-PLANT-01-pt`: top hit `NC_018117.1`, "Datura stramonium
      chloroplast, complete genome" at 99.60 % identity.
    - `INT-PLANT-01-mt`: no `blast_validate/` output — sample
      soft-fails the coverage gate before `BIN_TARGET` (per §4), as
      expected.
- BLAST_VALIDATE wall time per sample (`nextflow log … realtime`):
    - `INT-PLANT-01-pt` (against `refseq_pt`, the 582 MB DB): 40.1 s —
      well inside the 60-min CI budget, no runtime concern.
    - `INT-ANIMAL-01` (against `refseq_mt_metazoa`): 587 ms.
- No threshold surprises — both on-target samples returned an
  expected-species top hit at ≥ 98 % identity; no decision logic was
  added per the task's scope.
- Full `-profile integration` run: 40/40 processes succeeded in
  10 m 41 s (0.4 CPU-hours). `bash tests/integration/assertions.sh`
  green, including the new BLAST_VALIDATE block.
- `-profile stub -stub-run`: 17/17 processes green, no regressions.
