## 7. Reconciliation TODO against brief.md

Decisions made in this plan that supersede or refine brief.md v0.3:

- **§6 QC_FILTER → split into NANOPLOT_RAW + CHOPPER + FILTLONG + NANOPLOT_CLEAN.** Brief lists `chopper, NanoPlot`; this plan splits NanoPlot into pre/post passes and adds Filtlong for identity-weighted top-quality selection after chopper's threshold cut. Dorado handles adapter trimming at basecalling, so no separate trim stage is needed.
- **§6 add BANDAGE_NG, BLAST_VALIDATE, ORGANELLE_MAP.** Brief does not list visualisation or BLAST validation stages (ANNOTATE is now in brief.md §6).

## 8. Remaining open questions

Structural / architectural questions that need answering before or during
early development, distinct from the parameter tuning in §9. Bigger-picture
"which shape is the pipeline?" items land here; "what value does knob X
take?" items land in §9.

1. **Read-level fallback arm** for samples too degraded to assemble ([brief.md §8.2](brief.md)). Defer to post-P5 unless P5 testing surfaces a clear need. More likely to bite at skim depth.
2. **Assembly strategy — bait-then-refine?** Do we try a quick-and-dirty draft assembly with read baiting, and then re-map clean reads back to that to pick up any reads that baiting didn't capture? Then assemble those reads properly? Structural question (adds a whole pass), not a knob — parking here rather than in §9.
3. **Plant organelle annotator (stage 12 `ANNOTATE`, plant branch).** Tool choice deferred. Candidates: **[Chloe](https://github.com/ian-small/chloe)** (Julia, plastid-focused), **PGA**, **Mitofinder** (mt-only), or an in-house HMM-based annotator built on Prokka + custom HMM libraries. Note: any single tool may not cover both cp and mt; a split choice is on the table. Decision needed before P4. See [this paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC9503105/pdf/ijms-23-10804.pdf) for deep dive.
4. **Organelle map renderer (stage 14 `ORGANELLE_MAP`).** Tool choice deferred. Candidates: **[pyCirclize](https://github.com/moshi4/pyCirclize)** (Python, actively maintained, GenBank input), **[plotgenes](https://github.com/aaronphillips7493/plotgenes)** (from CLAW's author), Circos with a custom GenBank→config converter, or an in-house Jinja+SVG renderer that reads GFF/GenBank directly (would fit the report renderer container per [§6a.5](06a-reports.md#6a5-report-implementation-boilerplate)). Inline SVG output preferred so it drops into the report. Decision needed before P4.

## 9. Fine tuning (post-prototype benchmarking)

Items to benchmark **after** the first end-to-end prototype is running on
real data (roughly end of P3). Each entry is a knob or design choice where
we've picked a sensible default from prior art or first principles, but the
right value / choice will only become obvious with our own test data in
hand.

**Flye tuning reference.** For items 3–6 and 9 below (all Flye knobs), the
authoritative source for parameter behaviour, valid combinations, and
mode-specific caveats is the vendored
[Flye user guide](../reference-material/flye-user-guide.md). Consult it
before designing any Flye benchmark — it documents the `--meta` /
`--asm-coverage` incompatibility, read-mode implications for polishing
iterations, and `--genome-size` estimator behaviour that these experiments
depend on.

**Format:** `knob` → `current default` → `experiment` → `success metric`.

| # | Knob / choice | Current default | Experiment | Success metric | Cross-ref |
|---|---|---|---|---|---|
| 1 | **RECRUIT filter strictness** (plant cp) | `samtools view -F 4 -q 1` (any primary mapping) | Compare against ptGAUL-style PAF filter: `alignment_length/read_length ≥ 0.7` and `alignment_length ≥ 1 kb`. | Cleaner input to METAFLYE (fewer off-target reads) without a recall hit that pushes COVERAGE_GATE below MIN on borderline samples. Measure: recruited-read count, estimated coverage, final assembly N50 + edge count. | §2 stage 5, §3.6 |
| 2 | **COVERAGE_GATE MIN/MAX** (per kingdom / organelle) | plant cp 30×/500×; animal mt 30×/300× | Sweep MIN downward on marginal plant / animal samples to find the floor where assembly still recovers a canonical structure. Sweep MAX to find the point where metaFlye graph tangles reappear. | Widest MIN–MAX band that still produces canonical assemblies. Failed samples at MIN vs. successful ones just above tells us the real floor. | §2.1.1 |
| 3 | **Flye `--genome-size` hint** | plant `2m`, animal `20k` | Vary the plant hint (500 kb, 1 m, 2 m, 4 m) on a plant sample where mt is present. Does under-hinting fragment mt? Does over-hinting hurt cp assembly? | Assembly contiguity for whichever organelle is under scrutiny (N50, edge count). | §3.4 |
| 4 | **Flye `--asm-coverage`** (plant) | `MAX − 20` (= 480 with default MAX 500) | Vary in 50× steps around 480; also try without the flag. | Plastid edge-count + total assembly length converging on canonical 3-edge quadripartite. | §3.5 |
| 5 | **Flye read-mode** (`--nano-hq` vs `--nano-corr` after upstream correction) | `--nano-hq` on raw Dorado SUP R10.4.1 | Insert an error-correction stage (HERRO / Canu-correct / Ratatosk) between FILTLONG and RECRUIT; run Flye with `--nano-corr`. Compare against the default. | Assembly indel rate + panel ORF pass rate. If gain is modest (<5% ORF-pass improvement) the runtime cost isn't justified. | §2 stage 7, §8 (was open Q4) |
| 6 | **Polish depth** (Flye built-in ± medaka) | Flye `--iterations 1` (Flye's default under `--nano-hq`); MEDAKA stage deferred. | Baseline: Flye-default polish only. Benchmark on samples that fail ORF validation: (a) Flye default (current), (b) Flye `--iterations 2–3`, (c) add a single medaka pass on top of Flye's polish, (d) medaka iterative (2–3 passes). Only un-defer the MEDAKA stage if (c) or (d) show meaningful ORF gains over (a)/(b). | Panel ORF pass rate vs. wall-clock cost. Look for diminishing returns; Flye's built-in polish is free so it's the baseline. | §2 stage 7 + deferred stage 8, brief.md §8.3 |
| 7 | **CHOPPER thresholds** | min length + min mean Q from config (values TBD in P1) | Sweep min length (500 / 1000 / 2000 bp) and min Q (7 / 10 / 12) on plant and animal skim data. | Recruited-read yield after RECRUIT and estimated coverage post-gate. Aim for the loosest filter that doesn't push borderline samples below COVERAGE_GATE MIN. | §2 stage 2 |
| 8 | **FILTLONG `--keep_percent`** | 95 | Compare 90 / 95 / 99 / disabled. FILTLONG removes reads by identity-weighted quality percentile — how much low-Q tail does metaFlye actually tolerate at skim depth? | Assembly contiguity + panel ORF pass rate vs. recruited-read count. | §2 stage 3 |
| 9 | **Assembler choice** ([brief.md §8.1](brief.md)) | metaFlye | Benchmark against GetOrganelle, Oatk, and (for plants) ptGAUL itself as an end-to-end alternative. | Assembly completeness (canonical structure, panel gene recovery), indel rate after polish, within-kingdom mixture behaviour, runtime, packaging complexity. | §2 stage 7, brief.md §8.1 |
| 10 | **BIN_TARGET coverage-spike threshold** | `2× the median non-target coverage` (working assumption) | Once real data exists, plot per-contig coverage distributions and pick a threshold empirically. Currently a first-principles guess. | Correct target-contig selection rate on the P3 clean test datasets (all four kingdom-organelle combinations). | §2 stage 10 |
| 11 | **Miniprot per-locus representative count** | 5–10 proteins per locus per kingdom | Vary between 1 (single canonical protein), 5–10 (current), and 50+ (broad clade coverage). Measure sensitivity on divergent test taxa. | Panel locus recovery rate on the most divergent test sample per kingdom. Runtime scaling for the top end. | §4.3 |

Each entry becomes a P-something task once the prototype is stable enough
to run these experiments repeatably. Grouping suggestion:

- **P-tune-A (read prep + recruit):** items 1, 7, 8. All touch stages 2–5 and can share test data.
- **P-tune-B (assembly parameters):** items 2, 3, 4, 5, 6, 9. Shared assembly benchmark harness.
- **P-tune-C (downstream):** items 10, 11. Each is stage-local.

Priority ranking will fall out of what breaks in P3–P5 testing; nothing
here blocks P0–P4.
