# Task 59 — Collapse off-target contigs into one bar on the Assembly coverage chart

**Phase:** P4a (spec 06-phases.md) — report presentation only.
**Goal:** On the per-sample report's Assembly tab, show all off-target
contigs as **one** summary bar in the *Per-contig mean coverage* chart,
instead of one bar per contig. Target and secondary contigs keep their
own bars.

**Depends on:** nothing outstanding. Builds on task 52's bucket mapping
(`bin/contig_bucket.py`) and the chart layout from task 54.

**Touches:** `bin/report/report.py` (`_coverage_chart_data()`),
`scripts/report/templates/components/assembly.html` (the
`#contig-coverage-chart` script block), `scripts/tests/test_report.py`,
`spec/06a-reports.md` (Assembly tab row). No pipeline, Nextflow or
container changes.

---

## 1. Overview

The Assembly tab has a horizontal bar chart: one bar per assembled
contig, with bar length showing the mean read depth Flye reports for
that contig. Bars are coloured by what `BIN_TARGET` decided the contig
is: green for target (the organelle we were asked for), orange for
secondary (a plausible second target candidate), and grey for
off-target (everything else, including nuclear contigs, contaminants
and, on plant samples, the sibling organelle).

When an assembly is clean this works well. But metaFlye assembles
everything that was recruited, and when recruitment pulls in a lot of
non-organelle material the assembly can contain dozens of off-target
contigs. The chart's height grows by 30 px per contig, so it fills the
page with grey bars. That buries the one or two bars the reader needs,
and pushes the row (statistics table | chart | assembly graph) far
down the page.

Real examples from samples already run:

| Sample | Contigs | Target | Secondary | Off-target |
|---|---|---|---|---|
| R3-barcode20 (round 3) | 94 | 0 | 0 | 94 |
| AZ-INSECT-barcode13 / R3-barcode13 | 35 | 1 | 0 | 34 |
| INT-ANIMAL-01 (integration fixture) | 10 | 1 | 0 | 9 |
| R2-INSECT-barcode06 | 5 | 1 | 1 | 3 |

R3-barcode20 renders a chart about 2,800 px tall (94 × 30 px). Its
off-target contigs are nuclear sequence pulled in by microsatellite
recruitment (see the Recruitment item in tasks/todo.md). One grey bar
per contig tells the reader nothing that the contig statistics modal
does not already list.

What the reader actually needs to know about off-target contigs is:

1. **How many there were** — a sign of how dirty the recruited pool
   was.
2. **What depth they had compared with the target** — off-target
   contigs at or above the target's depth are a warning sign. For
   example, in INT-PLANT-01-mt the off-target bucket includes
   sibling-organelle (plastid) contigs at 168× against a mitochondrial
   target at 5–7×.

One summary bar can show both. The per-contig figures are still in
the *Contig statistics* modal, so no information is lost from the
report.

## 2. Constraints that shape the design

- **Auditability (CONSTITUTION rule 18).** A summary must not hide a
  warning sign. A single mean would hide one off-target contig at
  high depth among many at low depth. The summary bar therefore always
  shows the **maximum** as well as the central value (§3.2), so the
  most worrying contig is still visible.
- **Single source of truth (rule 19, task 52).** Off-target means
  whatever `contig_bucket.contig_bucket()` returns as `off-target`. Do
  not add a second classification rule here. The grey colour and the
  "Off-target" label come from `BUCKET_COLOURS` / `BUCKET_LABELS`.
- **Logic in Python, not JS (rule 14).** Do the aggregation in
  `_coverage_chart_data()` so it is unit-tested. The template only
  draws what it receives.
- **Spec 06a: values visible without hovering by default.** The Validation-tab
  charts are a documented exception to this. This chart is not, so
  the contig count must be in the bar's **axis label**, not just in
  the hover text.
- **Only the off-target bucket is collapsed.** Target and secondary
  contigs are the ones the reader cares about. They keep individual
  bars labelled by contig ID. This holds even on a plant_mt sample with
  four target contigs.
- **The Validation tab is out of scope.** Its *Estimated coverage* chart
  already shows only selected target contigs.

## 3. Design

### 3.1 When to collapse

Always collapse when there is at least one off-target contig with a
non-null coverage. Do not use a threshold such as "collapse only above N
contigs". A threshold makes the chart look different from sample to
sample for no biological reason, and it adds a magic number to tune.

- 0 off-target contigs → no off-target bar (same as today).
- 1 off-target contig → one bar, labelled as the bucket (§3.3). Not
  labelled by contig ID, so the chart is consistent across samples.
  Range whiskers are not drawn because min = max.
- ≥ 2 → one bar with range whiskers.

Contigs with `coverage: None` are still dropped before aggregation,
which matches the current `_coverage_chart_data()` behaviour and its
existing test.

### 3.2 What the bar shows

- **Bar length:** the **length-weighted mean coverage** of the
  off-target contigs, i.e. Σ(coverage × length) / Σ(length). This is
  the depth of the off-target pool taken as a whole, which is roughly
  "off-target bases assembled / off-target bp". It is used instead of a
  plain mean so that many tiny high-depth repeat fragments do not
  outweigh a few long contigs. If any contributing contig has no
  `length`, fall back to the unweighted mean and set a flag in the
  output (example: `weighted: false`). The hover text then says
  "mean", not "length-weighted mean". Do not fail silently.
- **Whiskers:** an asymmetric horizontal error bar running from the
  **minimum** to the **maximum** off-target coverage. The maximum is
  the auditability guarantee in §2.
- **Hover text:** number of contigs, total off-target bp, the
  length-weighted (or plain) mean, median, min and max. Each
  figure is formatted the same way the contig modal formats coverage
  (`N×`).

Example of the chart data the template would receive (illustrative
shape only, not a fixed schema):

```
bars: [
  {label: "contig_10", bucket: "target",  coverage: 259.0},
  {label: "Off-target (34 contigs)", bucket: "off-target",
   coverage: 7.4, min: 2.0, max: 19.0, median: 6.0,
   n: 34, total_bp: 412345, weighted: true},
]
```

Replacing the current parallel `x` / `y` / `colors` / `buckets` lists
with a list of bar records is acceptable and probably clearer. Choose
whichever reads better, and update the template and tests to match.
Drop the `colors` list if nothing reads it. The template already
colours by bucket via `bucket_legend`.

### 3.3 Labels and ordering

- The off-target bar's axis label is `Off-target (N contigs)`
  (`(1 contig)` when N = 1). Build it from
  `contig_bucket.BUCKET_LABELS['off-target']` so the wording cannot
  drift from the legend. Keep the format string as a module constant
  in `report.py`.
- Bar order follows `BUCKET_ORDER`: target, then secondary, then the
  off-target bar. Within target and secondary, keep the current contig
  order. Plotly horizontal bars draw the first category at the bottom.
  Check the rendered result and make sure the target bar appears **at
  the top** and the off-target summary **at the bottom**, reversing
  the category axis if necessary.
- Chart height: keep the `max(300, bars × 30)` formula but base it on
  the number of **bars**, not contigs. The chart then shrinks for
  samples with many off-target contigs.

### 3.4 Caption

Extend the caption under the chart (currently "Mean coverage per
contig, as estimated by the assembler (Flye).") with one sentence. It
should say that off-target contigs are summarised as a single bar
showing their length-weighted mean, with whiskers from the lowest to
the highest contig, and that per-contig figures are in *Contig
statistics*.

### 3.5 Things that must not change

- The *Contig statistics* modal still lists every contig individually.
- `target_bp` and the *Total assembly length (target)* figure are
  unchanged.
- The assembly graph colour key and node colouring are unchanged.
- `metadata.json` is unchanged. This is a view-layer change only, so
  COLLATE and the bundle schema are not touched.

## 4. Work items

1. `bin/report/report.py`: rewrite `_coverage_chart_data()` per §3.1–3.3,
   with helpers as needed (e.g. a small aggregator that takes the
   off-target contig list and returns the summary record). Constants
   for the label format.
2. `scripts/report/templates/components/assembly.html`: draw the new
   bar records. The off-target bar gets `error_x` (asymmetric,
   `array` = max − value, `arrayminus` = value − min) when n ≥ 2, plus
   a `hovertemplate` / `customdata` carrying the §3.2 figures. Size the
   chart by bar count. Update the caption (§3.4).
3. `spec/06a-reports.md`: in the Assembly tab row, after the per-contig
   coverage chart is mentioned, add a clause saying that off-target
   contigs are collapsed into a single summary bar (length-weighted
   mean, min–max whiskers, contig count in the label), with
   per-contig figures kept in the modal (task 59).
4. Render reports and check them visually (§5.3).

## 5. Tests

### 5.1 Unit — `scripts/tests/test_report.py` via `scripts/pytest.sh`

Replace or extend `test_coverage_chart_data_skips_none_coverage` and add
cases covering every branch of the new code (rule 14 targets 100 %
branch coverage on `bin/` logic):

- No off-target contigs → no summary bar. Target and secondary bars
  are unchanged and in order.
- Exactly one off-target contig → summary bar with n = 1, the singular
  label, min = max = value, and no whiskers flagged.
- Several off-target contigs with lengths → the length-weighted mean is
  correct for a hand-computed example (pick numbers where the weighted
  and plain means clearly differ). Also check min, max, median,
  total_bp and the plural label.
- An off-target contig missing `length` → falls back to the plain mean,
  with `weighted` false.
- Off-target contig with `coverage: None` → excluded from n and from
  every statistic. If it was the only one, there is no summary bar.
- One off-target contig at higher depth than the target → `max` on the
  summary equals that value (the auditability guarantee).
- Sibling-organelle and unknown classifications fall into the summary
  bar. This should follow from `contig_bucket()`. Assert it so a future
  bucket change is caught.
- Real fixture: `assembly_view(INT_ANIMAL_METADATA)` → one target bar
  (`contig_10`) plus one off-target summary with n = 9.

Add a render-level assertion through the existing `_render()` helper:
the rendered HTML contains the `Off-target (9 contigs)` label for
INT-ANIMAL-01, and contains no off-target contig IDs inside the chart
script block. Off-target IDs still appear in the modal, so scope the
assertion to the chart data rather than the whole page.

Run flake8 on `bin/report/report.py` and `scripts/tests/test_report.py`.

### 5.2 Integration

None needed. This changes presentation only, so no fixture output,
`metadata.json` field or `assertions.sh` check changes. The `-stub-run`
and integration RUN_REPORT steps will exercise the template as part of
the normal run.

### 5.3 Manual visual check (record in Outcomes)

Re-render the per-sample report from existing outputs, without
re-running the pipeline, for:

- tests/manual/output/round_3/R3-barcode20: 94 off-target contigs and
  no target. Expect a single bar and a short chart.
- tests/manual/output/round_3/R3-barcode13: 1 target + 34 off-target.
- tests/integration/output/INT-PLANT-01-mt: 4 target contigs at 5–7× plus
  2 sibling-organelle contigs, one at 168×. Expect the max whisker to
  reach about 168× and dwarf the target bars. That is the intended
  signal and should not be clipped.
- tests/integration/output/INT-PLANT-01-pt: 1 target + 1 secondary and
  no off-target. The chart should be unchanged apart from ordering.

Give the user the rendered report paths.

## 6. Non-goals

- Collapsing secondary or multiple target contigs.
- Changing the Validation tab's coverage chart, the contig modal or the
  assembly graph.
- A log-scale axis. The INT-PLANT-01-mt case may tempt you to add one.
  If the visual check suggests it is needed, note it in Outcomes as a
  follow-up rather than adding it here.
- Fixing the recruitment problem that produces many off-target contigs
  in the first place (tasks/todo.md, Recruitment).

## 7. Acceptance criteria

- On every sample with ≥ 1 off-target contig, the Assembly coverage
  chart shows exactly one grey off-target bar labelled with the contig
  count. The bar shows the length-weighted mean (or a flagged plain
  mean), and whiskers run from min to max when n ≥ 2.
- Target and secondary contigs keep individual, ID-labelled bars,
  ordered target → secondary → off-target from the top.
- The chart's height scales with bar count. R3-barcode20 renders a
  short chart.
- The *Contig statistics* modal still lists every contig.
- The unit tests in §5.1 pass, `report.py` has full branch coverage of
  the new code, and flake8 is clean.
- spec/06a-reports.md describes the collapsed bar.

---

## Outcomes

- `_coverage_chart_data()` now returns `{'bars': [...]}` (bar records,
  ordered target → secondary → off-target); the parallel
  `x`/`y`/`colors`/`buckets` lists and `colors` were dropped. Added
  `_off_target_summary()`, `OFF_TARGET_BUCKET` and
  `OFF_TARGET_LABEL_FORMAT` in `report.py`.
- Template: one Plotly trace per bucket, `error_x` whiskers when n ≥ 2,
  y-axis reversed so the target is at the top, height = bars × 30.
  Short value label drawn on bars; full figures in `hovertext`.
- Minor deviations: when the weighted mean falls back to the plain mean
  (missing length), `total_bp` is `None` and omitted from hover. Zero
  total length also falls back to the plain mean.
- Tests: 191 pass in `test_report.py`, `report.py` at 100 % branch
  coverage, flake8 clean.
- Visual check (renders in `tests/manual/output/task59_renders/`):
  - R3-barcode13: 1 target + "Off-target (34 contigs)" bar, as expected.
  - INT-PLANT-01-mt: 4 target bars (5–7×) + off-target bar with whisker
    reaching ~168×, not clipped.
  - INT-PLANT-01-pt: target + secondary bars only.
  - R3-barcode20: its report is in the terminal (soft-fail) state both
    before and after, so no Assembly chart is rendered; the 94-contig
    case is covered by unit tests only.
- On R3-barcode13 the off-target mean (~6×) is tiny next to the 259×
  target, so the value label overlaps the whisker. Follow-up requested
  and done: the x-axis switches to log scale (axis title says so) when
  min/max of all plotted values (bars and whiskers) is < 0.05
  (`LOG_SCALE_RATIO`, `_use_log_scale()`); also when a value is 0 (ratio 0); not when all values are 0.
  R3-barcode13 and INT-PLANT-01-mt go log; INT-PLANT-01-pt stays linear.
