# Task 15 — P1 stage 6: `COVERAGE_GATE`

**Phase:** P1 (from [spec §6](../spec/06-phases.md)).
**Goal:** Replace the P0 stub for stage 6 with a real coverage-estimation
+ soft-fail / passthrough / subsample gate driven by
[`bin/coverage_gate.py`](../bin/coverage_gate.py) (custom-logic component
C2 in [spec §2.2](../spec/02-stages.md#22-custom-logic-components)).

Per [spec §2.1](../spec/02-stages.md#21-coverage-gate-2-stage-6) the gate
estimates coverage as `total_recruited_bases / nominal_organelle_size`,
where nominal size comes from the per-target limits table, and dispatches
one of three branches (soft-fail / passthrough / subsample). The process
**always exits 0** — coverage decisions are data on `sample_status.json`,
never Nextflow errors ([spec §2.1.4](../spec/02-stages.md#214-cross-sample-failure-isolation)).

**Prerequisite:** [task 8 — RECRUIT](completed/8_recruit.md) is real (it
is; commit `a11a657` onward), so `RECRUIT.out.reads` carries a real
gzipped FASTQ into this stage.

**Exit criteria:**

- `-profile integration` produces, for every sample:
    - `results/<sample_id>/coverage_gate/<sample_id>.gated.fastq.gz`
    - `results/<sample_id>/coverage_gate/sample_status.json`
    - `results/<sample_id>/coverage_gate/coverage.json`
  Contents reflect real base counts and a real gate decision — not the
  stub's hardcoded `estimated_cov: 150.0`.
- Progressive-assertion RECRUIT block in
  [`tests/integration/assertions.sh`](../tests/integration/assertions.sh)
  uncommented and green: `sample_status.status == "ok"` on the three
  Tier 2 fixtures.
- COVERAGE_GATE runs in a dedicated container that has Python 3.12 +
  `seqkit` + `seqtk` on `PATH` (see §3).
- `bin/coverage_gate.py` has unit tests covering the three decision
  branches (`< MIN`, `MIN ≤ cov ≤ MAX`, `> MAX`) plus the empty-input
  edge case — no Nextflow required.
- `errorStrategy 'ignore'` is retained but the script never triggers it
  under normal operation: `set -euo pipefail` at the top guarantees a
  seqkit/seqtk crash surfaces as a Nextflow "ignored" status rather
  than a silent empty output.
- Both `-profile stub -stub-run` (fast CI) and `-profile integration`
  (nightly) go green end-to-end.

**Not in scope:**

- Sweeping MIN/MAX limits ([open question 2](../spec/07-open-questions.md)) —
  ship the spec defaults; tuning is a separate benchmark task.
- Mapping-based coverage estimation. Spec §2.1.2 pins bases-÷-nominal-size;
  a mapping estimate would need an assembly, which we don't have yet.
- Downstream conditional wiring of METAFLYE / COLLATE on soft-fail — the
  branch already exists in [main.nf:87](../main.nf#L87) from the stub
  scaffold. This task confirms it works with real status JSON; no
  channel-topology change.
- COLLATE's minimal-bundle path for soft-failed samples — its own task.

**Cross-cutting rules (from [spec §1a](../spec/01-pipeline-flow.md#1a-engineering-constraints)):**

- Container pinned by SHA, never `latest`.
- No host tools; `seqkit`, `seqtk`, and Python 3.12 all come from the
  container.
- `bin/coverage_gate.py` is auto-staged onto every process `PATH` by
  Nextflow's `bin/` convention — invoke it by bare name.

---

## 1. Per-target limits — `nextflow.config`

Add the coverage-limits table alongside the existing `seqtk_seed`
default in [`nextflow.config`](../nextflow.config):

```groovy
params {
    // ... existing ...
    seqtk_seed              = 42

    // Coverage-gate limits — spec §2.1.1 defaults.
    // Overridable per-run with e.g. --coverage_limits.animal_mt.min_cov 20
    coverage_limits = [
        animal_mt: [ nominal_size:  17000, min_cov: 30, max_cov: 300 ],
        plant_pt:  [ nominal_size: 150000, min_cov: 30, max_cov: 500 ],
        plant_mt:  [ nominal_size: 400000, min_cov: 30, max_cov: 300 ],
    ]
}
```

The Python script receives these as CLI flags (see §2), not by reaching
into `params` — keeps the script testable in isolation and keeps
Nextflow the single owner of parameter resolution.

## 2. `bin/coverage_gate.py` (new — custom logic C2)

Stdlib-only Python. Invokes `seqkit stats -T` for base counting and
`seqtk sample -s $seed` for subsampling — both required to be on
`PATH` inside the container (§3).

```python
#!/usr/bin/env python3
"""Coverage gate — spec §2.1.

Reads a recruited FASTQ, estimates coverage against a target-specific
nominal organelle size, and emits one of three outcomes:

  * status="low_coverage"     when estimated_cov < min_cov
  * status="ok" (passthrough) when min_cov <= estimated_cov <= max_cov
  * status="ok" (subsampled)  when estimated_cov > max_cov

Always exits 0 — coverage decisions are data, not errors (spec §2.1.4).
Non-zero exits are reserved for unexpected tool failures (seqkit/seqtk
crash, malformed input) so Nextflow's errorStrategy='ignore' can log
them without failing the run.
"""

import argparse
import json
import subprocess
import sys
from pathlib import Path


def total_bases(fastq: Path) -> int:
    """Sum of read lengths via `seqkit stats -T` (tab-delimited)."""
    out = subprocess.run(
        ["seqkit", "stats", "-T", str(fastq)],
        capture_output=True, text=True, check=True,
    ).stdout.splitlines()
    # Header + one data row; column 5 is "sum_len".
    header, row = out[0].split("\t"), out[1].split("\t")
    return int(row[header.index("sum_len")])


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--reads", type=Path, required=True)
    p.add_argument("--sample-id", required=True)
    p.add_argument("--nominal-size", type=int, required=True)
    p.add_argument("--min-cov", type=int, required=True)
    p.add_argument("--max-cov", type=int, required=True)
    p.add_argument("--seed", type=int, required=True)
    p.add_argument("--out-fastq", type=Path, required=True)
    args = p.parse_args()

    bases = total_bases(args.reads)
    est = bases / args.nominal_size if args.nominal_size else 0.0

    status = {
        "sample_id": args.sample_id,
        "estimated_cov": round(est, 2),
        "min_required": args.min_cov,
        "max_allowed": args.max_cov,
        "total_recruited_bases": bases,
        "nominal_organelle_size": args.nominal_size,
    }
    coverage = {
        "pre_subsample_cov": round(est, 2),
        "post_subsample_cov": round(est, 2),
        "seed": args.seed,
        "subsampled": False,
        "fraction": 1.0,
    }

    if est < args.min_cov:
        status["status"] = "low_coverage"
        # Emit an empty gated fastq so downstream channel typing holds.
        args.out_fastq.write_bytes(b"")
    elif est > args.max_cov:
        frac = args.max_cov / est
        subprocess.run(
            ["bash", "-c",
             f"seqtk sample -s {args.seed} {args.reads} {frac} "
             f"| gzip > {args.out_fastq}"],
            check=True,
        )
        post = total_bases(args.out_fastq) / args.nominal_size
        status["status"] = "ok"
        coverage["post_subsample_cov"] = round(post, 2)
        coverage["subsampled"] = True
        coverage["fraction"] = round(frac, 4)
    else:
        # Passthrough — copy bytes, don't re-gzip.
        args.out_fastq.write_bytes(args.reads.read_bytes())
        status["status"] = "ok"

    Path("sample_status.json").write_text(json.dumps(status, indent=2))
    Path("coverage.json").write_text(json.dumps(coverage, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
```

Rationale for the shape:

- **Stdlib only.** Keeps the container thin (no `pip install`) and
  simplifies unit testing.
- **`seqkit stats -T` not raw `zcat | wc`.** Handles compressed +
  uncompressed input uniformly and is the documented spec choice.
- **Passthrough copies bytes.** Re-gzipping would waste CPU and could
  change output byte-identity for no gain; the recruit output is
  already a valid gzip.
- **Empty output on soft-fail.** Downstream branch in
  [main.nf:87–90](../main.nf#L87-L90) routes soft-fails past METAFLYE
  entirely, so the empty file is never opened. Keeping it present
  preserves the channel tuple shape.
- **Always exits 0 on coverage decisions.** Only `subprocess.run(...,
  check=True)` on seqkit/seqtk can raise — those are legitimate
  unexpected failures.

## 3. Container build

Use `mulled-build` (the same approach as RECRUIT in
[task 2](completed/2_containers.md)) to produce a conda-based mulled image
with Python 3.12 + `seqkit` + `seqtk`:

```bash
which mulled-build
# If not installed:
pipx install galaxy-tool-util
mulled-build build 'python=3.12,seqkit=2.13,seqtk=1.5'
```

`mulled-build` prints the generated `quay.io/biocontainers/mulled-v2-…` tag.
Retag with the tool-version convention used for RECRUIT and push:

```bash
docker tag \
    quay.io/biocontainers/mulled-v2-<hash>:<hash>-0 \
    neoformit/daff-wf5-coverage-gate:python3.12_seqkit2.13_seqtk1.5
docker push neoformit/daff-wf5-coverage-gate:python3.12_seqkit2.13_seqtk1.5
```

Pin by digest in [`conf/containers.config`](../conf/containers.config),
removing the `TODO P1` comment:

```groovy
withName: 'COVERAGE_GATE' {
    container = 'neoformit/daff-wf5-coverage-gate:python3.12_seqkit2.13_seqtk1.5@sha256:<pinned>'
}
```

## 4. `modules/local/coverage_gate.nf` (replace stub)

```groovy
// Stage 6 — coverage estimation + subsample/soft-fail gate.
// C2 custom logic — see spec §2.1 and §2.2.
// Always exits 0 on a coverage decision; errorStrategy 'ignore' guards
// only against unexpected seqkit/seqtk crashes.

process COVERAGE_GATE {
    tag            "${meta.sample_id}"
    label          'process_low'
    errorStrategy  'ignore'
    publishDir     "${params.outdir}/${meta.sample_id}/coverage_gate",
                   mode: 'copy', enabled: params.publish_intermediates

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.sample_id}.gated.fastq.gz"),
          path("sample_status.json"), path("coverage.json"), emit: gated

    stub:
    """
    touch ${meta.sample_id}.gated.fastq.gz
    echo '{"status": "ok", "estimated_cov": 150.0}' > sample_status.json
    echo '{"pre_subsample_cov": 150.0, "post_subsample_cov": 150.0}' > coverage.json
    """

    script:
    def limits = params.coverage_limits[meta.assembly_target]
    """
    coverage_gate.py \\
        --reads ${reads} \\
        --sample-id ${meta.sample_id} \\
        --nominal-size ${limits.nominal_size} \\
        --min-cov ${limits.min_cov} \\
        --max-cov ${limits.max_cov} \\
        --seed ${params.seqtk_seed} \\
        --out-fastq ${meta.sample_id}.gated.fastq.gz
    """
}
```

Notes:

- No `set -euo pipefail` — the Python script owns exit-code semantics
  and always returns 0 on a valid decision.
- `meta.assembly_target` is guaranteed to be one of
  `{animal_mt, plant_pt, plant_mt}` by
  [`bin/parse_samplesheet.py`](../bin/parse_samplesheet.py); no fallback
  needed.
- `process_low` label — the script does a couple of `seqkit stats`
  calls and at most one `seqtk sample` pass. CPU-cheap, memory-cheap.

## 5. Unit tests — `scripts/tests/test_coverage_gate.py`

Independent of Nextflow. Tests each decision branch by constructing a
tiny synthetic FASTQ, invoking `coverage_gate.py` as a subprocess with
fixed `--nominal-size` / `--min-cov` / `--max-cov`, and asserting on
`sample_status.json` + `coverage.json`.

Cases:

| # | Input bases | nominal | min | max | Expected `status` | Subsampled? |
|---|---|---|---|---|---|---|
| 1 | 500 kb   | 17 kb  | 30 | 300 | `ok`            | no  |
| 2 | 100 kb   | 17 kb  | 30 | 300 | `low_coverage`  | no  |
| 3 | 10 Mb    | 17 kb  | 30 | 300 | `ok`            | yes |
| 4 | 0        | 17 kb  | 30 | 300 | `low_coverage`  | no  |

Case 3 additionally asserts that `coverage.post_subsample_cov` is
within ±10% of `max_cov` — the fraction is float-derived so exact
equality isn't guaranteed.

Fixture construction: build the synthetic FASTQ inline with
`gzip.open(..., 'wb')` and 100-bp reads named `read_1`, `read_2`, …
— avoids checking in any test binaries.

## 6. Integration-test wiring

Uncomment the RECRUIT block in
[`tests/integration/assertions.sh:53-66`](../tests/integration/assertions.sh#L53-L66):

```bash
for sample in INT-ANIMAL-01 INT-PLANT-01-pt INT-PLANT-01-mt; do
    status_file="$OUTDIR/$sample/sample_status.json"
    if [[ ! -s "$status_file" ]]; then
        echo "FAIL: $sample/sample_status.json missing"
        FAILED=1
        continue
    fi
    status=$(jq -r .status "$status_file")
    if [[ "$status" != "ok" ]]; then
        echo "FAIL: $sample coverage gate status=$status (expected ok)"
        FAILED=1
    fi
done
```

**Prerequisite for `status == "ok"` on all three fixtures:** each
Tier 2 fixture must recruit enough bases to clear MIN. From
[task 8 §9](completed/8_recruit.md#9--findings--recruitment-yield-on-ci-fixtures)
the tiny CI fixtures fell well below 30×; the integration fixtures
were sized deliberately to clear the gate ([task 11 §6](completed/11_integration_tests.md)
generation recipe). Confirm on the first real run.

If a fixture falls short, options in order of preference:

1. Regenerate that fixture with a larger `seqtk sample` size (target
   ≥ 60× on nominal size).
2. Lower `min_cov` for that assembly target in
   [`conf/integration.config`](../conf/integration.config) via a
   `params.coverage_limits.<target>.min_cov` override — document why.

Publish path for `sample_status.json`: the module's `publishDir` writes
under `<outdir>/<sample_id>/coverage_gate/`. The assertion above
expects it at `<outdir>/<sample_id>/sample_status.json` for symmetry
with other per-sample artefacts (COLLATE later re-publishes it to the
sample root as part of the metadata bundle). Either update the
assertion to match the actual publish path, or wire the module to
publish `sample_status.json` at the sample root as well — pick one
before uncommenting.

## 7. `.github/workflows/tests.yml` (fast CI)

No change needed. The stub block emits the same hardcoded
`sample_status.json` shape, and `-stub-run` still touches the three
output files; the stub-e2e assertions in
[`tests/output/`](../tests/output/) continue to pass.

## 8. Verification

1. `pytest scripts/tests/test_coverage_gate.py` — all 4 cases pass.
2. `flake8 bin/coverage_gate.py scripts/tests/test_coverage_gate.py
   --max-line-length=100` — clean.
3. `nextflow run . -profile stub -stub-run` — 17/17 processes hit
   their stub, no regressions in `tests/output/STUB-01/`.
4. `nextflow run . -profile integration` locally (requires the
   fetched fixtures + refdata bundle from
   [task 14](completed/14_complete_integration_ci.md)):
   - Each sample's `sample_status.json` has `status: "ok"` and a
     plausible `estimated_cov > 30`.
   - `INT-ANIMAL-01` at ~2000 recruited reads × ~5 kb should sit
     around 500–600× on `animal_mt` (17 kb nominal) → subsampled to
     ~300×; verify `coverage.subsampled == true` and
     `post_subsample_cov` in `[270, 330]`.
   - Plant samples pass through without subsampling.
5. `bash tests/integration/assertions.sh` — the newly uncommented
   RECRUIT block reports OK for all three samples.
6. CI green on both `Tests` and `Integration` workflows on the PR.

## 9. Deliverables checklist

- [x] `bin/coverage_gate.py` — real implementation, stdlib-only,
      always-exit-0 on decisions.
- [x] `scripts/tests/test_coverage_gate.py` — 4 branch cases + edge case.
- [x] Container built via `mulled-build build 'python=3.12,seqkit=2.13,seqtk=1.5'`,
      retagged as `neoformit/daff-wf5-coverage-gate:python3.12_seqkit2.13_seqtk1.5`,
      and pushed to Docker Hub.
- [x] [`conf/containers.config`](../conf/containers.config) — SHA-pinned
      container for `COVERAGE_GATE`, `TODO P1` comment removed.
- [x] [`nextflow.config`](../nextflow.config) — `params.coverage_limits`
      table added.
- [x] [`modules/local/coverage_gate.nf`](../modules/local/coverage_gate.nf)
      — real `script:` block invoking `coverage_gate.py` (stub retained
      unchanged for `-profile stub`).
- [x] [`tests/integration/assertions.sh`](../tests/integration/assertions.sh)
      — RECRUIT/coverage-gate block uncommented and passing on the
      nightly run.
- [ ] Fast CI (`Tests`) + `Integration` workflows both green on the PR.

## 10. Notes / non-issues

- **Why no channel-topology change.** The `.branch { ok / failed }`
  split at [main.nf:87–90](../main.nf#L87-L90) uses
  `status_json.text.contains('"status": "ok"')`, which already handles
  both stub and real JSON. No downstream wiring change is required for
  this task.
- **`errorStrategy 'ignore'` still correct.** Reserved for
  seqkit/seqtk crashes (rare, out-of-scope failures). A soft-fail is
  *not* an ignored error — it's a valid zero-exit outcome that
  produces a real `sample_status.json`.
- **No mapping-based coverage.** Spec §2.1.2 pins bases-÷-nominal.
  Reconsider if benchmarking (see [open question 2](../spec/07-open-questions.md))
  shows the base-count proxy misfires on high-duplication samples.
- **`--asm-coverage` handoff to METAFLYE.** Plant plastid Flye takes a
  `--asm-coverage $((MAX-20))` per spec §3.5. That's METAFLYE's own
  task to consume, not COVERAGE_GATE's — no coupling here beyond the
  MAX value in `params.coverage_limits`.

## Outcomes

- Used `seqkit=2.13` (not `2.10.1` as originally drafted) — user confirmed
  this version at build time.
- `mulled-build` image:
  `neoformit/daff-wf5-coverage-gate:python3.12_seqkit2.13_seqtk1.5@sha256:efd1248269e148c4a380a9731cc26a339de1831509b2a6e183d989135d89cb24`
- Added a second unconditional `publishDir` for `sample_status.json` (pattern-scoped)
  so it always lands in `<outdir>/<sample_id>/coverage_gate/` regardless of
  `publish_intermediates`. This is required for the integration assertions and
  for COLLATE's eventual re-bundle.
- Integration results: INT-ANIMAL-01 at 325× was subsampled to 294× (post_subsample_cov
  within 2% of MAX=300); both plant samples passed through without subsampling.
- All 4 unit test cases pass inside the container; all 11 integration assertions pass.
