# CONSTITUTION

**Project:** wf5 — Organelle Genome Assembly + Barcode Recovery for Biosecurity Identification
**Status:** v1 · **Companion:** [brief.md](brief.md), [spec/](spec/)

This file distils the non-negotiable principles that govern every design
decision in this workflow. It is deliberately short. When the brief, spec,
or task documents disagree, the constitution wins; the other documents
are updated to reconcile.

Principles fall into three tiers:

- **Hard constraints** — client requirements. Not open to revision.
- **Design principles** — the shape of the solution. Revised only with
  explicit reconciliation of every downstream doc.
- **Engineering rules** — how we build and ship. Revised with a pull
  request + CI evidence that the change is safe.

---

## Hard constraints (client)

### 1. Scope is organelle assembly + barcode extraction, nothing else.

The pipeline assembles organelle genomes from ONT reads and extracts
taxonomic barcode loci. It stops at barcode extraction; taxonomic
assignment, database search, and phylogenetics belong to Taxodactyl.
Fungi, bacteria, nuclear rDNA/ITS, and cross-kingdom recovery from a
single submission are out of scope
([brief.md §2](brief.md)).

### 2. Kingdom is declared at intake and is a hard pre-assembly gate.

Every sample row declares an `assembly_target` from a fixed enum
(`animal_mt`, `plant_pt`, `plant_mt`). Rows without a valid target are
rejected at parse. This declaration selects the reference index,
coverage limits, genetic-code tables, and BLAST DB for the entire run
of that sample ([brief.md §3.1](brief.md),
[spec §1a](spec/01-pipeline-flow.md#1a-engineering-constraints)).

### 3. One sample row → one organelle assembly.

No within-sample fan-out. A plant sample that wants both plastid and
mitogenome assembled submits two rows sharing the same reads. This
keeps reference selection, size hints, coverage limits, and genetic
codes unambiguous per-run and per-stage
([spec §1a](spec/01-pipeline-flow.md#1a-engineering-constraints)).

### 4. Non-recruited reads are discarded.

Off-target material is not retained. If a misdeclared kingdom is
suspected, resubmit with corrected metadata
([brief.md §3.3](brief.md)).

---

## Design principles

### 5. Host exclusion is positive recruitment, not negative depletion.

The kingdom gate recruits reads *toward* the target's organelle
reference panel; it does not deplete reads *against* a host reference.
Positive selection needs no host reference (the submitter may not know
the precise host) and cannot discard a genuine target read for
incidentally resembling the host. Recruitment is coarse enrichment;
precise target/background separation happens post-assembly at per-contig
binning ([brief.md §3.2](brief.md)).

### 6. Assemble first, then extract.

Two logically distinct stages: assemble the organelle, then extract loci
from the assembly using a config-driven panel. Keeping these decoupled
lets the barcode panel evolve with Taxodactyl without re-running
assembly, and lets the assembled + annotated organelle be a useful
standalone output ([brief.md §3.4, §3.6](brief.md)).

### 7. Negative clarity — three failures, three signals.

A degraded sample must never produce output indistinguishable from a
confident negative. Three distinct negative states are individually
reported:

- **low_coverage** — coverage gate soft-fail
  ([spec §2.1](spec/02-stages.md#21-coverage-gate-2-stage-6))
- **no_assembly** — recruited but nothing assembled
- **no_barcode** — assembled but no locus extractable

Each surfaces in both the per-sample `report.html` and the run-level
`run-report.html` ([brief.md §3.8](brief.md),
[spec §6a](spec/06a-reports.md)).

### 8. Cross-sample failure isolation.

A single failed sample must not block any sibling sample in the same
invocation. Coverage-gate outcomes are *data* (a status marker file),
not a Nextflow error. Every per-sample process is wired so its failure
mode is a marker + skip-to-COLLATE, not a run abort
([spec §2.1.4](spec/02-stages.md#214-cross-sample-failure-isolation)).

### 9. Locus panel is versioned config, not hardcoded.

The extraction panel is parsed at runtime from a versioned config
(`assets/loci.json`, schema `wf5/loci-panel/v1`), grouped by origin
(`animal_mt`, `plant_pt`, `plant_mt`), values are canonical NCBI gene
symbols. The panel evolves with Taxodactyl's accepted-loci set; the two
must not drift ([brief.md §3.7](brief.md),
[spec §4.3](spec/04-reference-data.md#43-protein-panel-for-miniprot_extract)).

### 10. Reference bundle is versioned, immutable, and provenance-tracked.

All reference data (recruit indices, BLAST DBs, protein panel) is built
by a single scripted pipeline into a versioned directory
(`refs/v<YYYY.MM>/`) with a `manifest.json` recording SHA256, source
URL, and release tag for every input. The bundle version is emitted
into every run's `metadata.json` and `run_manifest.json` so every result
is traceable to the exact references used
([spec §4.4](spec/04-reference-data.md#44-consolidated-build-script)).
Reference data is stored in public Azure bucket for download by workflow
consumers such as CI runs and Azure Batch nodes.

---

## Engineering rules

### 11. Every process runs in a container. Every tag is pinned.

No stage may rely on host-installed tools or conda envs. Every
`modules/local/*.nf` declares a `container` directive (or resolves one
via `conf/containers.config`). Tags are pinned by SHA or explicit
version — **never `latest`**. CI fails if any process lacks a container
or uses `latest`. Applies equally to Python + Jinja report stages
([spec §1a](spec/01-pipeline-flow.md#1a-engineering-constraints),
[spec §5b](spec/05-test-data.md#5b-ci)).

### 12. Prefer biocontainers; bespoke images are lean and deps-only.

Reach for a published biocontainer or nf-core image first. When a
bespoke image is unavoidable, build the leanest possible base with the
interpreter + pinned deps only. **Source code is not baked in** — Python
lives under `bin/` and is auto-staged onto the container `PATH` by
Nextflow at runtime. `images.yml` rebuilds only when
`requirements.txt` or the `Dockerfile` changes; day-to-day code edits
ship with the next `nextflow run`
([spec §2.2](spec/02-stages.md#22-custom-logic-components)).

### 13. Files cross process boundaries as stageable `Path`s, not strings.

Bare `${params.<file>}` in a `script:` block is forbidden — it produces
a string the executor has nothing to stage, breaking silently on remote
executors. Materialise every reference into a value channel at the top
of `main.nf`, or use `${file(params.x)}` inline. CI's stageable-file
lint enforces this
([spec §1a](spec/01-pipeline-flow.md#1a-engineering-constraints)).

### 14. Custom logic lives in `bin/`, is unit-tested, runs in its own container.

If more than ~10 lines of `awk`/`bash` accumulate inside a `.nf`
`script:` block, promote it to a `bin/*.py` script with a `pytest`
module under `scripts/tests/`. Each custom-logic component (C1–C7)
targets **100% branch coverage**. Off-the-shelf tool wrappers rely on
`-stub-run` for wiring and nightly integration for tool behaviour
([spec §2.2](spec/02-stages.md#22-custom-logic-components),
[spec §5a](spec/05-test-data.md#5a-tests)). Where multiple biocontainer
packages are required in one Nextflow process, a mulled container image should
be built according to
([spec §1.1](spec/01-pipeline-flow.md#1a-engineering-constraints).

### 15. Three test surfaces, three cadences.

- **`pytest` on `bin/*.py`** — every push.
- **`-stub-run` end-to-end** — every push (~60 s: DSL syntax, channel
  wiring, container pulls, params).
- **`-profile integration`** — nightly + on-demand (real tools, real
  fixtures, biology assertions).

No fourth framework. If a regression class emerges that none of the
three catch, revisit before adding tooling
([spec §5a](spec/05-test-data.md#5a-tests)).

### 16. Provenance is a first-class output.

Every per-sample `metadata.json` records tool versions, reference-bundle
version, chosen genetic-code table, and any decisions made under
uncertainty (e.g. animal-mt clade trial winner). Every run's
`run_manifest.json` records the samplesheet snapshot, reference-bundle
version, and pipeline commit. This is the biosecurity audit trail —
losing it defeats the workflow's purpose
([brief.md §5](brief.md),
[spec §2.2 C7](spec/02-stages.md#22-custom-logic-components)).

### 17. Reproducibility across laptop → HPC → cloud.

The workflow must run identically on a developer laptop, an HPC login
node, and a cloud batch executor (AWS/Azure Batch). This
is why rules 11–13 exist. Any decision that would work locally but
break on a remote executor is not a valid decision.

### 18. When in doubt, favour the choice that preserves auditability.

This is a biosecurity workflow. A subtly wrong answer that ships as
"result" is worse than a loud failure. Prefer explicit failure over
silent fallback, prefer versioned config over hardcoded value, prefer
data marker over pipeline exception, prefer pinned image over floating
tag. Every design tie-breaker resolves in favour of the auditor
reading the report six months from now.

### 19. Maintenance

All design decisions should be made with the goal of reducing long-term
maintenance of the project. This is why we try to rely on Biocontainers for
dependencies, and should only rely on third party dependencies (e.g. Python
libraries) when they offer significant returns on code volume, complexity
and/or quality. Tests are especially challenging, as they should offer as much
coverage as possible, while avoiding fragility due to changes in reference
data. Drift in test data is a primary source of long-term maintenance burden.

---

## Amendment procedure

- **Hard constraints (1–4):** change requires written client sign-off.
  Update the brief first, then the constitution, then the spec.
- **Design principles (5–10):** change requires a pull request that
  updates the constitution, the affected spec sections, and any task
  documents built on the old principle. No orphaned references.
- **Engineering rules (11–18):** change requires a pull request with
  CI evidence (or a documented CI update) that the new rule is
  enforced and the old failure mode is closed.

The constitution is the shortest document in the project. Keep it that
way — everything longer than a paragraph belongs in the spec.
