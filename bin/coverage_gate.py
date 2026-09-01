#!/usr/bin/env python3
"""Coverage gate — spec §2.1.

Reads a recruited FASTQ, estimates coverage against a target-specific
nominal organelle size, and resolves one of four outcomes against two
floors — a terminal hard floor and an advisory warn floor (task 35):

  * status="fail"                 when estimated_cov < hard_min_cov
                                   (no assembly attempted; empty gated
                                   FASTQ)
  * status="low_coverage"         when hard_min_cov <= estimated_cov
                                   < warn_cov (assembly attempted, but
                                   the sample is flagged as degraded;
                                   full read passthrough)
  * status="ok" (passthrough)     when warn_cov <= estimated_cov
                                   <= max_cov
  * status="ok" (subsampled)      when estimated_cov > max_cov

RECRUIT's positive selection does not separate the two plant organelles
(CONSTITUTION principle 5), so on `plant_mt` the recruited pool is
routinely dominated by plastid carry-over. Estimating from every
recruited base then over-states mitochondrial depth and passes samples
that should soft-fail (spec §2.1.5, task 25). Where the declared target
has a sibling organelle panel, reads are therefore split by panel first
and only target-assigned bases feed the estimate.

Always exits 0 — coverage decisions are data, not errors (spec §2.1.4).
Non-zero exits are reserved for unexpected tool failures (seqkit/seqtk
crash, malformed input) so Nextflow's errorStrategy='ignore' can log
them without failing the run.
"""

import argparse
import json
import subprocess
import sys
import tempfile
from pathlib import Path

# Sibling module in bin/, staged as a directory onto the container PATH.
from intervals import merge_intervals

# Sibling organelle panels scored alongside the declared panel. Mirrors
# SIBLING_PANELS in bin_target.py (C3) — the two stages must agree on
# what counts as a sibling (spec §3.7.1).
SIBLING_PANELS = {
    'plant_mt': ['plant_pt'],
    'plant_pt': ['plant_mt'],
    'animal_mt': [],
}

# How estimated_cov was derived, recorded for the auditor (rule 18).
BASIS_TARGET_ASSIGNED = 'target_assigned'
BASIS_TOTAL_RECRUITED = 'total_recruited'

MINIMAP2_PRESET = 'map-ont'


def total_bases(fastq: Path) -> int:
    """Sum of read lengths via `seqkit stats -T` (tab-delimited)."""
    out = subprocess.run(
        ["seqkit", "stats", "-T", str(fastq)],
        capture_output=True, text=True, check=True,
    ).stdout.splitlines()
    # Header + one data row; column 5 is "sum_len".
    header, row = out[0].split("\t"), out[1].split("\t")
    return int(row[header.index("sum_len")])


def map_to_panel(
    reads: Path,
    mmi: Path,
    paf: Path,
    threads: int,
) -> dict:
    """
    Map reads against one panel; return {read_id: (read_len, aligned_bp)}.

    `aligned_bp` is the union of query intervals over every alignment
    block for that read — the same merged-interval metric C3 applies to
    contigs (spec §3.7.2), so the two stages classify on like for like.
    """
    with open(paf, 'w') as fh:
        subprocess.run(
            ["minimap2", "-x", MINIMAP2_PRESET, "-t", str(threads),
             str(mmi), str(reads)],
            stdout=fh, stderr=subprocess.DEVNULL, check=True,
        )

    lengths = {}
    intervals = {}
    with open(paf) as fh:
        for line in fh:
            fields = line.split('\t')
            read_id = fields[0]
            lengths[read_id] = int(fields[1])
            intervals.setdefault(read_id, []).append(
                (int(fields[2]), int(fields[3])))

    return {
        read_id: (
            lengths[read_id],
            sum(end - start for start, end in merge_intervals(blocks)),
        )
        for read_id, blocks in intervals.items()
    }


def split_by_panel(panel_hits: dict, assembly_target: str) -> dict:
    """
    Assign each read to the panel it aligns to best; total the bases.

    `panel_hits` maps panel name → the dict returned by map_to_panel.
    A read is target-assigned where its merged aligned fraction against
    the declared panel strictly exceeds every sibling's. Ties fall to
    the sibling: a read equally explained by both organelles is not
    evidence of target depth (rule 18).

    Returns {'target_assigned_bases', 'sibling_assigned_bases',
    'target_merged_aligned_bases', 'reads_target', 'reads_sibling'}.
    """
    siblings = [p for p in panel_hits if p != assembly_target]
    declared = panel_hits[assembly_target]

    every_read = set()
    for hits in panel_hits.values():
        every_read.update(hits)

    target_bases = sibling_bases = target_aligned = 0
    reads_target = reads_sibling = 0

    for read_id in every_read:
        read_len, declared_aligned = declared.get(read_id, (0, 0))
        if not read_len:
            read_len = next(
                panel_hits[p][read_id][0]
                for p in siblings if read_id in panel_hits[p]
            )
        best_sibling = max(
            (panel_hits[p].get(read_id, (0, 0))[1] for p in siblings),
            default=0,
        )
        if declared_aligned > best_sibling:
            target_bases += read_len
            target_aligned += declared_aligned
            reads_target += 1
        else:
            sibling_bases += read_len
            reads_sibling += 1

    return {
        'target_assigned_bases': target_bases,
        'sibling_assigned_bases': sibling_bases,
        'target_merged_aligned_bases': target_aligned,
        'reads_target_assigned': reads_target,
        'reads_sibling_assigned': reads_sibling,
    }


def write_passthrough(out_fastq: Path, reads: Path) -> None:
    """
    Copy the recruited reads through to the gated FASTQ unmodified.

    Shared by the `low_coverage` and `ok` passthrough branches so there
    is exactly one copy of this line — reusing the `fail` branch's
    empty-file write on either of them would silently route a sample
    to the assembler with no reads (spec §2, task 35).
    """
    out_fastq.write_bytes(reads.read_bytes())


def estimate_target_bases(
    reads: Path,
    ref_dir: Path,
    assembly_target: str,
    threads: int,
) -> tuple[int, list, dict]:
    """
    Total the recruited bases attributable to the declared organelle.

    Returns (target_bases, sibling_panels_scored, diagnostics).
    `target_bases` is None where the split could not be made — no
    sibling panel is defined for the target, an index is missing from
    the bundle, or minimap2 failed. The caller then falls back to the
    whole recruited pool rather than failing the sample
    (CONSTITUTION principle 8).
    """
    siblings = SIBLING_PANELS.get(assembly_target, [])
    if not siblings:
        return None, [], {}

    required = [assembly_target] + siblings
    missing = [p for p in required if not (ref_dir / f'{p}.mmi').exists()]
    if missing:
        print(
            f'WARNING: reference panels {missing} missing from {ref_dir}'
            ' — sibling-organelle discrimination is disabled; coverage'
            ' is estimated from all recruited bases and may over-state'
            f' {assembly_target} depth',
            file=sys.stderr,
        )
        return None, [], {}

    with tempfile.TemporaryDirectory() as tmp:
        try:
            panel_hits = {
                panel: map_to_panel(
                    reads, ref_dir / f'{panel}.mmi',
                    Path(tmp) / f'{panel}.paf', threads,
                )
                for panel in required
            }
        except subprocess.CalledProcessError as exc:
            print(
                f'WARNING: minimap2 failed ({exc}) — sibling-organelle'
                ' discrimination is disabled; coverage is estimated'
                ' from all recruited bases',
                file=sys.stderr,
            )
            return None, [], {}

    diagnostics = split_by_panel(panel_hits, assembly_target)
    return diagnostics['target_assigned_bases'], siblings, diagnostics


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--reads", type=Path, required=True)
    p.add_argument("--sample-id", required=True)
    p.add_argument("--assembly-target", required=True)
    p.add_argument("--ref-dir", type=Path, required=True,
                   help="Reference bundle directory holding ${panel}.mmi")
    p.add_argument("--nominal-size", type=int, required=True)
    p.add_argument("--hard-min-cov", type=int, required=True)
    p.add_argument("--warn-cov", type=int, required=True)
    p.add_argument("--max-cov", type=int, required=True)
    p.add_argument("--seed", type=int, required=True)
    p.add_argument("--threads", type=int, default=1)
    p.add_argument("--out-fastq", type=Path, required=True)
    args = p.parse_args()

    bases = total_bases(args.reads)
    target_bases, siblings_scored, split = estimate_target_bases(
        args.reads, args.ref_dir, args.assembly_target, args.threads)

    if target_bases is None:
        basis = BASIS_TOTAL_RECRUITED
        gated_bases = bases
    else:
        basis = BASIS_TARGET_ASSIGNED
        gated_bases = target_bases

    est = gated_bases / args.nominal_size if args.nominal_size else 0.0

    status = {
        "sample_id": args.sample_id,
        "estimated_cov": round(est, 2),
        "hard_min_required": args.hard_min_cov,
        "warn_threshold": args.warn_cov,
        "below_warn_floor": est < args.warn_cov,
        "max_allowed": args.max_cov,
        "total_recruited_bases": bases,
        "nominal_organelle_size": args.nominal_size,
        "coverage_basis": basis,
        "sibling_panels_scored": siblings_scored,
        **split,
    }
    if bases:
        status["sibling_organelle_fraction"] = round(
            split.get('sibling_assigned_bases', 0) / bases, 4)
    coverage = {
        "pre_subsample_cov": round(est, 2),
        "post_subsample_cov": round(est, 2),
        "seed": args.seed,
        "subsampled": False,
        "fraction": 1.0,
        "coverage_basis": basis,
    }

    if est < args.hard_min_cov:
        status["status"] = "fail"
        # Skip assembly entirely — emit an empty gated fastq so
        # downstream channel typing holds.
        args.out_fastq.write_bytes(b"")
    elif est < args.warn_cov:
        status["status"] = "low_coverage"
        # Degraded, but assembly is still attempted — full passthrough,
        # exactly as the "ok" in-band branch below.
        write_passthrough(args.out_fastq, args.reads)
    elif est > args.max_cov:
        frac = args.max_cov / est
        subprocess.run(
            ["bash", "-c",
             f"seqtk sample -s {args.seed} {args.reads} {frac} "
             f"| gzip > {args.out_fastq}"],
            check=True,
        )
        # Subsampling is uniform over the whole pool, so the realised
        # reduction applies to the target-assigned fraction too.
        realised = total_bases(args.out_fastq) / bases if bases else 0.0
        status["status"] = "ok"
        coverage["post_subsample_cov"] = round(est * realised, 2)
        coverage["subsampled"] = True
        coverage["fraction"] = round(frac, 4)
    else:
        status["status"] = "ok"
        write_passthrough(args.out_fastq, args.reads)

    Path("sample_status.json").write_text(json.dumps(status, indent=2))
    Path("coverage.json").write_text(json.dumps(coverage, indent=2))
    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
