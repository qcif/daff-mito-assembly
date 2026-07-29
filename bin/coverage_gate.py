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
