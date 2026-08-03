#!/usr/bin/env python3
"""Shared interval arithmetic — spec §3.7.2.

The merged-interval metric is the common currency of three consumers:

  * C2 (`coverage_gate.py`) merges per-read query intervals to decide
    which recruited reads are target-assigned (spec §2.1.5);
  * C3 (`bin_target.py`) merges per-contig query intervals to compute
    aligned fraction against each panel (spec §3.7.2);
  * the build-time panel masker (`scripts/refdata/mask_panel.py`)
    merges plastid-hit intervals per mitogenome before masking them.

All three must agree on what "merged aligned length" means, so the
routine lives here once rather than being copied per consumer
(CONSTITUTION rule 19). `bin/` is staged onto the container PATH as a
directory, so the runtime consumers import this as a sibling — the same
mechanism `bin_target.py` already uses for `plastid_canonicalise`. The
build-time consumer runs on the host and adds `bin/` to `sys.path`.
"""


def merge_intervals(intervals: list) -> list:
    """
    Merge overlapping and adjacent [start, end) intervals.

    Returns a sorted, non-overlapping list. Adjacent intervals
    (end == next start) are merged; empty input returns [].
    """
    if not intervals:
        return []
    ordered = sorted(intervals)
    merged = [list(ordered[0])]
    for start, end in ordered[1:]:
        if start <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return [(s, e) for s, e in merged]


def merged_length(intervals: list) -> int:
    """Total length covered by `intervals` once merged."""
    return sum(end - start for start, end in merge_intervals(intervals))
