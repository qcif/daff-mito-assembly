"""Two-axis annotation confidence classification (spec §6a.2,
task 43b §5.4).

`annotation.cds_scores[]` carries `pident` and `qcovhsp` per feature —
two independent axes, never combined into a single score (rule 1 below
is the invariant that actually closes the risk this module exists for:
no composite/mean/weighted value is computed anywhere in this module,
the template, or the result class). Each axis is tested against its
own threshold; the *pair* of pass/fail outcomes selects one of four
affirmatively-worded quadrant labels, so a low-identity/high-coverage
feature — a complete gene in an under-referenced clade — reads as a
good result rather than as the absence of a warning.

Thresholds are calibrated against `INT-ANIMAL-01`'s real spread (task
43b §5.4): COX3 (39.5/96) and ATP6 (39.3/100) must land in
COMPLETE_UNDER_REFERENCED; COX1 (68.75/13) must land in
TRUNCATED_WELL_REFERENCED.
"""

from typing import Optional

PIDENT_THRESHOLD = 50.0
QCOVHSP_THRESHOLD = 70.0

COMPLETE_WELL_REFERENCED = "complete_well_referenced"
TRUNCATED_WELL_REFERENCED = "truncated_well_referenced"
COMPLETE_UNDER_REFERENCED = "complete_under_referenced"
FRAGMENTARY = "fragmentary"
UNSCORED = "unscored"

LABELS = {
    COMPLETE_WELL_REFERENCED: "Complete gene, well-referenced clade",
    TRUNCATED_WELL_REFERENCED: "Truncated fragment, well-referenced clade",
    COMPLETE_UNDER_REFERENCED: "Complete gene, under-referenced clade",
    FRAGMENTARY: "Fragmentary, under-referenced clade",
    UNSCORED: "Unscored — no reference hit for this gene",
}

DESCRIPTIONS = {
    COMPLETE_WELL_REFERENCED: (
        "High identity and high coverage against the reference panel."
    ),
    TRUNCATED_WELL_REFERENCED: (
        "High identity but low coverage — a truncated fragment of a "
        "well-referenced gene."
    ),
    COMPLETE_UNDER_REFERENCED: (
        "Low identity but high coverage — a complete gene in an "
        "under-referenced clade. This is not a failure."
    ),
    FRAGMENTARY: (
        "Low identity and low coverage — treat this call as tentative."
    ),
    UNSCORED: (
        "No blastp hit against this gene's own reference panel — "
        "explicitly unscored, not zero confidence."
    ),
}


def classify(pident: Optional[float], qcovhsp: Optional[float]) -> str:
    """Classify one `cds_scores` entry's two independent confidence
    axes into an affirmatively-labelled quadrant. `pident` and
    `qcovhsp` are never averaged, multiplied, weighted or reduced to a
    single number — each is tested against its own threshold and the
    pair of outcomes selects a label, nothing else."""
    if pident is None or qcovhsp is None:
        return UNSCORED
    well_referenced = pident >= PIDENT_THRESHOLD
    complete = qcovhsp >= QCOVHSP_THRESHOLD
    if well_referenced and complete:
        return COMPLETE_WELL_REFERENCED
    if well_referenced and not complete:
        return TRUNCATED_WELL_REFERENCED
    if not well_referenced and complete:
        return COMPLETE_UNDER_REFERENCED
    return FRAGMENTARY


def label(quadrant: str) -> str:
    return LABELS.get(quadrant, quadrant)


def description(quadrant: str) -> str:
    return DESCRIPTIONS.get(quadrant, "")
