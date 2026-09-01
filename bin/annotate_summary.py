#!/usr/bin/env python3
"""ANNOTATE (C8) — spec §2 stage 13a, §2.2 C8, task 31.

Merges ``MINIPROT_CDS``'s ``cds.gff`` with a specialist annotator's
non-CDS (tRNA/rRNA) features into one GFF3, keyed by contig name, and
writes ``annotation_summary.json``. Today the only specialist
annotator is MITOS2, run on ``animal_mt`` by the calling Nextflow
process; the plant branches have none yet (task 32, task 33) and stay
CDS-only.

**CDS features are copied verbatim from ``cds.gff``** — coordinates
are never adjusted, re-derived or re-sorted. The comprehensive protein
panel (task 29) carries several representative sequences per gene, so
more than one ``cds.gff`` hit commonly maps to the same genomic locus;
overlapping hits for the same gene are clustered and only the
highest-identity hit per cluster is emitted (mirroring
``validate_barcodes.py``'s locus clustering). Non-overlapping clusters
for the same gene are genuine distinct loci (inverted-repeat
duplication, plant-mt repeats) and are never collapsed.

**MITOS2's own protein-coding calls are compared against the
``cds.gff`` winners as an independent cross-check (``cds_crosscheck``)
and normally discarded** — emitting both would let the report show a
MITOS2 ``cox1`` at different coordinates from the miniprot ``cox1``
the barcode was cut from, exactly the incoherence the unified miniprot
pass (task 30) removed. **One bucket is the exception:** a call with
no miniprot counterpart at all (``annotator_only``) is rescued into a
real feature — attributed to its annotator, never conflated with a
miniprot CDS — because miniprot systematically misses short,
fast-evolving genes (``ATP8``) on under-referenced taxa (task 39).
``coordinate_conflicts`` stays untouched; miniprot remains
authoritative wherever both methods called a gene.

Stdlib-only: this module runs inside the MITOS2 biocontainer, which
carries no ``biopython``/``pandas`` (task 31 §2).

Always exits 0 (CONSTITUTION rule 8) — a poor or missing annotation is
supplementary information, not a pipeline error.
"""

import argparse
import json
import re
import sys
from pathlib import Path

GFF_COLUMNS = [
    "seqid", "source", "type", "start", "end",
    "score", "strand", "phase", "attributes",
]

STATUS_NO_ASSEMBLY = "no_assembly"
STATUS_OK_CDS_ONLY = "ok_cds_only"
STATUS_ANNOTATOR_FAILED = "annotator_failed"
STATUS_OK = "ok"
STATUS_NO_FEATURES = "no_features"

SCHEMA = "wf5/annotation-summary/v2"

MINIPROT_VERSION = "0.18"
MITOS_VERSION = "2.1.10"


def parse_attrs(field: str) -> dict:
    attrs = {}
    for kv in field.strip().split(";"):
        if not kv or "=" not in kv:
            continue
        key, value = kv.split("=", 1)
        attrs[key] = value
    return attrs


def parse_gff_line(line: str, line_no: int, path: Path):
    """Parse one tab-delimited GFF3 line into a dict, or None + a
    stderr warning if malformed (spec principle 8 — never abort)."""
    fields = line.split("\t")
    if len(fields) != 9:
        print(
            f"WARNING: skipping malformed GFF line {line_no} in "
            f"{path}: expected 9 columns, got {len(fields)}",
            file=sys.stderr,
        )
        return None
    try:
        seqid, source, ftype, start, end, score, strand, phase, attr = \
            fields
        record = dict(zip(GFF_COLUMNS, fields))
        record["start"] = int(start)
        record["end"] = int(end)
        record["attrs"] = parse_attrs(attr)
        record["raw"] = line
        return record
    except (ValueError, IndexError) as exc:
        print(
            f"WARNING: skipping malformed GFF line {line_no} in "
            f"{path}: {exc}",
            file=sys.stderr,
        )
        return None


# ── cds.gff (miniprot) parsing ──────────────────────────────────────

def parse_cds_gff(path: Path) -> list:
    """Parse mRNA + CDS rows into one record per miniprot hit.

    Mirrors ``validate_barcodes.parse_gff`` but keeps the raw mRNA
    line too (needed to reproduce it verbatim in the merged GFF) and
    does not need CIGAR/PAF context, since C8 does no translation.
    """
    mrna = {}
    with open(path) as fh:
        for line_no, raw_line in enumerate(fh, 1):
            line = raw_line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            record = parse_gff_line(line, line_no, path)
            if record is None:
                continue
            if record["type"] == "mRNA":
                attrs = record["attrs"]
                if "ID" not in attrs or "Target" not in attrs:
                    print(
                        f"WARNING: skipping malformed GFF line "
                        f"{line_no} in {path}: mRNA missing ID/Target",
                        file=sys.stderr,
                    )
                    continue
                target = attrs["Target"].split()
                protein_id = target[0]
                gene = protein_id.rsplit("_", 1)[-1]
                mrna[attrs["ID"]] = {
                    "seqid": record["seqid"],
                    "strand": record["strand"],
                    "gene": gene,
                    "identity": float(attrs.get("Identity", 0.0)),
                    "mrna_line": line,
                    "cds_lines": [],
                    "cds": [],
                }
            elif record["type"] == "CDS":
                parent = record["attrs"].get("Parent")
                if parent in mrna:
                    mrna[parent]["cds_lines"].append(line)
                    mrna[parent]["cds"].append(
                        (record["start"], record["end"]))
    return [rec for rec in mrna.values() if rec["cds"]]


def cluster_cds_by_gene(records: list) -> dict:
    """Group cds.gff records by gene symbol, then merge genomically
    overlapping hits (redundant panel representatives of the same
    locus), keeping the highest-identity hit per cluster.
    Non-overlapping hits are genuine distinct loci and stay separate
    (task 31 §4 point 6 — never de-duplicate by gene name)."""
    by_gene = {}
    for rec in records:
        rec["start"] = min(s for s, _e in rec["cds"])
        rec["end"] = max(e for _s, e in rec["cds"])
        by_gene.setdefault(rec["gene"], []).append(rec)

    winners_by_gene = {}
    for gene, recs in by_gene.items():
        recs = sorted(recs, key=lambda r: (r["seqid"], r["start"]))
        clusters = []
        for rec in recs:
            if (clusters
                    and clusters[-1][-1]["seqid"] == rec["seqid"]
                    and rec["start"]
                    <= max(r["end"] for r in clusters[-1])):
                clusters[-1].append(rec)
            else:
                clusters.append([rec])
        winners_by_gene[gene] = [
            max(cluster, key=lambda r: r["identity"])
            for cluster in clusters
        ]
    return winners_by_gene


# ── MITOS2 result.gff parsing ───────────────────────────────────────

def find_mitos_result_gffs(mitos_dir: Path) -> list:
    """MITOS2 writes ``result.gff`` directly under ``--outdir`` for a
    single-record input, or under one numbered subdirectory per input
    FASTA record (``<outdir>/<i>/result.gff``) for multi-record input.
    Either way, the seqid column already carries the real input record
    id (``mitos.gfffile.gffwriter`` writes the accession it is passed,
    which is the record's own id) — no index-based remapping needed,
    just gather every ``result.gff`` found."""
    return sorted(mitos_dir.rglob("result.gff"))


def parse_mitos_gff(mitos_dir: Path) -> tuple:
    """Returns (non_cds_records, cds_call_records, exons_by_parent).

    MITOS2 emits three rows per RNA gene (``ncRNA_gene``, ``tRNA`` or
    ``rRNA``, ``exon``) and two-or-more per protein gene (``gene``,
    one ``exon`` row per coding segment). ``tRNA``/``rRNA`` rows are
    the non-CDS features emitted into the merged GFF; ``gene`` rows
    are MITOS2's own CDS calls, used for the cross-check (§3/§4) and,
    for the ``annotator_only`` bucket, rescue (task 39). ``exon`` rows
    carry the reading frame the ``gene`` row lacks (phase is ``.`` on
    ``gene``) — collected here keyed by their raw ``Parent`` so rescue
    can assemble a gene's CDS block from its own exon(s).
    """
    non_cds = []
    cds_calls = []
    exons_by_parent = {}
    for gff_path in find_mitos_result_gffs(mitos_dir):
        with open(gff_path) as fh:
            for line_no, raw_line in enumerate(fh, 1):
                line = raw_line.rstrip("\n")
                if not line or line.startswith("#"):
                    continue
                record = parse_gff_line(line, line_no, gff_path)
                if record is None:
                    continue
                if record["type"] in ("tRNA", "rRNA"):
                    record["gene"] = record["attrs"].get("Name", "")
                    non_cds.append(record)
                elif record["type"] == "gene":
                    record["gene"] = record["attrs"].get("Name", "")
                    cds_calls.append(record)
                elif record["type"] == "exon":
                    parent = record["attrs"].get("Parent")
                    if parent:
                        exons_by_parent.setdefault(parent, []).append(record)
    return non_cds, cds_calls, exons_by_parent


# ── gene-name normalisation + completeness ──────────────────────────

ANTICODON_RE = re.compile(r"\(.*?\)")
FRAGMENT_SUFFIX_RE = re.compile(r"_\d+$")


def normalise_gene_symbol(raw: str) -> str:
    """Fold a gene symbol for cross-source comparison: uppercase,
    strip a MITOS2 anticodon parenthetical (``trnF(gaa)`` ->
    ``TRNF``), strip a MITOS2 fragment-index suffix (``atp8_0`` ->
    ``ATP8`` — task 39 §2.1; anchored to a trailing ``_<digits>`` so
    it never touches a real trailing digit like ``nad4l``'s, which
    has no underscore before it), strip an ``MT-`` locus prefix, drop
    spaces/dashes/underscores. The raw name in GFF attributes is left
    untouched — only comparison uses this folded form (task 31 §4
    point 4)."""
    s = ANTICODON_RE.sub("", raw)
    s = FRAGMENT_SUFFIX_RE.sub("", s).upper()
    if s.startswith("MT-"):
        s = s[3:]
    return re.sub(r"[\s_-]", "", s)


def build_alias_index(gene_set: dict) -> dict:
    """normalised alias -> canonical protein_coding gene name."""
    index = {}
    aliases = gene_set.get("protein_coding_aliases", {})
    for canonical in gene_set["protein_coding"]:
        index[normalise_gene_symbol(canonical)] = canonical
        for alias in aliases.get(canonical, []):
            index[normalise_gene_symbol(alias)] = canonical
    return index


def protein_coding_completeness(
    winners_by_gene: dict, gene_set: dict, rescued_genes=(),
):
    alias_index = build_alias_index(gene_set)
    found = set(rescued_genes)
    for gene in winners_by_gene:
        canonical = alias_index.get(normalise_gene_symbol(gene))
        if canonical is not None:
            found.add(canonical)
    canonical_set = set(gene_set["protein_coding"])
    missing = sorted(canonical_set - found)
    fraction = (
        len(found) / len(canonical_set) if canonical_set else None
    )
    return fraction, sorted(found), missing


# ── cross-check ──────────────────────────────────────────────────────

def overlaps(a: dict, b: dict) -> bool:
    return (
        a["seqid"] == b["seqid"]
        and a["start"] <= b["end"]
        and b["start"] <= a["end"]
    )


def canonicalise(raw: str, alias_index: dict) -> str:
    """Fold `raw` to its gene-sets canonical symbol where a mapping
    exists (so e.g. MITOS2's ``cob`` and miniprot's ``CYTB`` compare
    equal), else fall back to the plain normalised form."""
    normalised = normalise_gene_symbol(raw)
    return alias_index.get(normalised, normalised)


def cds_crosscheck(
    winners_by_gene: dict, mitos_cds_calls: list, alias_index: dict,
) -> tuple:
    """Returns (buckets, annotator_only_calls) — buckets is the dict
    written to ``annotation_summary.json`` minus ``annotator_only``
    (filled in later, once rescue has decided what stays held back);
    ``annotator_only_calls`` is the raw MITOS2 ``gene`` records with no
    miniprot counterpart under any name, passed on to
    ``rescue_annotator_only`` (task 39)."""
    miniprot_by_norm = {}
    for gene, winners in winners_by_gene.items():
        miniprot_by_norm.setdefault(
            canonicalise(gene, alias_index), []).extend(winners)

    annotator_by_norm = {}
    for call in mitos_cds_calls:
        annotator_by_norm.setdefault(
            canonicalise(call["gene"], alias_index), []).append(call)

    agreed, coordinate_conflicts = [], []
    annotator_only_calls = []
    for norm_name, calls in annotator_by_norm.items():
        miniprot_hits = miniprot_by_norm.get(norm_name)
        if not miniprot_hits:
            annotator_only_calls.extend(calls)
            continue
        if any(
            overlaps(call, hit)
            for call in calls for hit in miniprot_hits
        ):
            agreed.append(calls[0]["gene"])
        else:
            coordinate_conflicts.append(calls[0]["gene"])

    miniprot_only = [
        winners[0]["gene"]
        for norm_name, winners in miniprot_by_norm.items()
        if norm_name not in annotator_by_norm
    ]

    buckets = {
        "agreed": sorted(agreed),
        "miniprot_only": sorted(miniprot_only),
        "coordinate_conflicts": sorted(coordinate_conflicts),
    }
    return buckets, annotator_only_calls


# ── rescue (task 39) ─────────────────────────────────────────────────

REASON_OFF_PANEL = "off_panel"
REASON_OVERLAP = "overlap"
REASON_NO_EXON_DATA = "no_exon_data"

# §2.4 — a rescue candidate covering this much of the shorter of
# itself and a same-strand miniprot CDS is the same locus, not a
# neighbour. Deliberately coarse: measured overlaps separate into
# ~1-2% (adjacent genes, e.g. the canonical 7 bp ATP8/ATP6 junction)
# and 64-100% (genuine collisions), so anything in ~0.1-0.6 behaves
# identically on real data.
RESCUE_MAX_OVERLAP_FRACTION = 0.5


def overlap_fraction(a: dict, b: dict) -> float:
    """Overlap length as a fraction of the *shorter* of the two
    features, 0.0 when they don't overlap at all. The shorter feature
    is the denominator so the measure catches a small feature buried
    in a large one and a large feature swallowing a small one alike;
    a candidate-length denominator misses the second case."""
    if not overlaps(a, b):
        return 0.0
    shared = min(a["end"], b["end"]) - max(a["start"], b["start"]) + 1
    shortest = min(
        a["end"] - a["start"] + 1,
        b["end"] - b["start"] + 1,
    )
    return shared / shortest


def cluster_annotator_calls(calls: list, alias_index: dict) -> list:
    """MITOS2 splits an uncertain call into fragment-suffixed rows
    (``atp8_0``/``atp8_1``) that commonly overlap on the same locus —
    alternative candidate predictions, not genuine extra copies (task
    39 §2.1). Cluster overlapping same-strand calls **within one
    canonical gene name** and keep the longest-spanning one per
    cluster; non-overlapping clusters are genuine distinct loci and
    stay separate, mirroring ``cluster_cds_by_gene``'s treatment of
    redundant miniprot hits.

    Clustering is per-gene, never across genes: two *different* genes
    that happen to abut on the same strand are two genes, and a
    position-only cluster would drop one of them with no feature and
    no ``annotator_only`` record — the silent loss principle 18
    forbids.
    """
    by_gene = {}
    for call in calls:
        by_gene.setdefault(
            canonicalise(call["gene"], alias_index), []).append(call)

    kept = []
    for gene_calls in by_gene.values():
        gene_calls = sorted(
            gene_calls, key=lambda c: (c["seqid"], c["start"]))
        clusters = []
        for call in gene_calls:
            if (clusters
                    and clusters[-1][-1]["seqid"] == call["seqid"]
                    and clusters[-1][-1]["strand"] == call["strand"]
                    and call["start"]
                    <= max(c["end"] for c in clusters[-1])):
                clusters[-1].append(call)
            else:
                clusters.append([call])
        kept.extend(
            max(cluster, key=lambda c: c["end"] - c["start"])
            for cluster in clusters
        )
    return kept


def rescue_conflicts_with_miniprot(
    candidate: dict, winners_by_gene: dict,
) -> bool:
    """§2.4 — reject a rescue candidate that occupies the same locus
    as an existing miniprot CDS on the same strand, judged by
    ``RESCUE_MAX_OVERLAP_FRACTION`` rather than by any overlap at all
    (a few bp of junction overlap is normal between adjacent genes in
    a compact organelle genome, and is not evidence of duplication).
    Opposite-strand overlap is left alone entirely: dense organelle
    genomes legitimately carry genes on complementary strands at
    overlapping coordinates."""
    return any(
        winner["strand"] == candidate["strand"]
        and overlap_fraction(candidate, winner)
        >= RESCUE_MAX_OVERLAP_FRACTION
        for winners in winners_by_gene.values()
        for winner in winners
    )


def rescue_disqualification(
    call: dict, winners_by_gene: dict, alias_index: dict,
    protein_coding_set: set,
):
    """Why this call cannot be rescued, or ``None`` if it can.

    Off-panel takes precedence over overlap deliberately. A call
    outside the canonical set would never be rescued whatever its
    coordinates, so reporting ``overlap`` for it would point an
    auditor at a positional conflict that is not the real
    disqualifier — ``INT-ANIMAL-01``'s LAGLIDADG homing endonuclease
    happens to sit inside miniprot's ``COX1``, but the reason it is
    held back is that it is not a metazoan mitochondrial gene.
    """
    if canonicalise(call["gene"], alias_index) not in protein_coding_set:
        return REASON_OFF_PANEL
    if rescue_conflicts_with_miniprot(call, winners_by_gene):
        return REASON_OVERLAP
    return None


def transcript_parent_for(call: dict) -> str:
    """MITOS2 links a gene row's exon(s) by naming convention, not a
    shared feature row: ``ID=gene_atp6`` pairs with exon rows carrying
    ``Parent=transcript_atp6`` (task 39 §2.2) — never an explicit
    ``transcript`` feature."""
    gene_id = call["attrs"].get("ID", "")
    return "transcript_" + gene_id[len("gene_"):]


def build_rescued_feature(call: dict, exons: list, canonical_gene: str):
    """Assemble one rescued CDS feature from a MITOS2 gene call and
    its exon rows, in the shape ``write_gff`` already renders for
    miniprot winners — an mRNA line plus one CDS line per exon,
    preserving each exon's own phase (task 39 §2.2). Source is
    ``mitos``, distinguishing it from a miniprot feature in the GFF
    itself (task 39 §3)."""
    exons_sorted = sorted(exons, key=lambda e: e["start"])
    seqid = call["seqid"]
    strand = call["strand"]
    start = min(e["start"] for e in exons_sorted)
    end = max(e["end"] for e in exons_sorted)
    rescue_id = f"mitos_rescue_{canonical_gene}_{start}"
    mrna_line = "\t".join([
        seqid, "mitos", "mRNA", str(start), str(end), ".", strand, ".",
        f"ID={rescue_id};Name={canonical_gene};"
        f"OriginalName={call['gene']}",
    ])
    cds_lines = [
        "\t".join([
            seqid, "mitos", "CDS", str(exon["start"]), str(exon["end"]),
            ".", strand, exon["phase"],
            f"ID={rescue_id}_cds{i};Parent={rescue_id}",
        ])
        for i, exon in enumerate(exons_sorted, 1)
    ]
    return {
        "seqid": seqid, "strand": strand, "gene": canonical_gene,
        "start": start, "end": end,
        "mrna_line": mrna_line, "cds_lines": cds_lines,
    }


def rescue_annotator_only(
    annotator_only_calls: list, winners_by_gene: dict,
    alias_index: dict, protein_coding_set: set, exons_by_parent: dict,
) -> tuple:
    """Decide, per the four task 39 §2 guards, which ``annotator_only``
    calls become real CDS features. Returns (rescued, held_back) —
    ``rescued`` is a list of features shaped for ``write_gff``,
    ``held_back`` is ``[{"gene": ..., "reason": ...}]`` for the calls
    that stay in ``annotator_only``."""
    rescued, held_back = [], []

    # Both guards run on every *raw* call, before fragment clustering.
    # Clustering first would let a colliding fragment win its cluster
    # on length and discard a clean sibling before the guard ever saw
    # it, losing a gene that was rescuable (§2.4).
    survivors = []
    for call in annotator_only_calls:
        reason = rescue_disqualification(
            call, winners_by_gene, alias_index, protein_coding_set)
        if reason is None:
            survivors.append(call)
        else:
            held_back.append({"gene": call["gene"], "reason": reason})

    for call in cluster_annotator_calls(survivors, alias_index):
        canonical = canonicalise(call["gene"], alias_index)
        exons = exons_by_parent.get(transcript_parent_for(call), [])
        if not exons:
            # Boundary guard on external tool output: MITOS2's
            # documented shape always pairs a protein gene row with
            # at least one exon row, but nothing in this process
            # enforces that on MITOS2's side.
            held_back.append(
                {"gene": call["gene"], "reason": REASON_NO_EXON_DATA})
            continue
        rescued.append(build_rescued_feature(call, exons, canonical))
    return rescued, held_back


# ── output writers ───────────────────────────────────────────────────

def write_gff(
    out_path: Path, winners_by_gene: dict, non_cds: list, rescued=None,
):
    lines = ["##gff-version 3"]
    cds_blocks = [
        winner
        for winners in winners_by_gene.values()
        for winner in winners
    ] + list(rescued or [])
    cds_blocks.sort(key=lambda w: (w["seqid"], w["start"]))
    for winner in cds_blocks:
        lines.append(winner["mrna_line"])
        lines.extend(winner["cds_lines"])

    non_cds_sorted = sorted(
        non_cds, key=lambda r: (r["seqid"], r["start"]))
    for record in non_cds_sorted:
        lines.append(record["raw"])

    out_path.write_text("\n".join(lines) + "\n")


def run(
    cds_gff: Path, sample_id: str, assembly_target: str,
    gene_sets_path: Path, mitos_dir, genetic_code_annotate,
    genetic_code_cds, reference_data, annotator_exit, out_gff: Path,
    out_summary: Path,
) -> None:
    gene_sets = json.loads(gene_sets_path.read_text())["sets"]
    gene_set = gene_sets.get(assembly_target)

    cds_bytes = cds_gff.stat().st_size if cds_gff.exists() else 0

    if cds_bytes == 0:
        write_gff(out_gff, {}, [])
        summary = {
            "$schema": SCHEMA,
            "sample_id": sample_id,
            "assembly_target": assembly_target,
            "status": STATUS_NO_ASSEMBLY,
            "reason": "empty cds.gff — no target contig assembled upstream",
            "cds_source": {"miniprot": 0, "mitos": 0},
            "non_cds_source": None,
            "annotator_exit_code": None,
            "tool_versions": {"miniprot": MINIPROT_VERSION},
            "reference_data": None,
            "genetic_code_annotate": None,
            "genetic_code_cds": None,
            "genetic_code_agreement": None,
            "contigs_annotated": [],
            "feature_counts": {
                "gene": 0, "CDS": 0, "tRNA": 0, "rRNA": 0,
            },
            "protein_coding_completeness": None,
            "protein_coding_genes_missing": None,
            "cds_crosscheck": None,
            "cds_rescued": None,
            "canonical_gene_set": gene_set["name"] if gene_set else None,
        }
        out_summary.write_text(json.dumps(summary, indent=2) + "\n")
        return

    records = parse_cds_gff(cds_gff)
    winners_by_gene = cluster_cds_by_gene(records)
    n_cds = sum(len(w) for w in winners_by_gene.values())

    non_cds = []
    mitos_cds_calls = []
    exons_by_parent = {}
    non_cds_source = None
    found_result_gffs = []
    tool_versions = {"miniprot": MINIPROT_VERSION}
    if mitos_dir is not None:
        non_cds_source = "mitos2"
        tool_versions["mitos2"] = MITOS_VERSION
        non_cds, mitos_cds_calls, exons_by_parent = \
            parse_mitos_gff(mitos_dir)
        found_result_gffs = find_mitos_result_gffs(mitos_dir)

    n_trna = sum(1 for r in non_cds if r["type"] == "tRNA")
    n_rrna = sum(1 for r in non_cds if r["type"] == "rRNA")

    if n_cds == 0 and not non_cds:
        status = STATUS_NO_FEATURES
        reason = "no features parsed from cds.gff or the non-CDS annotator"
    elif non_cds_source is None:
        status = STATUS_OK_CDS_ONLY
        reason = (
            "no non-CDS annotator for this assembly_target — "
            "see task 32 (plant_pt) / task 33 (plant_mt)"
        )
    elif not non_cds:
        status = STATUS_ANNOTATOR_FAILED
        exit_text = (
            f"exited {annotator_exit}" if annotator_exit is not None
            else "exit code unknown"
        )
        outcome = (
            "wrote no result.gff" if not found_result_gffs
            else "wrote result.gff with no tRNA/rRNA features"
        )
        reason = (
            f"{non_cds_source} {exit_text} and {outcome} "
            f"under {mitos_dir.name}/"
        )
    else:
        status = STATUS_OK
        reason = None

    alias_index = build_alias_index(gene_set) if gene_set else {}

    crosscheck, rescued, cds_rescued = None, [], None
    if non_cds_source is not None:
        buckets, annotator_only_calls = cds_crosscheck(
            winners_by_gene, mitos_cds_calls, alias_index)
        protein_coding_set = (
            set(gene_set["protein_coding"]) if gene_set else set()
        )
        rescued, held_back = rescue_annotator_only(
            annotator_only_calls, winners_by_gene, alias_index,
            protein_coding_set, exons_by_parent)
        buckets["annotator_only"] = sorted(
            held_back, key=lambda h: h["gene"])
        crosscheck = buckets
        cds_rescued = sorted(r["gene"] for r in rescued)

    completeness, found, missing = (
        protein_coding_completeness(
            winners_by_gene, gene_set, rescued_genes=cds_rescued or ())
        if gene_set else (None, None, None)
    )

    contigs = sorted({
        w["seqid"] for winners in winners_by_gene.values()
        for w in winners
    } | {r["seqid"] for r in rescued} | {r["seqid"] for r in non_cds})

    genetic_code_agreement = (
        genetic_code_annotate == genetic_code_cds
        if genetic_code_annotate is not None
        and genetic_code_cds is not None
        else None
    )

    write_gff(out_gff, winners_by_gene, non_cds, rescued)

    n_cds_rescued = len(rescued)
    summary = {
        "$schema": SCHEMA,
        "sample_id": sample_id,
        "assembly_target": assembly_target,
        "status": status,
        "reason": reason,
        "cds_source": {"miniprot": n_cds, "mitos": n_cds_rescued},
        "non_cds_source": non_cds_source,
        "annotator_exit_code": annotator_exit,
        "tool_versions": tool_versions,
        "reference_data": reference_data,
        "genetic_code_annotate": genetic_code_annotate,
        "genetic_code_cds": genetic_code_cds,
        "genetic_code_agreement": genetic_code_agreement,
        "contigs_annotated": contigs,
        "feature_counts": {
            "gene": n_cds + n_cds_rescued + n_trna + n_rrna,
            "CDS": n_cds + n_cds_rescued,
            "tRNA": n_trna,
            "rRNA": n_rrna,
        },
        "protein_coding_completeness": completeness,
        "protein_coding_genes_missing": missing,
        "cds_crosscheck": crosscheck,
        "cds_rescued": cds_rescued,
        "canonical_gene_set": gene_set["name"] if gene_set else None,
    }
    out_summary.write_text(json.dumps(summary, indent=2) + "\n")


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cds-gff", type=Path, required=True)
    p.add_argument("--sample-id", required=True)
    p.add_argument("--assembly-target", required=True)
    p.add_argument("--gene-sets", type=Path, required=True)
    p.add_argument("--mitos-dir", type=Path, default=None)
    p.add_argument("--genetic-code-annotate", default="",
                   help="NCBI genetic-code table MITOS2 was run with")
    p.add_argument(
        "--genetic-code-json", type=Path, default=None,
        help="MINIPROT_CDS's genetic_code.json — carries the table "
             "actually selected and used for this sample's cds.gff, "
             "an independent value from genetic-code-annotate "
             "(task 38 §4)")
    p.add_argument("--reference-data", default="",
                   help="Non-CDS annotator reference-data tag, "
                        "e.g. refseq89m")
    p.add_argument(
        "--annotator-exit", type=int, default=None,
        help="Non-CDS annotator's process exit code, e.g. MITOS2's "
             "$MITOS_EXIT — omitted when no annotator was configured")
    p.add_argument("--out-gff", type=Path, required=True)
    p.add_argument("--out-summary", type=Path, required=True)
    args = p.parse_args()

    genetic_code_annotate = (
        int(args.genetic_code_annotate)
        if args.genetic_code_annotate else None
    )
    genetic_code_cds = (
        json.loads(args.genetic_code_json.read_text())["selected_table"]
        if args.genetic_code_json is not None else None
    )
    reference_data = args.reference_data or None

    run(
        args.cds_gff, args.sample_id, args.assembly_target,
        args.gene_sets, args.mitos_dir, genetic_code_annotate,
        genetic_code_cds, reference_data, args.annotator_exit,
        args.out_gff, args.out_summary,
    )
    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
