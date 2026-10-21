"""Shared parser for the merged annotation GFF — task 40.

Both halves of the confidence-scoring pipeline (``translate_annotation
_cds.py`` and ``annotation_scores.py``) need to enumerate exactly the
same set of CDS features from ``ANNOTATE``'s (C8) merged GFF, so they
agree on one definition of "one feature" regardless of provenance.
``bin/`` is staged onto the container PATH as a directory; both
scripts import this as a sibling module — the mechanism
``intervals.py`` already establishes (CONSTITUTION rule 19).

A rescued (MITOS2-provenance) feature's mRNA line carries an explicit
``Name=`` attribute (the canonical gene symbol — see
``annotate_summary.build_rescued_feature``). A miniprot winner's line
is copied verbatim from ``cds.gff`` (task 39 §3) and carries no
``Name=``, only ``Target=<protein_id> ...``; its gene symbol is
derived the same way ``annotate_summary.parse_cds_gff`` derives it.
"""


def parse_attrs(field: str) -> dict:
    attrs = {}
    for kv in field.strip().split(";"):
        if not kv or "=" not in kv:
            continue
        key, value = kv.split("=", 1)
        attrs[key] = value
    return attrs


def gene_symbol(attrs: dict) -> str:
    if "Name" in attrs:
        return attrs["Name"]
    target = attrs.get("Target", "").split()
    return target[0].rsplit("_", 1)[-1] if target else ""


def parse_annotation_gff(path) -> list:
    """Parse mRNA + CDS rows into one record per CDS feature.

    Malformed lines are skipped, not fatal (CONSTITUTION principle 8).
    """
    features = {}
    with open(path) as fh:
        for raw_line in fh:
            line = raw_line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) != 9:
                continue
            seqid, source, ftype, start, end, _score, strand, \
                _phase, attr = fields
            attrs = parse_attrs(attr)
            if ftype == "mRNA":
                fid = attrs.get("ID")
                if not fid:
                    continue
                try:
                    features[fid] = {
                        "id": fid,
                        "gene": gene_symbol(attrs),
                        "seqid": seqid,
                        "source": source,
                        "strand": strand,
                        "start": int(start),
                        "end": int(end),
                        "cds": [],
                    }
                except ValueError:
                    continue
            elif ftype == "CDS":
                parent = attrs.get("Parent")
                if parent in features:
                    try:
                        features[parent]["cds"].append(
                            (int(start), int(end)))
                    except ValueError:
                        continue
    return [f for f in features.values() if f["cds"]]


# ── Locus clustering (C5 + C8) ──────────────────────────────────────

def exons_overlap(exons_a: list, exons_b: list) -> bool:
    """True when any exon of A overlaps any exon of B. Exons are
    inclusive (start, end) GFF pairs."""
    return any(
        s1 <= e2 and s2 <= e1
        for s1, e1 in exons_a for s2, e2 in exons_b
    )


def _winner_key(rec: dict) -> tuple:
    """Highest identity wins; ties go to the leftmost, then the lowest
    ID, so the choice does not depend on input order."""
    return (
        -rec["identity"], rec["start"], rec["end"],
        rec.get("mrna_id", ""),
    )


def _components(recs: list) -> list:
    """Connected components of ``exons_overlap`` over same-strand
    records on one contig, found by sweeping the sorted exons: records
    whose exons fall in one merged run are joined."""
    parent = list(range(len(recs)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    exons = sorted(
        (s, e, i) for i, rec in enumerate(recs) for s, e in rec["cds"])
    run_end = None
    run_owner = None
    for start, end, i in exons:
        if run_end is not None and start <= run_end:
            parent[find(i)] = find(run_owner)
            run_end = max(run_end, end)
        else:
            run_end, run_owner = end, i

    components = {}
    for i, rec in enumerate(recs):
        components.setdefault(find(i), []).append(rec)
    return list(components.values())


def cluster_loci(records: list, gene_key) -> dict:
    """Collapse redundant panel hits into one winner per gene locus.

    Hits share a locus when they have the same gene (per ``gene_key``),
    contig and strand and an exon of one overlaps an exon of the other,
    transitively. A hit lying in another's intron is a separate locus,
    as is the same gene elsewhere (an inverted-repeat copy). Each
    record needs ``gene``, ``seqid``, ``strand``, ``cds`` (exon pairs)
    and ``identity``; ``start`` and ``end`` are set on it. Returns
    {gene: [winner, ...]}.
    """
    groups = {}
    for rec in records:
        rec["start"] = min(s for s, _e in rec["cds"])
        rec["end"] = max(e for _s, e in rec["cds"])
        key = (gene_key(rec), rec["seqid"], rec["strand"])
        groups.setdefault(key, []).append(rec)

    winners_by_gene = {}
    for (gene, _seqid, _strand), recs in groups.items():
        winners_by_gene.setdefault(gene, []).extend(
            min(component, key=_winner_key)
            for component in _components(recs)
        )
    for winners in winners_by_gene.values():
        winners.sort(key=lambda r: (r["seqid"], r["start"], r["end"]))
    return winners_by_gene
