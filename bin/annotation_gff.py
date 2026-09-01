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
