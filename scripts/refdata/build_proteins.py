#!/usr/bin/env python3
"""
Build refs/<version>/proteins/<origin>/<gene>.faa.

Two modes:

  barcode (legacy) — for each (origin, gene) entry in assets/loci.json
    (schema wf5/loci-panel/v1), queries NCBI RefSeq protein via
    POST-based esearch → efetch (history/WebEnv, to avoid HTTP 414 on
    long ID batches) and writes one FASTA per locus. Superseded by
    `full` (task 29) but kept for reference / re-derivation.

  full (task 29, re-selected task 48) — builds the organelle's
    comprehensive protein-coding complement from CDS translations
    already present in the reference bundle's RefSeq organelle records
    (validate/refseq_pt.fa, validate/refseq_mt_{metazoa,viridiplantae}
    .fa), rather than issuing fresh per-gene NCBI queries. Reads a CDS
    jsonl file produced by scripts/refdata/parse_gbff_cds.py (run
    inside the neoformit/daff-wf5-scripts:test image, which has
    biopython) and assets/organelle_gene_sets.json for the canonical
    gene list + alias table, then selects representatives per gene by
    a deterministic, rank-aware breadth-first traversal of each
    candidate's GenBank taxonomy lineage — within one or more
    mutually-exclusive per-origin `--stratum` budgets, or a single
    `--max-per-locus`-sized stratum for origins given no `--stratum`
    (task 48 §3-§4).

Origin → NCBI taxid restriction (barcode mode only):
    animal_mt → txid33208 (Metazoa)
    plant_pt  → txid33090 (Viridiplantae)
    plant_mt  → txid33090 (Viridiplantae)

Requires: Python 3.10+, urllib (stdlib). No third-party dependencies —
the biopython-dependent GBFF parsing step lives in the separate
parse_gbff_cds.py, run inside a container (see that script's docstring).

Usage:
    python3 build_proteins.py barcode \
        --loci-config  assets/loci.json \
        --out          refs/v2026.07/proteins \
        [--ncbi-api-key KEY] \
        [--max-per-locus 50]

    python3 build_proteins.py full \
        --gene-sets       assets/organelle_gene_sets.json \
        --cds-jsonl       animal_mt=animal_mt_cds.jsonl \
        --cds-jsonl       plant_pt=plant_pt_cds.jsonl \
        --cds-jsonl       plant_mt=plant_mt_cds.jsonl \
        --out             refs/v2026.09_1/proteins \
        --refseq-release  236 \
        --stratum         animal_mt:Arthropoda=1000 \
        --stratum         animal_mt:Metazoa=100 \
        [--min-per-locus 5] [--max-per-locus 10]
"""

import argparse
import json
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path

EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# Origin key → NCBI [Organism] filter.
ORIGIN_TAXID = {
    "animal_mt": "txid33208[Organism]",   # Metazoa
    "plant_pt":  "txid33090[Organism]",   # Viridiplantae (plastid genes)
    "plant_mt":  "txid33090[Organism]",   # Viridiplantae (mitochondrial genes)
}

MIN_RECORDS = 3


def parse_args():
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest="mode", required=True)

    barcode = sub.add_parser("barcode", help="legacy per-locus NCBI fetch")
    barcode.add_argument("--loci-config", required=True)
    barcode.add_argument("--out", required=True)
    barcode.add_argument("--ncbi-api-key", default="")
    barcode.add_argument("--max-per-locus", type=int, default=50)

    full = sub.add_parser(
        "full", help="comprehensive panel from bundled RefSeq CDS records")
    full.add_argument("--gene-sets", required=True)
    full.add_argument(
        "--cds-jsonl", action="append", required=True,
        metavar="ORIGIN=PATH",
        help="repeatable: origin=path to a parse_gbff_cds.py jsonl output")
    full.add_argument("--out", required=True)
    full.add_argument("--refseq-release", required=True)
    full.add_argument("--min-per-locus", type=int, default=MIN_RECORDS)
    full.add_argument("--max-per-locus", type=int, default=10)
    full.add_argument(
        "--stratum", action="append", default=[],
        metavar="ORIGIN:TERM=BUDGET",
        help=(
            "repeatable: a rank-aware selection stratum for an origin, "
            "e.g. animal_mt:Arthropoda=1000. Strata for the same origin "
            "apply in declaration order and are mutually exclusive — a "
            "record matching an earlier stratum's lineage term is not "
            "eligible for a later one (task 48 §4.2). An origin with no "
            "--stratum falls back to a single unrestricted stratum sized "
            "by --max-per-locus."
        ),
    )

    return p.parse_args()


def load_loci(path: str) -> dict:
    with open(path) as fh:
        return json.load(fh)


def build_query(gene: str, origin_filter: str) -> str:
    """
    Build an NCBI protein query string for a canonical gene symbol +
    origin taxid restriction. NCBI [Gene Name] is case-insensitive.
    """
    return f'"{gene}"[Gene Name] AND {origin_filter} AND refseq[Filter]'


def _post_with_retry(endpoint: str, params: dict, timeout: int = 60,
                     retries: int = 4) -> bytes:
    """POST to an eutils endpoint with exponential backoff on 5xx errors."""
    body = urllib.parse.urlencode(params).encode()
    for attempt in range(retries):
        req = urllib.request.Request(
            f"{EUTILS_BASE}/{endpoint}",
            data=body,
            headers={"Content-Type": "application/x-www-form-urlencoded"},
        )
        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return resp.read()
        except urllib.error.HTTPError as exc:
            if exc.code >= 500 and attempt < retries - 1:
                wait = 2 ** attempt
                print(
                    f"  → {exc.code} on {endpoint}, retry in {wait}s "
                    f"({attempt + 1}/{retries})",
                    file=sys.stderr,
                )
                time.sleep(wait)
                continue
            raise
    raise RuntimeError(f"exhausted retries against {endpoint}")


def esearch(query: str, api_key: str = ""):
    """Return (webenv, querykey, count) for an esearch query."""
    params = {
        "db": "protein",
        "term": query,
        "retmax": "0",
        "usehistory": "y",
        "retmode": "json",
    }
    if api_key:
        params["api_key"] = api_key
    data = json.loads(_post_with_retry("esearch.fcgi", params))
    return (
        data["esearchresult"].get("webenv", ""),
        data["esearchresult"].get("querykey", ""),
        int(data["esearchresult"].get("count", 0)),
    )


def efetch_fasta(webenv: str, query_key: str, retmax: int,
                 api_key: str = "") -> str:
    """Fetch protein FASTA for a WebEnv history set."""
    params = {
        "db": "protein",
        "query_key": query_key,
        "WebEnv": webenv,
        "rettype": "fasta",
        "retmode": "text",
        "retmax": str(retmax),
    }
    if api_key:
        params["api_key"] = api_key
    return _post_with_retry("efetch.fcgi", params, timeout=120).decode()


def fetch_proteins(query: str, max_records: int,
                   api_key: str = "") -> tuple[str, int]:
    """
    Return (fasta_text, record_count) for an NCBI protein query.
    Uses esearch → efetch (WebEnv history) to avoid huge GET URLs.
    """
    delay = 0.12 if api_key else 0.40
    webenv, query_key, total = esearch(query, api_key=api_key)
    time.sleep(delay)
    if total == 0:
        return "", 0
    fetch_n = min(total, max_records)
    fasta = efetch_fasta(webenv, query_key, fetch_n, api_key=api_key)
    time.sleep(delay)
    record_count = fasta.count(">")
    return fasta, record_count


def stage_dir(out_root: Path):
    """Wipe/create a `.staging` sibling of out_root; return its Path."""
    import shutil
    staging = Path(str(out_root) + ".staging")
    if staging.exists():
        shutil.rmtree(staging)
    staging.mkdir(parents=True)
    return staging


def promote(staging: Path, out_root: Path):
    """Atomically replace out_root with the built staging directory."""
    import shutil
    if out_root.exists():
        shutil.rmtree(out_root)
    out_root.parent.mkdir(parents=True, exist_ok=True)
    staging.rename(out_root)


def main_barcode(args):
    loci_config = load_loci(args.loci_config)
    out_root = Path(args.out)
    staging = stage_dir(out_root)

    errors = []

    # Iterate origin keys in the panel; skip metadata keys like $schema.
    origins = [k for k in loci_config if not k.startswith("$")]
    for origin in origins:
        if origin not in ORIGIN_TAXID:
            print(
                f"[proteins] WARNING: unknown origin '{origin}' — skipping",
                file=sys.stderr,
            )
            continue
        genes = loci_config[origin]
        origin_dir = staging / origin
        origin_dir.mkdir(exist_ok=True)
        origin_filter = ORIGIN_TAXID[origin]

        for gene in genes:
            out_faa = origin_dir / f"{gene}.faa"
            print(
                f"[proteins] {origin}/{gene} ...",
                end=" ",
                flush=True,
            )
            try:
                query = build_query(gene, origin_filter)
                fasta, count = fetch_proteins(
                    query,
                    max_records=args.max_per_locus,
                    api_key=args.ncbi_api_key,
                )
            except Exception as exc:
                print(f"ERROR: {exc}", file=sys.stderr)
                errors.append(f"{origin}/{gene}: {exc}")
                continue

            if count < MIN_RECORDS:
                msg = (
                    f"{origin}/{gene}: only {count} records returned "
                    f"(min={MIN_RECORDS}); query: {query}"
                )
                print(f"FAIL ({count} records)", file=sys.stderr)
                errors.append(msg)
                continue

            out_faa.write_text(fasta)
            print(f"{count} records")

    if errors:
        print("\n[proteins] ERRORS — build incomplete:", file=sys.stderr)
        for e in errors:
            print(f"  {e}", file=sys.stderr)
        sys.exit(1)

    promote(staging, out_root)
    print(f"[proteins] done: {out_root}")


# ── full-complement mode (task 29) ──────────────────────────────────────

def load_gene_sets(path: str) -> dict:
    with open(path) as fh:
        return json.load(fh)["sets"]


def normalise_symbol(raw: str) -> str:
    """
    Fold a GenBank /gene qualifier for alias matching: uppercase, strip
    an 'MT-' locus prefix (vertebrate mitogenome convention, e.g.
    MT-CO1), and drop spaces/dashes/underscores.
    """
    s = raw.upper()
    if s.startswith("MT-"):
        s = s[3:]
    return re.sub(r"[\s_-]", "", s)


def build_alias_index(gene_set: dict) -> dict:
    """canonical GenBank-symbol-folded alias -> canonical panel gene name."""
    index = {}
    aliases = gene_set.get("protein_coding_aliases", {})
    for canonical in gene_set["protein_coding"]:
        index[normalise_symbol(canonical)] = canonical
        for alias in aliases.get(canonical, []):
            index[normalise_symbol(alias)] = canonical
    return index


def genus_of(organism: str) -> str:
    return organism.split()[0] if organism else "unknown"


# ── rank-aware breadth-first selection (task 48 §3) ─────────────────────
#
# Spends a stratum's budget as evenly as possible across the highest
# taxonomic ranks first, then recurses: before taking a second beetle,
# take a mite. Deterministic throughout — every ordering decision has
# an explicit sort key, with accession as the final tie-break, so two
# builds from the same inputs produce byte-identical output regardless
# of input order or dict/set iteration order (CONSTITUTION rule 18).

UNRESOLVED = object()  # sentinel bucket for ragged/short lineages


def _rel_lineage(rec: dict, root_term: str | None) -> list:
    """The record's lineage relative to a stratum root: everything
    after `root_term`, or the full lineage if `root_term` is None
    (the unrestricted fallback stratum) or not found."""
    lineage = rec.get("lineage") or []
    if root_term is not None and root_term in lineage:
        return lineage[lineage.index(root_term) + 1:]
    return list(lineage)


def _partition_by_rank(items: list, depth: int) -> list:
    """Group (rec, rel_lineage) items by the lineage term at `depth`.
    Records whose relative lineage is shorter than depth land in a
    single UNRESOLVED bucket, ordered last. Real groups are ordered
    largest-candidate-pool-first, ties broken by term name — both
    deterministic, independent of dict insertion order."""
    groups = defaultdict(list)
    unresolved = []
    for item in items:
        rec, rel = item
        if depth < len(rel):
            groups[rel[depth]].append(item)
        else:
            unresolved.append(item)

    ordered = sorted(groups.items(), key=lambda kv: (-len(kv[1]), kv[0]))
    if unresolved:
        ordered.append((UNRESOLVED, unresolved))
    return ordered


def _take_one_per_genus(items: list, budget: int, used_genera: set,
                        used_accessions: set) -> list:
    """Leaf selection: at most one record per genus, the longest
    translation wins within a genus, accession breaks ties, and the
    translation string itself is the final tie-break — real RefSeq
    data has duplicate (accession, gene) CDS records (e.g. isoform
    reannotations) with equal-length but non-identical translations,
    so accession alone does not always disambiguate. Genera are
    offered in largest-candidate-pool-first, name-second order, so the
    genera actually chosen when budget is scarce are deterministic."""
    by_genus = defaultdict(list)
    for rec, _rel in items:
        genus = genus_of(rec["organism"])
        if genus in used_genera or rec["accession"] in used_accessions:
            continue
        by_genus[genus].append(rec)
    if not by_genus or budget <= 0:
        return []

    best_per_genus = {
        genus: min(recs, key=lambda r: (-len(r["translation"]),
                                        r["accession"],
                                        r["translation"]))
        for genus, recs in by_genus.items()
    }
    ordered_genera = sorted(
        best_per_genus, key=lambda g: (-len(by_genus[g]), g))

    chosen = []
    for genus in ordered_genera[:budget]:
        rec = best_per_genus[genus]
        chosen.append(rec)
        used_genera.add(genus)
        used_accessions.add(rec["accession"])
    return chosen


def _distinct_genus_count(group_items: list) -> int:
    return len({genus_of(rec["organism"]) for rec, _rel in group_items})


def _fair_shares(active: list, remaining: int, offset: int) -> list:
    """
    Split `remaining` across `active` groups by max-min water-filling
    against each group's *distinct genus count* as its capacity, not a
    blind even split.

    A blind even split looked appealing (and matches §3's illustrative
    pseudocode) but real GenBank lineages are not the tidy
    class/order/family/genus ladder the pseudocode sketches — NCBI's
    taxonomy inserts many extra unranked clades, and Arthropoda's path
    to e.g. Coleoptera is ~10 lineage entries deep, mostly near-binary
    splits (Altocrustacea/Communostraca, Neoptera/Palaeoptera,
    Endopterygota/Paraneoptera, ...). Splitting evenly at every one of
    those halves the budget regardless of which side holds the
    diversity, so by depth 10 a blind 50/50 rule leaves a large,
    genuinely diverse clade with less than one representative's worth
    of budget — caught on the first live `animal_mt` build, where
    Coleoptera (438 distinct genera among the candidates) received zero
    of a 1,000-strong Arthropoda budget.

    Water-filling fixes this while still honouring "before a second
    beetle, take a mite" exactly: a small-capacity group is given
    *only* what it can use (its full genus capacity, no more), freeing
    the rest for groups that can use more, rather than pinning every
    group to an equal fraction regardless of capacity. Groups are
    processed smallest-capacity-first; once no remaining group's
    capacity is below the current equal share, whatever budget is left
    is split evenly across them, with the integer remainder going to
    the largest-capacity groups first (§3's own remainder rule) and
    `offset` rotating ties across repeated calls (the caller's
    under-fill redistribution loop) so a fixed sort order cannot
    perpetually starve the same tail of groups.
    """
    n = len(active)
    caps = [max(_distinct_genus_count(gi), 1) for _, gi in active]
    shares = [0] * n
    pending = list(range(n))
    pool = remaining
    while pending and pool > 0:
        equal_share = pool // len(pending)
        finished = [i for i in pending if caps[i] <= equal_share]
        if not finished:
            extra = pool - equal_share * len(pending)
            order = sorted(
                pending, key=lambda i: (-caps[i], (i - offset) % n))
            bonus = set(order[:extra])
            for i in pending:
                shares[i] = equal_share + (1 if i in bonus else 0)
            break
        for i in finished:
            shares[i] = caps[i]
            pool -= caps[i]
        pending = [i for i in pending if i not in finished]
    return shares


def select(items: list, budget: int, used_genera: set,
           used_accessions: set, depth: int = 0) -> list:
    """Rank-aware breadth-first traversal. `items` are (rec,
    rel_lineage) pairs already restricted to one stratum. Splits
    `budget` across the groups found at `depth` in proportion to each
    group's distinct genus count (`_fair_shares`), recurses one rank
    deeper into each, then reallocates any unspent budget from
    under-filled groups to groups that still have candidates — until
    the budget is spent or no group can make further progress."""
    if budget <= 0 or not items:
        return []

    groups = _partition_by_rank(items, depth)
    if len(groups) <= 1:
        return _take_one_per_genus(
            items, budget, used_genera, used_accessions)

    active = groups
    chosen = []
    remaining = budget
    offset = 0
    while remaining > 0 and active:
        shares = _fair_shares(active, remaining, offset)
        next_active = []
        spent = 0
        for (key, group_items), group_budget in zip(active, shares):
            if group_budget <= 0:
                next_active.append((key, group_items))
                continue
            picked = select(
                group_items, group_budget, used_genera, used_accessions,
                depth + 1)
            chosen.extend(picked)
            spent += len(picked)
            picked_acc = {r["accession"] for r in picked}
            leftover = [
                it for it in group_items
                if it[0]["accession"] not in picked_acc
            ]
            if leftover:
                next_active.append((key, leftover))
        remaining -= spent
        if spent == 0:
            # No group made progress this round (all remaining
            # candidates are exhausted or genus-capped out) — further
            # rounds would loop forever offering the same records.
            break
        active = next_active
        offset = offset % len(active) if active else 0
    return chosen


def select_representatives(candidates: list, strata: list,
                           min_per_locus: int) -> tuple:
    """
    Split `candidates` across ordered, mutually-exclusive strata
    (`[(lineage_term_or_None, budget), ...]`) and select representatives
    within each via `select()`. A record matching an earlier stratum's
    lineage term is never eligible for a later one, even if it was not
    itself selected — that is what keeps e.g. `Arthropoda` and
    `Metazoa` disjoint (§4.2).

    Returns `(representatives, stratum_stats)`. `[]` for both if fewer
    than `min_per_locus` candidates exist in total (the existing
    `--min-per-locus` floor, now applied across strata rather than to a
    single pool).
    """
    if len(candidates) < min_per_locus:
        return [], []

    # Real RefSeq data has duplicate (accession, gene) CDS records —
    # e.g. isoform reannotations of the same genome — with different
    # translations. A plain `{accession: rec}` dict comprehension would
    # silently keep whichever duplicate happens to appear last in
    # `candidates`, which depends on parse/file order, not content —
    # caught on the first live build (task 48: `animal_mt/ATP6` picked
    # a different one of two NC_009093.1 translations depending on
    # jsonl line order). Resolve duplicates the same deterministic way
    # as the genus-level tie-break: longest translation wins, then the
    # translation string itself as the final tie-break.
    remaining = {}
    best_key = {}
    for rec in candidates:
        acc = rec["accession"]
        key = (-len(rec["translation"]), rec["translation"])
        if acc not in remaining or key < best_key[acc]:
            remaining[acc] = rec
            best_key[acc] = key
    chosen = []
    stats = []
    prior_terms = []
    for term, budget in strata:
        if term is None:
            eligible_acc = list(remaining.keys())
        else:
            eligible_acc = [
                acc for acc, rec in remaining.items()
                if term in (rec.get("lineage") or [])
            ]
        eligible = [remaining[acc] for acc in eligible_acc]
        items = [(rec, _rel_lineage(rec, term)) for rec in eligible]

        picked = select(items, budget, set(), set())
        picked.sort(key=lambda r: r["accession"])
        chosen.extend(picked)

        for acc in eligible_acc:
            del remaining[acc]

        label = term if term is not None else "all"
        if prior_terms:
            label = f"{label}!{'!'.join(prior_terms)}"
        if term is not None:
            prior_terms.append(term)

        stats.append({
            "label": label,
            "budget": budget,
            "candidates": len(eligible),
            "selected": len(picked),
            "genera": len({genus_of(r["organism"]) for r in picked}),
        })
    return chosen, stats


def write_faa(path: Path, records: list):
    """
    Write one FASTA record per representative. The sequence ID is
    `<accession>_<gene>`, not the bare accession: a given source genome
    is often picked as the representative for more than one gene (it is
    diversity-per-gene, not diversity-per-genome), so bare accessions
    collide once every per-gene file is concatenated into a single
    query at runtime — and miniprot requires unique query IDs.
    """
    lines = []
    for rec in records:
        seq_id = f"{rec['accession']}_{rec['gene']}"
        header = f">{seq_id} {rec['gene']} [{rec['organism']}]"
        lines.append(header)
        lines.append(rec["translation"])
    path.write_text("\n".join(lines) + "\n")


def parse_strata(raw_strata: list) -> dict:
    """
    Parse repeatable `--stratum ORIGIN:TERM=BUDGET` flags into
    `{origin: [(term, budget), ...]}`, preserving declaration order
    per origin (§4.2 — order determines the mutual-exclusion chain).
    """
    strata = defaultdict(list)
    for raw in raw_strata:
        origin, sep1, rest = raw.partition(":")
        term, sep2, budget_s = rest.partition("=")
        if not (sep1 and sep2 and origin and term and budget_s):
            print(
                f"ERROR: --stratum expects ORIGIN:TERM=BUDGET, got "
                f"{raw!r}",
                file=sys.stderr,
            )
            sys.exit(1)
        try:
            budget = int(budget_s)
        except ValueError:
            print(
                f"ERROR: --stratum budget must be an integer, got "
                f"{raw!r}",
                file=sys.stderr,
            )
            sys.exit(1)
        strata[origin].append((term, budget))
    return strata


def main_full(args):
    gene_sets = load_gene_sets(args.gene_sets)
    out_root = Path(args.out)
    staging = stage_dir(out_root)
    strata_by_origin = parse_strata(args.stratum)

    cds_sources = {}
    for entry in args.cds_jsonl:
        origin, _, path = entry.partition("=")
        if not origin or not path:
            print(
                f"ERROR: --cds-jsonl expects ORIGIN=PATH, got {entry!r}",
                file=sys.stderr,
            )
            sys.exit(1)
        cds_sources[origin] = path

    manifest_genes = {}
    strata_config = {}
    errors = []

    for origin, gene_set in gene_sets.items():
        if origin not in cds_sources:
            print(
                f"[proteins] WARNING: no --cds-jsonl for '{origin}' — "
                "skipping",
                file=sys.stderr,
            )
            continue

        origin_strata = strata_by_origin.get(origin) or [
            (None, args.max_per_locus)
        ]
        strata_config[origin] = [
            {"term": term, "budget": budget}
            for term, budget in origin_strata
        ]

        alias_index = build_alias_index(gene_set)
        by_canonical = defaultdict(list)
        with open(cds_sources[origin]) as fh:
            for line in fh:
                rec = json.loads(line)
                canonical = alias_index.get(normalise_symbol(rec["gene"]))
                if canonical is not None:
                    by_canonical[canonical].append(rec)

        origin_dir = staging / origin
        origin_dir.mkdir(exist_ok=True)
        origin_genes = {}

        for gene in gene_set["protein_coding"]:
            candidates = by_canonical.get(gene, [])
            reps, stratum_stats = select_representatives(
                candidates, origin_strata, args.min_per_locus)
            print(
                f"[proteins] {origin}/{gene}: {len(candidates)} candidates "
                f"-> {len(reps)} representatives",
                file=sys.stderr,
            )
            if not reps:
                msg = (
                    f"{origin}/{gene}: only {len(candidates)} candidate "
                    f"CDS translations found (min={args.min_per_locus})"
                )
                errors.append(msg)
                continue
            write_faa(origin_dir / f"{gene}.faa", reps)
            origin_genes[gene] = {
                "candidates": len(candidates),
                "representatives": len(reps),
                "strata": {
                    s["label"]: {
                        "budget": s["budget"],
                        "selected": s["selected"],
                        "genera": s["genera"],
                    }
                    for s in stratum_stats
                },
                "accessions": [r["accession"] for r in reps],
            }

        manifest_genes[origin] = origin_genes

    if errors:
        print("\n[proteins] ERRORS — build incomplete:", file=sys.stderr)
        for e in errors:
            print(f"  {e}", file=sys.stderr)
        sys.exit(1)

    provenance = {
        "refseq_release": args.refseq_release,
        "source": (
            "CDS translations parsed from bundled RefSeq organelle "
            "GenBank flat files — see scripts/refdata/parse_gbff_cds.py"
        ),
        "min_per_locus": args.min_per_locus,
        "max_per_locus": args.max_per_locus,
        "strata_config": strata_config,
        "genes": manifest_genes,
    }
    (staging / "provenance.json").write_text(
        json.dumps(provenance, indent=2) + "\n")

    promote(staging, out_root)
    print(f"[proteins] done: {out_root}")


def main():
    args = parse_args()
    if args.mode == "barcode":
        main_barcode(args)
    else:
        main_full(args)


if __name__ == "__main__":
    main()
