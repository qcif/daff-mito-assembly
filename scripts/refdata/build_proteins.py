#!/usr/bin/env python3
"""
Build refs/<version>/proteins/<origin>/<gene>.faa.

Two modes:

  barcode (legacy) — for each (origin, gene) entry in assets/loci.json
    (schema wf5/loci-panel/v1), queries NCBI RefSeq protein via
    POST-based esearch → efetch (history/WebEnv, to avoid HTTP 414 on
    long ID batches) and writes one FASTA per locus. Superseded by
    `full` (task 29) but kept for reference / re-derivation.

  full (task 29) — builds the organelle's comprehensive protein-coding
    complement from CDS translations already present in the reference
    bundle's RefSeq organelle records (validate/refseq_pt.fa,
    validate/refseq_mt_{metazoa,viridiplantae}.fa), rather than issuing
    fresh per-gene NCBI queries. Reads a CDS jsonl file produced by
    scripts/refdata/parse_gbff_cds.py (run inside the
    neoformit/daff-wf5-scripts:test image, which has biopython) and
    assets/organelle_gene_sets.json for the canonical gene list +
    alias table, then picks 5-10 phylogenetically-spread
    representatives per gene (diversity proxy: distinct genus, taken
    from the GenBank /organism qualifier).

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
        --out             refs/v2026.09/proteins \
        --refseq-release  236 \
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


def pick_representatives(records: list, min_n: int, max_n: int) -> list:
    """
    Pick up to max_n records maximising genus spread: iterate distinct
    genera in first-seen order, taking one record per pass, until max_n
    is reached or genera are exhausted (then top up from remaining
    records). Returns [] if fewer than min_n records are available.
    """
    if len(records) < min_n:
        return []
    by_genus = defaultdict(list)
    for rec in records:
        by_genus[genus_of(rec["organism"])].append(rec)

    chosen = []
    seen_acc = set()
    while len(chosen) < max_n and any(by_genus.values()):
        for genus in list(by_genus.keys()):
            if len(chosen) >= max_n:
                break
            bucket = by_genus[genus]
            if not bucket:
                del by_genus[genus]
                continue
            rec = bucket.pop(0)
            if rec["accession"] not in seen_acc:
                chosen.append(rec)
                seen_acc.add(rec["accession"])
            if not bucket:
                del by_genus[genus]
    return chosen


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


def main_full(args):
    gene_sets = load_gene_sets(args.gene_sets)
    out_root = Path(args.out)
    staging = stage_dir(out_root)

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
    errors = []

    for origin, gene_set in gene_sets.items():
        if origin not in cds_sources:
            print(
                f"[proteins] WARNING: no --cds-jsonl for '{origin}' — "
                "skipping",
                file=sys.stderr,
            )
            continue

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
            reps = pick_representatives(
                candidates, args.min_per_locus, args.max_per_locus)
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
