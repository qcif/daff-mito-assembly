#!/usr/bin/env python3
"""
Build refs/<version>/proteins/<origin>/<gene>.faa from assets/loci.json
(schema wf5/loci-panel/v1). See tasks/3_refdata.md §2.3.

For each (origin, gene) entry in the panel, queries NCBI RefSeq protein
via POST-based esearch → efetch (history/WebEnv, to avoid HTTP 414 on
long ID batches) and writes one FASTA per locus.

Origin → NCBI taxid restriction:
    animal_mt → txid33208 (Metazoa)
    plant_pt  → txid33090 (Viridiplantae)
    plant_mt  → txid33090 (Viridiplantae)

Requires: Python 3.10+, urllib (stdlib). No third-party dependencies.

Usage:
    python3 build_proteins.py \
        --loci-config  assets/loci.json \
        --out          refs/v2026.07/proteins \
        [--ncbi-api-key KEY] \
        [--max-per-locus 50]
"""

import argparse
import json
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
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
    p.add_argument("--loci-config", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--ncbi-api-key", default="")
    p.add_argument("--max-per-locus", type=int, default=50)
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


def main():
    args = parse_args()
    loci_config = load_loci(args.loci_config)
    out_root = Path(args.out)
    staging = Path(str(out_root) + ".staging")
    # Wipe any residue from a prior failed run — stale sub-dirs would
    # otherwise leak through the atomic promote below.
    import shutil
    if staging.exists():
        shutil.rmtree(staging)
    staging.mkdir(parents=True)

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

    # Atomic promote
    if out_root.exists():
        shutil.rmtree(out_root)
    out_root.parent.mkdir(parents=True, exist_ok=True)
    staging.rename(out_root)
    print(f"[proteins] done: {out_root}")


if __name__ == "__main__":
    main()
