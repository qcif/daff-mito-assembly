#!/usr/bin/env python3
"""
Split RefSeq mitochondrion FASTA by kingdom.

Reads a (concatenated, uncompressed) RefSeq mitochondrion FASTA, looks up
each accession's taxid via NCBI eutils, then uses taxonkit to determine
whether each record belongs to Metazoa (txid 33208) or Viridiplantae
(txid 33090). Writes two output FASTAs; records outside those kingdoms are
logged to stderr and discarded.

Usage:
    python3 split_refseq_mt.py \
        --fasta  refseq_mt.fa \
        --out-metazoa    refseq_mt_metazoa.fa \
        --out-viridiplantae refseq_mt_viridiplantae.fa \
        --taxonkit /home/cameron/.local/bin/taxonkit \
        --taxdump  /home/cameron/.taxonkit/new_v0.20 \
        [--ncbi-api-key KEY] \
        [--batch-size 500]
"""

import argparse
import subprocess
import sys
import time
import urllib.request
import urllib.parse
import json

EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

METAZOA_TXID = 33208
VIRIDIPLANTAE_TXID = 33090


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--fasta", required=True)
    p.add_argument("--out-metazoa", required=True)
    p.add_argument("--out-viridiplantae", required=True)
    p.add_argument("--taxonkit", default="taxonkit")
    p.add_argument("--taxdump", required=True)
    p.add_argument("--ncbi-api-key", default="")
    p.add_argument("--batch-size", type=int, default=500)
    return p.parse_args()


def iter_fasta_records(path):
    """Yield (header_line, sequence_block) tuples."""
    header = None
    seq_lines = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    yield header, "\n".join(seq_lines)
                header = line
                seq_lines = []
            else:
                seq_lines.append(line)
    if header is not None:
        yield header, "\n".join(seq_lines)


def extract_accession(header):
    """Return versioned accession (e.g. NC_001628.1) from a FASTA header."""
    return header[1:].split()[0]


def fetch_taxids_batch(accessions, api_key="", rate_delay=0.12):
    """
    Return {accession: taxid_int} for a list of accessions.
    Uses NCBI esummary via POST to avoid HTTP 414 on large batches.
    Sleeps between requests to respect rate limits
    (3 req/sec without key → 0.40s; 10 req/sec with key → 0.12s).
    """
    result = {}
    ids_str = ",".join(accessions)
    params = {"db": "nucleotide", "id": ids_str, "retmode": "json"}
    if api_key:
        params["api_key"] = api_key
    body = urllib.parse.urlencode(params).encode()
    req = urllib.request.Request(
        f"{EUTILS_BASE}/esummary.fcgi",
        data=body,
        headers={"Content-Type": "application/x-www-form-urlencoded"},
    )
    try:
        with urllib.request.urlopen(req, timeout=120) as resp:
            data = json.loads(resp.read())
        uid_list = data.get("result", {}).get("uids", [])
        for uid in uid_list:
            rec = data["result"][uid]
            acc = rec.get("accessionversion", rec.get("caption", ""))
            taxid = int(rec.get("taxid", 0))
            if acc and taxid:
                result[acc] = taxid
    except Exception as exc:
        print(f"  WARNING: esummary failed for batch: {exc}", file=sys.stderr)
    time.sleep(rate_delay)
    return result


def fetch_all_taxids(accessions, api_key="", batch_size=500):
    """Batch-fetch taxids for all accessions. Returns {accession: taxid}."""
    mapping = {}
    batches = [
        accessions[i:i+batch_size]
        for i in range(0, len(accessions), batch_size)
    ]
    print(
        f"[split_refseq_mt] fetching taxids for {len(accessions)} accessions "
        f"in {len(batches)} batches",
        file=sys.stderr,
    )
    for i, batch in enumerate(batches, 1):
        if i % 10 == 0 or i == len(batches):
            print(f"  batch {i}/{len(batches)}", file=sys.stderr)
        delay = 0.12 if api_key else 0.40
        mapping.update(
            fetch_taxids_batch(batch, api_key=api_key, rate_delay=delay)
        )
    return mapping


def get_kingdom_set(taxids, taxonkit_bin, taxdump_dir, target_txid):
    """
    Return the set of taxids that are descendants-of (or equal to) target_txid.
    Uses: taxonkit list --ids <target_txid> --indent "" | sort
    """
    result = subprocess.run(
        [
            taxonkit_bin,
            "list",
            "--data-dir", taxdump_dir,
            "--ids", str(target_txid),
            "--indent", "",
        ],
        capture_output=True,
        text=True,
        check=True,
    )
    return {int(t) for t in result.stdout.split() if t.strip()}


def main():
    args = parse_args()
    fasta_path = args.fasta

    print(
        f"[split_refseq_mt] reading headers from {fasta_path}",
        file=sys.stderr,
    )
    records = list(iter_fasta_records(fasta_path))
    print(f"[split_refseq_mt] {len(records)} records found", file=sys.stderr)

    accessions = [extract_accession(h) for h, _ in records]
    acc_to_taxid = fetch_all_taxids(
        accessions,
        api_key=args.ncbi_api_key,
        batch_size=args.batch_size,
    )

    missing = [a for a in accessions if a not in acc_to_taxid]
    if missing:
        print(
            f"[split_refseq_mt] WARNING: {len(missing)} accessions "
            f"had no taxid (first 5: {missing[:5]})",
            file=sys.stderr,
        )

    all_taxids = list({t for t in acc_to_taxid.values() if t})

    print("[split_refseq_mt] building Metazoa taxid set ...", file=sys.stderr)
    metazoa_set = get_kingdom_set(
        all_taxids, args.taxonkit, args.taxdump, METAZOA_TXID
    )
    print(
        f"[split_refseq_mt] Metazoa subtree: {len(metazoa_set)} taxids",
        file=sys.stderr,
    )

    print(
        "[split_refseq_mt] building Viridiplantae taxid set ...",
        file=sys.stderr,
    )
    plant_set = get_kingdom_set(
        all_taxids, args.taxonkit, args.taxdump, VIRIDIPLANTAE_TXID
    )
    print(
        f"[split_refseq_mt] Viridiplantae subtree: {len(plant_set)} taxids",
        file=sys.stderr,
    )

    n_metazoa = n_plant = n_other = 0
    with (
        open(args.out_metazoa, "w") as fh_meta,
        open(args.out_viridiplantae, "w") as fh_plant,
    ):
        for header, seq in records:
            acc = extract_accession(header)
            taxid = acc_to_taxid.get(acc, 0)
            if taxid in metazoa_set:
                fh_meta.write(f"{header}\n{seq}\n")
                n_metazoa += 1
            elif taxid in plant_set:
                fh_plant.write(f"{header}\n{seq}\n")
                n_plant += 1
            else:
                print(
                    f"[split_refseq_mt] discarding {acc} taxid={taxid} "
                    f"(not Metazoa or Viridiplantae)",
                    file=sys.stderr,
                )
                n_other += 1

    print(
        f"[split_refseq_mt] done — metazoa={n_metazoa}, "
        f"viridiplantae={n_plant}, discarded={n_other}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
