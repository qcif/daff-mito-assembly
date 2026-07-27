#!/usr/bin/env python3
"""C1: Validate wf5 sample sheet and emit a normalised JSON array.

Exit 0 on success; exits non-zero with all accumulated errors on stderr.
"""

import argparse
import csv
import json
import re
import sys
from datetime import date
from pathlib import Path

REQUIRED_COLUMNS = {'sample_id', 'kingdom', 'reads'}
OPTIONAL_COLUMNS = {
    'sample_info',
    'sample_type',
    'sample_receipt_date',
    'storage_location',
}
KNOWN_COLUMNS = REQUIRED_COLUMNS | OPTIONAL_COLUMNS
VALID_KINGDOMS = {'plant', 'animal'}
SAMPLE_ID_RE = re.compile(r'^[A-Za-z0-9_.-]+$')
VALID_EXTENSIONS = ('.fastq', '.fastq.gz', '.fq', '.fq.gz')


def _has_fastq_ext(path_str):
    p = path_str.lower()
    return any(p.endswith(ext) for ext in VALID_EXTENSIONS)


def _fail(errors):
    for err in errors:
        print(f'ERROR: {err}', file=sys.stderr)
    sys.exit(1)


def parse_args():
    parser = argparse.ArgumentParser(
        description='Validate wf5 sample sheet and emit normalised JSON.')
    parser.add_argument('--samplesheet', required=True)
    parser.add_argument('--data-dir', required=True, dest='data_dir')
    parser.add_argument('--out', required=True)
    return parser.parse_args()


def _validate_header(fieldnames, errors):
    fields = set(fieldnames)
    for col in sorted(REQUIRED_COLUMNS - fields):
        errors.append(f'Missing required column: {col!r}')
    for col in sorted(fields - KNOWN_COLUMNS):
        errors.append(f'Unknown column: {col!r} (typo?)')


def _resolve_reads(reads_raw, label, data_dir, errors):
    if not reads_raw:
        errors.append(f'{label}: reads column is empty')
        return []

    resolved = []
    for part in reads_raw.split('|'):
        part = part.strip()
        if not part:
            errors.append(
                f'{label}: empty entry in reads (check for stray pipes)'
            )
            continue
        if part.startswith('/'):
            errors.append(
                f'{label}: absolute path not allowed in reads'
                f' — use --data-dir root: {part!r}'
            )
            continue
        if not _has_fastq_ext(part):
            errors.append(
                f'{label}: unsupported extension in {part!r}'
                f' — expected .fastq, .fastq.gz, .fq, .fq.gz'
            )
            continue
        full = (data_dir / part).resolve()
        if not full.exists():
            errors.append(f'{label}: reads file not found: {full}')
            continue
        resolved.append(str(full))

    return resolved


def _validate_row(row, row_num, data_dir, seen_ids, errors, warnings):
    sample_id = row.get('sample_id', '').strip()
    kingdom_raw = row.get('kingdom', '').strip()
    reads_raw = row.get('reads', '').strip()
    date_val = row.get('sample_receipt_date', '').strip()
    label = f'row {row_num} (sample_id={sample_id!r})'

    if not SAMPLE_ID_RE.match(sample_id):
        errors.append(
            f'{label}: invalid sample_id {sample_id!r}'
            f' — must match [A-Za-z0-9_.-]+'
        )

    if sample_id in seen_ids:
        errors.append(
            f'{label}: duplicate sample_id {sample_id!r}'
            f' (first seen at row {seen_ids[sample_id]})'
        )
    else:
        seen_ids[sample_id] = row_num

    kingdom = kingdom_raw.lower()
    if kingdom not in VALID_KINGDOMS:
        errors.append(
            f'{label}: invalid kingdom {kingdom_raw!r}'
            f' — must be one of: {", ".join(sorted(VALID_KINGDOMS))}'
        )

    resolved = _resolve_reads(reads_raw, label, data_dir, errors)

    if date_val:
        try:
            date.fromisoformat(date_val)
        except ValueError:
            warnings.append(
                f'{label}: sample_receipt_date {date_val!r}'
                f' is not ISO 8601 (YYYY-MM-DD) — passing through raw'
            )

    return {
        'sample_id': sample_id,
        'kingdom': kingdom,
        'reads': resolved,
        'sample_info': row.get('sample_info', '').strip(),
        'sample_type': row.get('sample_type', '').strip(),
        'sample_receipt_date': date_val,
        'storage_location': row.get('storage_location', '').strip(),
    }


def main():
    args = parse_args()
    errors = []
    warnings = []
    rows = []
    data_dir = Path(args.data_dir)

    with open(args.samplesheet, newline='') as fh:
        reader = csv.DictReader(fh)

        if not reader.fieldnames:
            _fail(['Sample sheet is empty or has no header row'])

        _validate_header(reader.fieldnames, errors)
        if errors:
            _fail(errors)

        seen_ids = {}
        for i, row in enumerate(reader, start=2):
            rows.append(
                _validate_row(row, i, data_dir, seen_ids, errors, warnings)
            )

    for w in warnings:
        print(f'WARNING: {w}', file=sys.stderr)

    if errors:
        _fail(errors)

    with open(args.out, 'w') as fh:
        json.dump(rows, fh, indent=2)


if __name__ == '__main__':
    main()
