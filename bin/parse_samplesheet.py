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
from urllib.parse import urlparse

REQUIRED_COLUMNS = {'sample_id', 'assembly_target', 'reads'}
OPTIONAL_COLUMNS = {
    'sample_info',
    'sample_type',
    'sample_receipt_date',
    'storage_location',
}
KNOWN_COLUMNS = REQUIRED_COLUMNS | OPTIONAL_COLUMNS
VALID_ASSEMBLY_TARGETS = {'animal_mt', 'plant_pt', 'plant_mt'}
SAMPLE_ID_RE = re.compile(r'^[A-Za-z0-9_.-]+$')
ALLOWED_URI_SCHEMES = ('az://',)
URI_SCHEME_RE = re.compile(r'^([A-Za-z][A-Za-z0-9+.-]*)://')
READ_KIND_RELATIVE = 'relative'
READ_KIND_URI = 'uri'
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
    parser.add_argument(
        '--data-dir-set', action='store_true', dest='data_dir_set',
        help='--data_dir was supplied to the workflow (contents not read)')
    parser.add_argument('--out', required=True)
    return parser.parse_args()


def _validate_header(fieldnames, errors):
    fields = set(fieldnames)
    for col in sorted(REQUIRED_COLUMNS - fields):
        errors.append(f'Missing required column: {col!r}')
    for col in sorted(fields - KNOWN_COLUMNS):
        errors.append(f'Unknown column: {col!r} (typo?)')


def _classify_read(part):
    """Return (kind, error) for one reads entry.

    kind is 'relative' or 'uri'; error is None when the entry is usable.
    """
    scheme = URI_SCHEME_RE.match(part)
    if scheme:
        if scheme.group(1).lower() + '://' not in ALLOWED_URI_SCHEMES:
            return None, (
                f'unsupported URI scheme in {part!r}'
                f' — allowed: {", ".join(sorted(ALLOWED_URI_SCHEMES))}'
            )
        parsed = urlparse(part)
        if not parsed.netloc or not parsed.path.strip('/'):
            return None, (
                f'malformed URI {part!r}'
                f' — expected <scheme>://<container>/<path>'
            )
        return READ_KIND_URI, None
    if part.startswith(('/', '~')):
        return None, (
            f'local absolute path not allowed in reads: {part!r}'
            f' — use a path relative to --data_dir, or an az:// URI'
        )
    return READ_KIND_RELATIVE, None


def _parse_reads(reads_raw, label, errors, kinds):
    """Validate a reads cell, returning its entries verbatim.

    Each accepted entry's kind is recorded in kinds (kind -> labels).
    """
    if not reads_raw:
        errors.append(f'{label}: reads column is empty')
        return []

    entries = []
    for part in reads_raw.split('|'):
        part = part.strip()
        if not part:
            errors.append(
                f'{label}: empty entry in reads (check for stray pipes)'
            )
            continue
        kind, error = _classify_read(part)
        if error:
            errors.append(f'{label}: {error}')
            continue
        if not _has_fastq_ext(urlparse(part).path
                              if kind == READ_KIND_URI else part):
            errors.append(
                f'{label}: unsupported extension in {part!r}'
                f' — expected .fastq, .fastq.gz, .fq, .fq.gz'
            )
            continue
        kinds.setdefault(kind, []).append(label)
        entries.append(part)

    return entries


def _check_read_sources(kinds, data_dir_set, errors, warnings):
    """Sheet-level rules: no mixing, and --data_dir only where needed."""
    if len(kinds) > 1:
        errors.append(
            'reads mix relative paths (need --data_dir) and az:// URIs'
            ' (no --data_dir) — use one kind per sheet. Relative:'
            f' {", ".join(sorted(set(kinds[READ_KIND_RELATIVE])))}.'
            f' URI: {", ".join(sorted(set(kinds[READ_KIND_URI])))}'
        )
    elif READ_KIND_RELATIVE in kinds and not data_dir_set:
        errors.append(
            '--data_dir is required for relative reads paths. Rows:'
            f' {", ".join(sorted(set(kinds[READ_KIND_RELATIVE])))}'
        )
    elif READ_KIND_URI in kinds and data_dir_set:
        warnings.append(
            '--data_dir was supplied but every reads entry is a URI;'
            ' it will be unused'
        )


def _validate_row(row, row_num, seen_ids, errors, warnings, kinds):
    sample_id = row.get('sample_id', '').strip()
    assembly_target_raw = row.get('assembly_target', '').strip()
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

    assembly_target = assembly_target_raw.lower()
    if assembly_target not in VALID_ASSEMBLY_TARGETS:
        errors.append(
            f'{label}: invalid assembly_target {assembly_target_raw!r}'
            f' — must be one of:'
            f' {", ".join(sorted(VALID_ASSEMBLY_TARGETS))}'
        )

    reads = _parse_reads(reads_raw, label, errors, kinds)

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
        'assembly_target': assembly_target,
        'reads': reads,
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
    kinds = {}

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
                _validate_row(row, i, seen_ids, errors, warnings, kinds)
            )

    _check_read_sources(kinds, args.data_dir_set, errors, warnings)

    for w in warnings:
        print(f'WARNING: {w}', file=sys.stderr)

    if errors:
        _fail(errors)

    with open(args.out, 'w') as fh:
        json.dump(rows, fh, indent=2)


if __name__ == '__main__':  # pragma: no cover
    main()
