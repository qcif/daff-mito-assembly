#!/usr/bin/env python3
"""Verify SHA256 checksums of fetched integration fixtures against a pinned manifest."""

import hashlib
import sys
from pathlib import Path

BLOCK_SIZE = 1 << 20  # 1 MB


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        while chunk := fh.read(BLOCK_SIZE):
            h.update(chunk)
    return h.hexdigest()


def load_manifest(manifest_path: Path) -> dict[str, str]:
    manifest = {}
    for line in manifest_path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        digest, filename = line.split(None, 1)
        manifest[filename.strip()] = digest.strip()
    return manifest


def main() -> None:
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <fixtures_dir> <manifest.sha256>", file=sys.stderr)
        sys.exit(1)

    fixtures_dir = Path(sys.argv[1])
    manifest_path = Path(sys.argv[2])

    manifest = load_manifest(manifest_path)
    errors = []

    for filename, expected in manifest.items():
        path = fixtures_dir / filename
        if not path.exists():
            errors.append(f"MISSING: {filename}")
            continue
        actual = sha256_file(path)
        if actual != expected:
            errors.append(
                f"CHECKSUM MISMATCH: {filename}\n"
                f"  expected: {expected}\n"
                f"  got:      {actual}"
            )
        else:
            print(f"OK: {filename}")

    if errors:
        for e in errors:
            print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"All {len(manifest)} fixtures verified.")


if __name__ == '__main__':
    main()
