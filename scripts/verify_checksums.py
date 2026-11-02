#!/usr/bin/env python3
"""Verify MD5 checksums of pipeline input files.

The checksum file uses the standard ``md5sum`` format (``<md5>  <path>``);
blank lines and ``#`` comments are ignored. Unlike ``md5sum -c`` this check is
*closed*: every file passed via ``--files`` MUST be listed, otherwise it fails,
so an input can never silently skip verification.

Usage:
    python scripts/verify_checksums.py --checksums data/raw/checksums.md5 \
        --files data/raw/genome.fasta data/raw/gene.fasta
"""
import argparse
import hashlib
import os
import sys

from utils.runtime import event, get_logger, init_script, run_main

log = get_logger("verify_checksums")


def md5_of(path, chunk=1 << 20):
    digest = hashlib.md5()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(chunk), b""):
            digest.update(block)
    return digest.hexdigest()


def read_checksums(path):
    entries = {}
    with open(path) as handle:
        for lineno, line in enumerate(handle, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split(None, 1)
            if len(parts) != 2 or len(parts[0]) != 32:
                raise ValueError(f"{path}:{lineno}: not in 'md5  path' format")
            name = parts[1].lstrip("*")          # md5sum binary-mode marker
            entries[os.path.normpath(name)] = parts[0].lower()
    return entries


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--checksums", required=True)
    parser.add_argument("--files", nargs="+", required=True)
    args = parser.parse_args()
    init_script(log, args, None, {"checksums": args.checksums,
                                  **{f"file{i}": f for i, f in enumerate(args.files)}})

    try:
        expected = read_checksums(args.checksums)
    except (OSError, ValueError) as err:
        log.error("%s", err)
        sys.exit(1)

    failures = []
    for path in args.files:
        key = os.path.normpath(path)
        if key not in expected:
            failures.append(f"{path}: no entry in {args.checksums} "
                            f"(add one with: md5sum {path} >> {args.checksums})")
            continue
        actual = md5_of(path)
        if actual != expected[key]:
            failures.append(f"{path}: MD5 mismatch "
                            f"(expected {expected[key]}, got {actual})")
        else:
            event(log, "INFO", "checksum OK", path=path, md5=actual)

    if failures:
        for msg in failures:
            log.error("%s", msg)
        sys.exit(1)
    event(log, "INFO", "all checksums match", n_files=len(args.files),
          manifest=args.checksums)


if __name__ == "__main__":
    run_main(main, "verify_checksums")
