#!/usr/bin/env python3
"""Validate FASTA input files before the pipeline runs.

Checks, with actionable messages that name the file, line and column:
  * FASTA syntax: first record starts with '>', every header has an ID, IDs
    are unique, no record is empty, no whitespace inside a sequence line;
  * alphabet: target gene must be strictly ACGT (case-insensitive); the genome
    may also contain IUPAC ambiguity codes, up to
    ``validation.max_genome_ambiguous_fraction``;
  * length: every target-gene record >= ``validation.min_gene_length`` and the
    genome total >= ``validation.min_genome_length``.

All problems in a file are collected and reported together; exit code 1 if any
check fails. Paths and thresholds come from config/config.yaml.
"""
import argparse
import os
import sys

from utils.config import load_config
from utils.runtime import event, get_logger, init_script, run_main

log = get_logger("01_validate_inputs")

ACGT = frozenset("ACGT")
IUPAC = frozenset("ACGTRYSWKMBDHVN")
MAX_REPORTED = 8          # cap per-file error spam; the count is still exact


class Problems:
    """Collects error messages for one input file."""

    def __init__(self, path):
        self.path, self.items, self.total = path, [], 0

    def add(self, msg):
        self.total += 1
        if len(self.items) < MAX_REPORTED:
            self.items.append(msg)


def _hint(char):
    if char in "Uu":
        return " (RNA 'U' found - convert to DNA 'T')"
    if char.isdigit():
        return " (digits are not allowed - remove line/position numbers)"
    if char in "*-.":
        return f" (gap/stop symbol '{char}' - remove it or use a gap-free sequence)"
    return ""


def scan_fasta(path, allowed, problems):
    """Stream ``path``; returns [(id, length, n_non_acgt)] for each record."""
    records, seen_ids = [], {}
    cur = None                                   # [id, header_line, length, non_acgt]
    try:
        handle = open(path, encoding="utf-8")
    except OSError as err:
        problems.add(f"cannot open file: {err}")
        return records

    def close_record():
        if cur is not None:
            if cur[2] == 0:
                problems.add(f"line {cur[1]}: record '{cur[0]}' has no sequence "
                             "(header followed directly by another header or EOF)")
            records.append((cur[0], cur[2], cur[3]))

    with handle:
        try:
            for lineno, raw in enumerate(handle, start=1):
                line = raw.rstrip("\r\n")
                if not line.strip():
                    continue
                if line.startswith(">"):
                    close_record()
                    header = line[1:].split()
                    rec_id = header[0] if header else ""
                    if not rec_id:
                        problems.add(f"line {lineno}: header has no ID "
                                     "(expected '>some_id optional description')")
                    elif rec_id in seen_ids:
                        problems.add(f"line {lineno}: duplicate record ID '{rec_id}' "
                                     f"(first seen on line {seen_ids[rec_id]}); "
                                     "IDs must be unique")
                    else:
                        seen_ids[rec_id] = lineno
                    cur = [rec_id, lineno, 0, 0]
                    continue
                if cur is None:
                    problems.add(f"line {lineno}: sequence data before the first "
                                 "'>' header - is this really a FASTA file?")
                    cur = ["<no header>", lineno, 0, 0]
                if line != line.strip():
                    if line.rstrip() != line:
                        line = line.rstrip()     # trailing blanks are harmless
                    if line != line.strip():
                        problems.add(f"line {lineno}: leading whitespace in sequence line")
                        line = line.strip()
                if any(c.isspace() for c in line):
                    problems.add(f"line {lineno}: whitespace inside sequence line; "
                                 "remove spaces/tabs")
                    line = "".join(line.split())
                seq = line.upper()
                if set(seq) - allowed:
                    offenders = [(col, ch) for col, ch in enumerate(seq, start=1)
                                 if ch not in allowed]
                    for col, ch in offenders[:3]:
                        problems.add(f"line {lineno}, col {col}: invalid character "
                                     f"'{ch}'{_hint(ch)}; allowed: "
                                     f"{''.join(sorted(allowed))}")
                    if len(offenders) > 3:
                        problems.add(f"line {lineno}: {len(offenders) - 3} more "
                                     "invalid character(s) on this line")
                cur[2] += len(seq)
                cur[3] += len(seq) - sum(seq.count(b) for b in "ACGT")
        except UnicodeDecodeError as err:
            problems.add(f"not a UTF-8 text file ({err.reason} at byte {err.start}); "
                         "is it compressed or binary? Decompress it first")
            return records
    close_record()
    if not records and not problems.total:
        problems.add("no FASTA records found (file is empty or has no '>' header)")
    return records


def validate_gene(path, min_len):
    problems = Problems(path)
    records = scan_fasta(path, ACGT, problems)
    for rec_id, length, _ in records:
        if 0 < length < min_len and rec_id != "<no header>":
            problems.add(f"record '{rec_id}' is {length} bp, below the minimum "
                         f"target length of {min_len} bp "
                         "(validation.min_gene_length in config.yaml)")
    return records, problems


def validate_genome(path, min_len, max_amb_frac):
    problems = Problems(path)
    records = scan_fasta(path, IUPAC, problems)
    total = sum(r[1] for r in records)
    amb = sum(r[2] for r in records)
    if records and total < min_len:
        problems.add(f"genome totals {total:,} bp, below the minimum of {min_len:,} bp "
                     "(validation.min_genome_length in config.yaml) - wrong or "
                     "truncated file?")
    if total and amb / total > max_amb_frac:
        problems.add(f"{amb:,} of {total:,} bases ({amb / total:.2%}) are ambiguity "
                     f"codes, above the allowed {max_amb_frac:.2%} "
                     "(validation.max_genome_ambiguous_fraction)")
    elif amb:
        event(log, "WARNING", "genome contains ambiguity codes",
              path=path, ambiguous_bases=amb, fraction=round(amb / total, 6))
    return records, problems


def report(label, path, records, problems):
    if problems.total:
        event(log, "ERROR", f"{label} FAILED validation", path=path,
              n_problems=problems.total)
        for msg in problems.items:
            log.error("%s: %s", path, msg)
        if problems.total > len(problems.items):
            log.error("%s: ... and %d more problem(s)", path,
                      problems.total - len(problems.items))
        return False
    event(log, "INFO", f"{label} OK", path=path, n_records=len(records),
          total_bp=sum(r[1] for r in records))
    return True


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--config", default="config/config.yaml")
    parser.add_argument("--only", choices=["genome", "gene", "both"],
                        default="both",
                        help="which input to validate (default: both)")
    parser.add_argument("--genome", default=None,
                        help="override genome_fasta from the config")
    parser.add_argument("--gene", default=None,
                        help="override target_gene_fasta (per-gene wildcard)")
    args = parser.parse_args()

    cfg = load_config(args.config)
    val = cfg["validation"]
    genome = args.genome or cfg["genome_fasta"]
    gene = args.gene or cfg["target_gene_fasta"]
    inputs = {}
    if args.only in ("genome", "both"):
        inputs["genome"] = genome
    if args.only in ("gene", "both"):
        inputs["gene"] = gene
    init_script(log, args, cfg, inputs)

    ok = True
    if "genome" in inputs:
        recs, probs = validate_genome(genome, val["min_genome_length"],
                                      val["max_genome_ambiguous_fraction"])
        ok &= report("Genome", genome, recs, probs)
    if "gene" in inputs:
        recs, probs = validate_gene(gene, val["min_gene_length"])
        ok &= report("Target gene", gene, recs, probs)

    if not ok:
        sys.exit(1)
    log.info("Inputs look valid. Proceed to PAM scanning.")


if __name__ == "__main__":
    run_main(main, "01_validate_inputs")
