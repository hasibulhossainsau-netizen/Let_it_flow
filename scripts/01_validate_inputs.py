#!/usr/bin/env python3
"""
Validate FASTA input files before pipeline runs.
"""
from Bio import SeqIO
import sys, os

def validate_fasta(path, label, min_len=100):
    records = list(SeqIO.parse(path, "fasta"))
    if not records:
        print(f"[ERROR] {label}: No sequences found in {path}")
        sys.exit(1)
    total_len = sum(len(r.seq) for r in records)
    print(f"[OK] {label}: {len(records)} record(s), total {total_len:,} bp")
    for r in records:
        if len(r.seq) < min_len:
            print(f"  [WARN] Record {r.id} is suspiciously short ({len(r.seq)} bp)")
    return records

if __name__ == "__main__":
    genome = validate_fasta("data/raw/kpneumoniae_HS11286_genome.fasta",
                             "Genome", min_len=1_000_000)
    gene   = validate_fasta("data/raw/blaCTX-M-15_ref.fasta",
                             "blaCTX-M-15", min_len=700)
    print("\n[SUMMARY] Inputs look valid. Proceed to PAM scanning.")
