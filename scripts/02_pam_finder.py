#!/usr/bin/env python3
"""
PAM Finder for SpCas9 (NGG PAM).
Scans forward and reverse complement strands of the target gene.
Outputs a CSV of candidate gRNA spacers.

Usage:
    python scripts/02_pam_finder.py \
        --target data/raw/blaCTX-M-15_ref.fasta \
        --output data/processed/candidate_grnas_raw.csv
"""

import re
import argparse
import csv
from Bio import SeqIO
from Bio.Seq import Seq

# ── Constants ────────────────────────────────────────────────────────────────
SPACER_LEN  = 20          # bp upstream of PAM
PAM_PATTERN = r'(?=(.{20})GG)'   # NGG: any N + GG, with 20bp spacer
GC_MIN      = 0.40        # 40% minimum GC
GC_MAX      = 0.75        # 75% maximum GC
MAX_HOMOPOL = 5           # reject if any nucleotide runs >= 5 in a row
SEED_LEN    = 12          # length of PAM-proximal seed region

# ── Helper functions ──────────────────────────────────────────────────────────

def gc_content(seq: str) -> float:
    seq = seq.upper()
    return (seq.count('G') + seq.count('C')) / len(seq)

def seed_gc(spacer: str, seed_len: int = SEED_LEN) -> float:
    """GC content of the PAM-proximal seed region (last seed_len nt of spacer)."""
    seed = spacer[-seed_len:]
    return gc_content(seed)

def has_homopolymer(seq: str, max_run: int = MAX_HOMOPOL) -> bool:
    """Return True if the sequence contains a homopolymer run >= max_run."""
    return bool(re.search(r'([ACGT])\1{' + str(max_run - 1) + r',}', seq.upper()))

def extract_grnas(seq: str, strand: str, seq_id: str) -> list[dict]:
    """Extract all candidate gRNAs from a sequence string."""
    candidates = []
    for m in re.finditer(PAM_PATTERN, seq.upper()):
        spacer   = m.group(1)               # 20 bp spacer
        pam_pos  = m.start() + SPACER_LEN   # position of N in NGG
        pam      = seq[pam_pos:pam_pos + 3].upper()  # NGG
        gc       = gc_content(spacer)
        homo     = has_homopolymer(spacer)
        pass_gc  = GC_MIN <= gc <= GC_MAX
        pass_hom = not homo

        # Seed region GC
        sgc = seed_gc(spacer)
        good_seed = sgc >= 0.40

        candidates.append({
            "grna_id"      : f"{seq_id}_{strand}_{pam_pos}",
            "spacer"       : spacer,
            "pam"          : pam,
            "strand"       : strand,
            "position_nt"  : pam_pos,
            "gc_content"   : round(gc, 3),
            "pass_gc"      : pass_gc,
            "pass_homopol" : pass_hom,
            "pre_filter_pass": pass_gc and pass_hom,
            "seed_gc"      : round(sgc, 3),
            "good_seed"    : good_seed,
        })
    return candidates

def reverse_complement(seq: str) -> str:
    return str(Seq(seq).reverse_complement())

# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="PAM Finder — SpCas9 NGG")
    parser.add_argument("--target", required=True, help="Target gene FASTA")
    parser.add_argument("--output", required=True, help="Output CSV path")
    args = parser.parse_args()

    all_candidates = []

    for record in SeqIO.parse(args.target, "fasta"):
        fwd_seq = str(record.seq)
        rev_seq = reverse_complement(fwd_seq)
        seq_id  = record.id.split(".")[0]

        print(f"[INFO] Scanning {seq_id} ({len(fwd_seq)} bp) on both strands...")

        fwd_hits = extract_grnas(fwd_seq, "+", seq_id)
        rev_hits = extract_grnas(rev_seq, "-", seq_id)

        print(f"  Forward strand: {len(fwd_hits)} raw candidates")
        print(f"  Reverse strand: {len(rev_hits)} raw candidates")

        all_candidates.extend(fwd_hits)
        all_candidates.extend(rev_hits)

    # Write output
    if not all_candidates:
        print("[WARN] No gRNA candidates found. Check input sequence.")
        return

    fieldnames = list(all_candidates[0].keys())
    with open(args.output, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(all_candidates)

    passed = sum(1 for c in all_candidates if c["pre_filter_pass"])
    print(f"\n[SUMMARY] Total candidates: {len(all_candidates)}")
    print(f"[SUMMARY] Passed GC + homopolymer filter: {passed}")
    print(f"[OUTPUT] Written to: {args.output}")

if __name__ == "__main__":
    main()