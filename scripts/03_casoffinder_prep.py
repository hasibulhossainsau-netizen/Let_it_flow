#!/usr/bin/env python3
"""
Prepare Cas-OFFinder input file from filtered gRNA CSV.

Cas-OFFinder input format (tab-separated):
Line 1: /path/to/genome.fasta
Line 2: PAM pattern (e.g., NNNNNNNNNNNNNNNNNNNNNGG for 20bp + NGG)
Line 3+: <SPACER_SEQUENCE> <MISMATCH_COUNT>
"""

import pandas as pd
import argparse
import os

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input",   required=True, help="Filtered gRNA CSV")
    parser.add_argument("--genome",  required=True, help="Path to genome FASTA")
    parser.add_argument("--output",  required=True, help="Cas-OFFinder input .txt")
    parser.add_argument("--mm",      type=int, default=3, help="Max mismatches (default: 3)")
    args = parser.parse_args()

    df = pd.read_csv(args.input)
    df = df[df["pre_filter_pass"] == True].copy()

    pam_pattern = "N" * 20 + "NGG"   # 20bp spacer + NGG

    with open(args.output, "w") as f:
        f.write(os.path.abspath(args.genome) + "\n")
        f.write(pam_pattern + "\n")
        for _, row in df.iterrows():
            # Cas-OFFinder needs both spacer + PAM (total 23bp) to match the 23-char pattern
            spacer = row["spacer"].upper().replace("U", "T")
            pam = row["pam"].upper().replace("U", "T")
            full_target = spacer + pam
            f.write(f"{full_target}\t{args.mm}\n")

    print(f"[OK] Cas-OFFinder input written: {args.output}")
    print(f"     {len(df)} gRNA queries | max mismatches: {args.mm}")

if __name__ == "__main__":
    main()