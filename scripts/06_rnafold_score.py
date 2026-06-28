#!/usr/bin/env python3
"""
Run RNAfold on each gRNA spacer + truncated scaffold and extract MFE.
Less negative MFE = less structure = better accessibility = higher efficiency.

Usage:
    python scripts/06_rnafold_score.py \
        --grnas  data/processed/grnas_with_offtargets.csv \
        --output data/processed/grnas_rnafold.csv
"""

import subprocess
import re
import pandas as pd
import argparse
import tempfile
import os

# Partial sgRNA scaffold (first 20 nt; full structure modelling)
SCAFFOLD = "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCG"

def run_rnafold(spacer: str) -> tuple[float, str]:
    """
    Fold spacer+scaffold with RNAfold.
    Returns (MFE float, dot-bracket structure string).
    """
    rna_seq = spacer.replace("T", "U") + SCAFFOLD.replace("T", "U")

    with tempfile.NamedTemporaryFile(mode='w', suffix='.fa',
                                     delete=False) as tmp:
        tmp.write(f">grna\n{rna_seq}\n")
        tmp_path = tmp.name

    try:
        result = subprocess.run(
            ["RNAfold", "--noPS", tmp_path],
            capture_output=True, text=True, timeout=30
        )
        output = result.stdout.strip().split("\n")
        # Last line: structure  (MFE)
        last = output[-1] if output else ""
        match = re.search(r'\((-?\d+\.\d+)\)', last)
        mfe   = float(match.group(1)) if match else 0.0
        struct = last.split()[0] if last else ""
        return mfe, struct
    except Exception as e:
        print(f"[WARN] RNAfold failed for {spacer[:8]}...: {e}")
        return 0.0, ""
    finally:
        os.unlink(tmp_path)

def mfe_score(mfe: float) -> float:
    """
    Convert MFE to efficiency score contribution (0–1 scale).
    MFE closer to 0 = less structure = better. 
    Typical range: -30 to 0 kcal/mol.
    """
    clipped = max(mfe, -30.0)
    return round(1.0 - (abs(clipped) / 30.0), 3)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--grnas",  required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    df = pd.read_csv(args.grnas)
    df = df[df["pre_filter_pass"] == True].copy()

    print(f"[INFO] Running RNAfold on {len(df)} gRNAs...")
    mfe_vals, struct_vals, mfe_scores = [], [], []

    for i, (_, row) in enumerate(df.iterrows()):
        mfe, struct = run_rnafold(row["spacer"])
        mfe_vals.append(mfe)
        struct_vals.append(struct[:40] if struct else "")
        mfe_scores.append(mfe_score(mfe))
        if (i + 1) % 10 == 0:
            print(f"  Processed {i+1}/{len(df)}...")

    df["mfe_kcal"]   = mfe_vals
    df["structure"]  = struct_vals
    df["mfe_score"]  = mfe_scores   # 0–1 scale; higher = better

    df.to_csv(args.output, index=False)
    print(f"[OK] RNAfold scores written to: {args.output}")
    print(df[["spacer","mfe_kcal","mfe_score"]].head(10).to_string(index=False))

if __name__ == "__main__":
    main()
