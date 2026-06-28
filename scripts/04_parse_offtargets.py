#!/usr/bin/env python3
"""
Parse Cas-OFFinder results and compute per-gRNA off-target summary.
Columns produced: grna_spacer, n_offtargets_0mm, n_offtargets_1mm,
                  n_offtargets_2mm, n_offtargets_3mm, total_offtargets
"""

import pandas as pd
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--casoff",  required=True, help="Cas-OFFinder results file")
    parser.add_argument("--grnas",   required=True, help="Filtered gRNA CSV")
    parser.add_argument("--output",  required=True, help="Output CSV with OT counts")
    args = parser.parse_args()

    # Load Cas-OFFinder results (tab-separated, no header)
    cols = ["spacer", "chromosome", "position", "direction",
            "off_target_seq", "mismatches"]
    try:
        ot_df = pd.read_csv(args.casoff, sep="\t", header=None, names=cols)
    except Exception as e:
        print(f"[WARN] Could not parse Cas-OFFinder output: {e}")
        ot_df = pd.DataFrame(columns=cols)

    # Pivot: count off-targets per spacer per mismatch level
    if not ot_df.empty:
        # !!! এই লাইনটি যোগ করুন: ২৩ বেস থেকে কেটে শুধু প্রথম ২০ বেস রাখুন !!!
        ot_df["spacer"] = ot_df["spacer"].str[:20] 
        
        summary = (
            ot_df.groupby(["spacer", "mismatches"])
                 .size()
                 .unstack(fill_value=0)
                 .reset_index()
        )
        # Rename columns
        mm_cols = {0: "ot_0mm", 1: "ot_1mm", 2: "ot_2mm", 3: "ot_3mm"}
        summary.rename(columns=mm_cols, inplace=True)
        for col in mm_cols.values():
            if col not in summary.columns:
                summary[col] = 0
        summary["total_offtargets"] = (
            summary[["ot_0mm","ot_1mm","ot_2mm","ot_3mm"]].sum(axis=1)
        )
    else:
        summary = pd.DataFrame(columns=["spacer","ot_0mm","ot_1mm",
                                         "ot_2mm","ot_3mm","total_offtargets"])

    # IMPORTANT: exclude on-target site (1 hit at 0mm is expected — the gene itself)
    # If ot_0mm == 1, the only match is the target gene → acceptable
    summary["ot_0mm_adjusted"] = (summary.get("ot_0mm", 0) - 1).clip(lower=0)

    # Merge with gRNA table
    grna_df = pd.read_csv(args.grnas)
    merged  = grna_df.merge(summary, left_on="spacer", right_on="spacer", how="left")
    for col in ["ot_0mm","ot_1mm","ot_2mm","ot_3mm","total_offtargets","ot_0mm_adjusted"]:
        merged[col] = merged[col].fillna(0).astype(int)

    merged.to_csv(args.output, index=False)
    print(f"[OK] Off-target summary merged: {args.output}")
    print(merged[["spacer","ot_0mm_adjusted","ot_1mm","ot_2mm","ot_3mm",
                   "total_offtargets"]].head(10).to_string(index=False))

if __name__ == "__main__":
    main()
