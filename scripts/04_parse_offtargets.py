#!/usr/bin/env python3
"""
Parse Cas-OFFinder results and compute per-gRNA off-target summary.
Columns produced: grna_spacer, n_offtargets_0mm, n_offtargets_1mm,
                  n_offtargets_2mm, n_offtargets_3mm, total_offtargets
"""

import argparse

import pandas as pd

from utils.runtime import event, get_logger, init_script, run_main

log = get_logger("04_parse_offtargets")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--casoff",  required=True, help="Cas-OFFinder results file")
    parser.add_argument("--grnas",   required=True, help="Filtered gRNA CSV")
    parser.add_argument("--output",  required=True, help="Output CSV with OT counts")
    args = parser.parse_args()
    init_script(log, args, None, {"casoff": args.casoff, "grnas": args.grnas})

    # Load Cas-OFFinder results (tab-separated, no header)
    cols = ["spacer", "chromosome", "position", "direction",
            "off_target_seq", "mismatches"]
    try:
        ot_df = pd.read_csv(args.casoff, sep="\t", header=None, names=cols)
    except Exception as e:
        # An empty result file (zero hits) is legitimate; anything else is logged.
        event(log, "WARNING", "no Cas-OFFinder hits parsed (empty or unreadable)",
              error=str(e), path=args.casoff)
        ot_df = pd.DataFrame(columns=cols)

    # Pivot: count off-targets per spacer per mismatch level
    if not ot_df.empty:
        # Cas-OFFinder reports spacer+PAM (23 nt); keep only the 20-nt spacer.
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
    event(log, "INFO", "off-target summary merged", output=args.output,
          rows=len(merged), cas_offinder_hits=len(ot_df))
    log.info("first rows:\n%s",
             merged[["spacer","ot_0mm_adjusted","ot_1mm","ot_2mm","ot_3mm",
                     "total_offtargets"]].head(10).to_string(index=False))

if __name__ == "__main__":
    run_main(main, "04_parse_offtargets")
