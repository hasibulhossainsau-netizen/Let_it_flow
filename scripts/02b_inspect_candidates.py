#!/usr/bin/env python3
"""Keep only gRNA candidates that passed the pre-filter.

Usage:
    python scripts/02b_inspect_candidates.py \
        --input data/processed/candidate_grnas_raw.csv \
        --output data/processed/candidate_grnas_filtered.csv
"""
import argparse

import pandas as pd

from utils.dataframe import require_columns
from utils.runtime import event, get_logger, init_script, run_main

log = get_logger("02b_inspect_candidates")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", default="data/processed/candidate_grnas_raw.csv")
    parser.add_argument("--output",
                        default="data/processed/candidate_grnas_filtered.csv")
    args = parser.parse_args()

    init_script(log, args, None, {"input": args.input})
    df = pd.read_csv(args.input)
    require_columns(df, ["pre_filter_pass", "gc_content"], args.input)
    event(log, "INFO", "candidate table loaded", rows=len(df), cols=len(df.columns),
          pre_filter_counts={str(k): int(v) for k, v in
                             df["pre_filter_pass"].value_counts().items()})
    log.info("gc_content of passing candidates:\n%s",
             df[df["pre_filter_pass"]]["gc_content"].describe().to_string())

    df_pass = df[df["pre_filter_pass"]].copy()
    df_pass.to_csv(args.output, index=False)
    event(log, "INFO", "filtered candidates saved", kept=len(df_pass),
          output=args.output)


if __name__ == "__main__":
    run_main(main, "02b_inspect_candidates")
