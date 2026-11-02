#!/usr/bin/env python3
"""
Run RNAfold on each gRNA spacer + truncated scaffold and extract MFE.
Less negative MFE = less structure = better accessibility = higher efficiency.

Usage:
    python scripts/06_rnafold_score.py \
        --grnas  data/processed/grnas_with_offtargets.csv \
        --output data/processed/grnas_rnafold.csv
"""

import argparse

import pandas as pd

# ``utils`` must be imported first: it puts the repo root on sys.path so that
# ``interfaces`` (a top-level package) can be found.
from utils.config import load_config
from utils.dataframe import require_columns
from utils.runtime import event, get_logger, init_script, run_main

from interfaces import get_structure_predictor  # noqa: E402
from scoring_metrics import mfe_score  # noqa: E402


log = get_logger("06_rnafold_score")


def fold_spacer(predictor, spacer: str, scaffold: str):
    """Fold spacer+scaffold; returns (MFE, structure) or (0.0, "") on failure."""
    rna_seq = spacer.replace("T", "U") + scaffold.replace("T", "U")
    try:
        return predictor.predict(rna_seq)
    except RuntimeError as err:
        event(log, "WARNING", "structure prediction failed",
              spacer_prefix=spacer[:8], error=str(err))
        return 0.0, ""


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--grnas",  required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--config", default="config/config.yaml")
    args = parser.parse_args()

    cfg = load_config(args.config)
    struct_cfg = cfg["structure"]
    predictor = get_structure_predictor(cfg)
    init_script(log, args, cfg, {"grnas": args.grnas},
                tools={predictor.name: predictor})

    df = pd.read_csv(args.grnas)
    require_columns(df, ["spacer", "pre_filter_pass"], args.grnas)
    df = df[df["pre_filter_pass"] == True].copy()

    event(log, "INFO", "running structure prediction", n_grnas=len(df),
          predictor=predictor.name)
    mfe_vals, struct_vals, mfe_scores = [], [], []

    for i, (_, row) in enumerate(df.iterrows()):
        mfe, struct = fold_spacer(predictor, row["spacer"],
                                  struct_cfg["scaffold"])
        mfe_vals.append(mfe)
        struct_vals.append(struct[:40] if struct else "")
        mfe_scores.append(mfe_score(mfe, struct_cfg["mfe_floor"],
                                    struct_cfg["mfe_score_decimals"]))
        if (i + 1) % 10 == 0:
            event(log, "INFO", "progress", done=i + 1, total=len(df))

    df["mfe_kcal"]   = mfe_vals
    df["structure"]  = struct_vals
    df["mfe_score"]  = mfe_scores   # 0–1 scale; higher = better

    df.to_csv(args.output, index=False)
    n_failed = sum(1 for s in struct_vals if not s)
    event(log, "WARNING" if n_failed else "INFO", "RNAfold scores written",
          output=args.output, rows=len(df), failed_predictions=n_failed)
    log.info("first rows:\n%s",
             df[["spacer","mfe_kcal","mfe_score"]].head(10).to_string(index=False))

if __name__ == "__main__":
    run_main(main, "06_rnafold_score")
