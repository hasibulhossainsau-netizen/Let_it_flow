#!/usr/bin/env python3
"""Compute composite efficiency scores and rank gRNA candidates.

Usage:
    python scripts/rank_candidates.py --grnas in.csv --output out.csv \
        [--config config/config.yaml]
"""

import argparse

import pandas as pd

import scoring_metrics as sm
from utils.config import load_config
from utils.dataframe import fill_missing_columns, require_columns
from utils.runtime import event, get_logger, init_script, run_main

log = get_logger("rank_candidates")


def rank_candidates(df: pd.DataFrame, scoring: dict) -> pd.DataFrame:
    """Score, assign verdicts to, and rank ``df``; returns a new DataFrame."""
    df = df.copy()
    require_columns(df, ["gc_content"], "gRNA table")
    fill_missing_columns(df, scoring["missing_column_defaults"])

    df["gc_score"] = df["gc_content"].apply(
        lambda v: sm.gc_score(v, scoring["gc_curve"]))
    df["seed_score"] = df["seed_gc"].apply(
        lambda v: sm.seed_gc_score(v, scoring["seed_full_score_gc"]))
    df["ot_score"] = df["total_offtargets"].apply(
        lambda v: sm.offtarget_penalty_score(
            v, scoring["off_target_scores"],
            scoring["off_target_fallback_score"]))

    df["efficiency_score"] = [
        sm.composite_score(g, s, m, o, scoring["weights"],
                           scoring["score_decimals"])
        for g, s, m, o in zip(df["gc_score"], df["seed_score"],
                              df["mfe_score"], df["ot_score"])
    ]
    df["verdict"] = [
        sm.verdict(sc, ot, ot0, scoring["verdict"])
        for sc, ot, ot0 in zip(df["efficiency_score"], df["total_offtargets"],
                               df["ot_0mm_adjusted"])
    ]

    # Fully deterministic order: best score first; ties broken by grna_id
    # (ascending), so equal-scoring guides never depend on input row order.
    # Merge sort ("stable") additionally keeps any remaining exact ties in
    # their input order instead of an unspecified quicksort order.
    require_columns(df, ["grna_id"], "gRNA table")
    df = df.sort_values(["efficiency_score", "grna_id"],
                        ascending=[False, True], kind="stable")
    df.insert(0, "rank", range(1, len(df) + 1))
    return df


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--grnas", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--config", default="config/config.yaml")
    args = parser.parse_args()

    cfg = load_config(args.config)
    init_script(log, args, cfg, {"grnas": args.grnas})
    df = rank_candidates(pd.read_csv(args.grnas), cfg["scoring"])
    df.to_csv(args.output, index=False)

    top = df[df["verdict"] == "TOP CANDIDATE"]
    event(log, "INFO", "final ranked table written", output=args.output,
          rows=len(df), top_candidates=len(top),
          verdicts={str(k): int(v) for k, v in df["verdict"].value_counts().items()})
    log.info("top candidates:\n%s",
             top[["rank", "spacer", "pam", "gc_content", "mfe_kcal",
                  "total_offtargets", "efficiency_score", "verdict"]]
             .to_string(index=False))


if __name__ == "__main__":
    run_main(main, "rank_candidates")
