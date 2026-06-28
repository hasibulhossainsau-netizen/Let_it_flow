#!/usr/bin/env python3
"""
Composite efficiency scorer.
Combines GC content, seed GC, MFE score, and off-target penalty.
Produces the final ranked gRNA table.
"""

import pandas as pd
import argparse

# ── Scoring weights (sum to 1.0) ──────────────────────────────────────────────
W_GC_SCORE   = 0.25   # Optimal GC (40–65%  → score 1.0; outside → lower)
W_SEED_GC    = 0.20   # Seed region GC ≥ 40%
W_MFE        = 0.25   # RNAfold MFE score
W_OT_PENALTY = 0.30   # Off-target penalty (lower OTs → higher score)

def gc_score(gc: float) -> float:
    """Score GC content. Optimal = 0.50–0.65 → 1.0"""
    if 0.50 <= gc <= 0.65:
        return 1.0
    elif 0.40 <= gc < 0.50:
        return 0.7 + (gc - 0.40) * 3.0
    elif 0.65 < gc <= 0.75:
        return 1.0 - (gc - 0.65) * 5.0
    else:
        return 0.3

def seed_gc_score(seed_gc: float) -> float:
    return min(1.0, seed_gc / 0.5)

def offtarget_penalty_score(total_ot: int) -> float:
    """Returns 1.0 for 0 off-targets, decreasing with more."""
    if total_ot == 0:   return 1.00
    elif total_ot == 1: return 0.85
    elif total_ot == 2: return 0.65
    elif total_ot <= 5: return 0.45
    elif total_ot <= 10: return 0.25
    else:               return 0.05

def verdict(score: float, total_ot: int, ot_0mm_adj: int) -> str:
    """Final selection verdict."""
    if ot_0mm_adj > 0:
        return "REJECT — exact off-target"
    if total_ot > 10:
        return "REJECT — high off-target burden"
    if score >= 0.75:
        return "TOP CANDIDATE"
    elif score >= 0.55:
        return "CANDIDATE"
    elif score >= 0.35:
        return "MARGINAL"
    else:
        return "REJECT — low efficiency"

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--grnas",  required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    df = pd.read_csv(args.grnas)

    # Ensure required columns exist with defaults
    for col, default in [("seed_gc", 0.5), ("mfe_score", 0.5),
                          ("total_offtargets", 0), ("ot_0mm_adjusted", 0),
                          ("blast_offtarget_hits", 0)]:
        if col not in df.columns:
            df[col] = default

    df["gc_score"]   = df["gc_content"].apply(gc_score)
    df["seed_score"] = df["seed_gc"].apply(seed_gc_score)
    df["ot_score"]   = df["total_offtargets"].apply(offtarget_penalty_score)

    df["efficiency_score"] = (
        W_GC_SCORE   * df["gc_score"]   +
        W_SEED_GC    * df["seed_score"] +
        W_MFE        * df["mfe_score"]  +
        W_OT_PENALTY * df["ot_score"]
    ).round(4)

    df["verdict"] = df.apply(
        lambda r: verdict(r["efficiency_score"],
                          r["total_offtargets"],
                          r["ot_0mm_adjusted"]), axis=1
    )

    # Rank by efficiency score descending
    df.sort_values("efficiency_score", ascending=False, inplace=True)
    df.insert(0, "rank", range(1, len(df) + 1))

    df.to_csv(args.output, index=False)
    print(f"\n[OK] Final ranked table written: {args.output}")
    top = df[df["verdict"] == "TOP CANDIDATE"]
    print(f"[SUMMARY] Top candidates: {len(top)}")
    print(top[["rank","spacer","pam","gc_content","mfe_kcal",
               "total_offtargets","efficiency_score","verdict"]].to_string(index=False))

if __name__ == "__main__":
    main()
