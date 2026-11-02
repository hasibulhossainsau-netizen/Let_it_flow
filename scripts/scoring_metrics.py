"""Pure scoring functions for gRNA candidates.

Every function takes its tunable parameters explicitly (sourced from the
``scoring`` section of config/config.yaml) and has no I/O or global state.
"""

import numpy as np


def gc_score(gc: float, curve: dict) -> float:
    """Score spacer GC content on a 0-1 scale.

    ``curve`` holds optimal_min/optimal_max (score 1.0), lower_bound,
    upper_bound, lower_base, lower_slope, upper_slope and out_of_range.
    """
    if curve["optimal_min"] <= gc <= curve["optimal_max"]:
        return 1.0
    if curve["lower_bound"] <= gc < curve["optimal_min"]:
        return curve["lower_base"] + (gc - curve["lower_bound"]) * curve["lower_slope"]
    if curve["optimal_max"] < gc <= curve["upper_bound"]:
        return 1.0 - (gc - curve["optimal_max"]) * curve["upper_slope"]
    return curve["out_of_range"]


def seed_gc_score(seed_gc: float, full_score_gc: float) -> float:
    """Score seed-region GC; reaches 1.0 at ``full_score_gc``."""
    return min(1.0, seed_gc / full_score_gc)


def offtarget_penalty_score(total_ot: int, tiers: list, fallback: float) -> float:
    """Return the score of the first tier with ``max_hits >= total_ot``.

    ``tiers`` is a list of ``{"max_hits": int, "score": float}`` dicts in
    ascending ``max_hits`` order; ``fallback`` applies beyond the last tier.
    """
    for tier in tiers:
        if total_ot <= tier["max_hits"]:
            return tier["score"]
    return fallback


def mfe_score(mfe: float, floor: float, decimals: int) -> float:
    """Convert an RNAfold MFE (kcal/mol) to a 0-1 score (closer to 0 = better)."""
    clipped = max(mfe, floor)
    return round(1.0 - (abs(clipped) / abs(floor)), decimals)


def composite_score(gc_s: float, seed_s: float, mfe_s: float, ot_s: float,
                    weights: dict, decimals: int) -> float:
    """Weighted sum of the four component scores.

    Rounds with ``numpy.round`` (not built-in ``round``) to stay bit-identical
    to the original vectorised pandas computation on exact decimal ties.
    """
    total = (weights["gc"] * gc_s
             + weights["seed_gc"] * seed_s
             + weights["mfe"] * mfe_s
             + weights["off_target"] * ot_s)
    return float(np.round(total, decimals))


def verdict(score: float, total_ot: int, ot_0mm_adj: int, rules: dict) -> str:
    """Final selection verdict from the composite score and off-target counts."""
    if ot_0mm_adj > 0:
        return "REJECT — exact off-target"
    if total_ot > rules["reject_total_offtargets_above"]:
        return "REJECT — high off-target burden"
    if score >= rules["top_candidate_min"]:
        return "TOP CANDIDATE"
    if score >= rules["candidate_min"]:
        return "CANDIDATE"
    if score >= rules["marginal_min"]:
        return "MARGINAL"
    return "REJECT — low efficiency"
