"""Unit tests for scripts/scoring_metrics.py (pure scoring functions)."""
import pytest

import scoring_metrics as sm


@pytest.fixture
def scoring(shipped_config):
    return shipped_config["scoring"]


class TestGcScore:
    @pytest.mark.parametrize("gc", [0.50, 0.55, 0.65])
    def test_optimal_window_scores_one(self, scoring, gc):
        assert sm.gc_score(gc, scoring["gc_curve"]) == 1.0

    def test_lower_ramp(self, scoring):
        # lower_base + (gc - lower_bound) * lower_slope = 0.70 + 0.05 * 3
        assert sm.gc_score(0.45, scoring["gc_curve"]) == pytest.approx(0.85)

    def test_lower_bound_gets_lower_base(self, scoring):
        assert sm.gc_score(0.40, scoring["gc_curve"]) == pytest.approx(0.70)

    def test_upper_ramp(self, scoring):
        # 1 - (gc - optimal_max) * upper_slope = 1 - 0.05 * 5
        assert sm.gc_score(0.70, scoring["gc_curve"]) == pytest.approx(0.75)

    @pytest.mark.parametrize("gc", [0.0, 0.20, 0.399, 0.76, 1.0])
    def test_out_of_range_gets_floor(self, scoring, gc):
        assert sm.gc_score(gc, scoring["gc_curve"]) == scoring["gc_curve"]["out_of_range"]

    def test_scores_stay_in_unit_interval_inside_bounds(self, scoring):
        for i in range(40, 76):
            assert 0.0 <= sm.gc_score(i / 100, scoring["gc_curve"]) <= 1.0


class TestSeedGcScore:
    def test_proportional_below_full_score(self):
        assert sm.seed_gc_score(0.25, 0.5) == pytest.approx(0.5)

    @pytest.mark.parametrize("seed", [0.5, 0.6, 1.0])
    def test_capped_at_one(self, seed):
        assert sm.seed_gc_score(seed, 0.5) == 1.0

    def test_zero(self):
        assert sm.seed_gc_score(0.0, 0.5) == 0.0


class TestOfftargetPenaltyScore:
    @pytest.mark.parametrize("hits, expected", [
        (0, 1.00), (1, 0.85), (2, 0.65), (3, 0.45), (5, 0.45), (6, 0.25),
        (10, 0.25), (11, 0.05), (500, 0.05),
    ])
    def test_default_tiers(self, scoring, hits, expected):
        assert sm.offtarget_penalty_score(
            hits, scoring["off_target_scores"],
            scoring["off_target_fallback_score"]) == expected

    def test_custom_tiers_and_fallback(self):
        tiers = [{"max_hits": 2, "score": 0.9}]
        assert sm.offtarget_penalty_score(2, tiers, 0.1) == 0.9
        assert sm.offtarget_penalty_score(3, tiers, 0.1) == 0.1


class TestMfeScore:
    def test_zero_mfe_is_best(self):
        assert sm.mfe_score(0.0, -30.0, 3) == 1.0

    def test_half_floor(self):
        assert sm.mfe_score(-15.0, -30.0, 3) == 0.5

    def test_clipped_at_floor(self):
        assert sm.mfe_score(-30.0, -30.0, 3) == 0.0
        assert sm.mfe_score(-55.0, -30.0, 3) == 0.0

    def test_rounding(self):
        assert sm.mfe_score(-8.3, -30.0, 3) == 0.723

    def test_monotonic(self):
        values = [sm.mfe_score(m, -30.0, 6) for m in (0, -5, -10, -20, -30)]
        assert values == sorted(values, reverse=True)


class TestCompositeScore:
    def test_all_ones_is_one(self, scoring):
        assert sm.composite_score(1, 1, 1, 1, scoring["weights"], 4) == 1.0

    def test_all_zeros_is_zero(self, scoring):
        assert sm.composite_score(0, 0, 0, 0, scoring["weights"], 4) == 0.0

    def test_weighted_sum(self, scoring):
        # 0.25*0.8 + 0.20*0.5 + 0.25*0.4 + 0.30*1.0 = 0.70
        assert sm.composite_score(0.8, 0.5, 0.4, 1.0,
                                  scoring["weights"], 4) == pytest.approx(0.70)

    def test_rounds_to_requested_decimals(self):
        weights = {"gc": 1 / 3, "seed_gc": 1 / 3, "mfe": 1 / 3, "off_target": 0}
        assert sm.composite_score(1, 0, 0, 0, weights, 2) == 0.33


class TestVerdict:
    def test_exact_offtarget_always_rejects(self, scoring):
        assert sm.verdict(1.0, 0, 1, scoring["verdict"]) == "REJECT — exact off-target"

    def test_high_offtarget_burden_rejects(self, scoring):
        assert sm.verdict(1.0, 11, 0, scoring["verdict"]) == "REJECT — high off-target burden"

    @pytest.mark.parametrize("score, expected", [
        (0.95, "TOP CANDIDATE"), (0.75, "TOP CANDIDATE"),
        (0.74, "CANDIDATE"), (0.55, "CANDIDATE"),
        (0.54, "MARGINAL"), (0.35, "MARGINAL"),
        (0.34, "REJECT — low efficiency"), (0.0, "REJECT — low efficiency"),
    ])
    def test_score_thresholds(self, scoring, score, expected):
        assert sm.verdict(score, 0, 0, scoring["verdict"]) == expected
