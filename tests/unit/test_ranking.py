"""Ranking must be deterministic: independent of input row order."""
import pandas as pd
import pytest

from rank_candidates import rank_candidates


def make_table(n=12):
    """Many guides with IDENTICAL scores, so only the tie-breaker orders them."""
    return pd.DataFrame({
        "grna_id": [f"g_{'+-'[i % 2]}_{i * 37 % 101:03d}" for i in range(n)],
        "gc_content": [0.55] * n, "seed_gc": [0.6] * n, "mfe_score": [0.8] * n,
        "total_offtargets": [0] * n, "ot_0mm_adjusted": [0] * n,
        "blast_offtarget_hits": [0] * n,
    })


@pytest.fixture
def scoring(shipped_config):
    return shipped_config["scoring"]


def test_ties_are_broken_by_grna_id(scoring):
    ranked = rank_candidates(make_table(), scoring)
    assert ranked["efficiency_score"].nunique() == 1
    assert ranked["grna_id"].tolist() == sorted(ranked["grna_id"])


@pytest.mark.parametrize("seed", range(5))
def test_row_order_of_input_does_not_matter(scoring, seed):
    table = make_table()
    reference = rank_candidates(table, scoring)
    shuffled = rank_candidates(table.sample(frac=1, random_state=seed), scoring)
    pd.testing.assert_frame_equal(reference.reset_index(drop=True),
                                  shuffled.reset_index(drop=True))


def test_higher_score_wins_over_id(scoring):
    table = make_table(4)
    table.loc[table["grna_id"] == table["grna_id"].max(), "mfe_score"] = 1.0
    ranked = rank_candidates(table, scoring)
    assert ranked["grna_id"].iloc[0] == table["grna_id"].max()


def test_ranks_are_1_to_n_and_verdicts_assigned(scoring):
    ranked = rank_candidates(make_table(), scoring)
    assert ranked["rank"].tolist() == list(range(1, 13))
    assert ranked["verdict"].notna().all()


def test_exact_offtarget_is_rejected_and_ranked_by_score(scoring):
    table = make_table(3)
    table.loc[0, "ot_0mm_adjusted"] = 1
    ranked = rank_candidates(table, scoring)
    assert ranked.loc[ranked["grna_id"] == table.loc[0, "grna_id"], "verdict"].item() \
        == "REJECT — exact off-target"


def test_missing_grna_id_column_is_reported(scoring):
    with pytest.raises(ValueError, match="grna_id"):
        rank_candidates(make_table().drop(columns=["grna_id"]), scoring)


def test_input_is_not_mutated(scoring):
    table = make_table()
    before = table.copy()
    rank_candidates(table, scoring)
    pd.testing.assert_frame_equal(table, before)
