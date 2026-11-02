"""Tests for the post-step sanity checks (scripts/validate_outputs.py)."""
import numpy as np
import pandas as pd
import pytest

from conftest import load_script

vo = load_script("validate_outputs")


def candidates(n=4):
    return pd.DataFrame({
        "grna_id": [f"g{i}" for i in range(n)], "spacer": ["ACGT" * 5] * n,
        "pam": ["GGA"] * n, "strand": ["+"] * n, "position_nt": range(n),
        "gc_content": [0.5] * n, "pass_gc": [True] * n, "pass_homopol": [True] * n,
        "pre_filter_pass": [True] * n, "seed_gc": [0.5] * n, "good_seed": [True] * n,
    })


def errors_for(df, stage="candidates_raw", cfg=None):
    errors = []
    vo.check_table(df, stage, vo.CSV_STAGES[stage], cfg or {}, errors)
    return errors


def test_valid_table_has_no_errors():
    assert errors_for(candidates()) == []


def test_empty_table():
    assert errors_for(candidates(0))[0].startswith("table has a header but no rows")


def test_missing_column():
    assert "missing required column(s): ['seed_gc']" in errors_for(candidates().drop(columns=["seed_gc"]))[0]


def test_nan_in_key_column():
    df = candidates()
    df.loc[1, "gc_content"] = np.nan
    assert any("'gc_content' has 1 NaN" in e for e in errors_for(df))


def test_duplicate_ids():
    df = candidates()
    df.loc[1, "grna_id"] = "g0"
    assert any("duplicates" in e for e in errors_for(df))


def test_filtered_stage_rejects_failed_rows():
    df = candidates()
    df.loc[0, "pre_filter_pass"] = False
    assert any("pre_filter_pass" in e for e in errors_for(df, "candidates_filtered"))


def test_score_outside_unit_interval():
    df = candidates()
    df.loc[0, "gc_content"] = 1.5
    assert any("outside [0, 1]" in e for e in errors_for(df))


def test_rnafold_stage_flags_total_prediction_failure():
    df = candidates()
    for col in ("ot_0mm", "ot_1mm", "ot_2mm", "ot_3mm", "total_offtargets",
                "ot_0mm_adjusted", "blast_offtarget_hits"):
        df[col] = 0
    df["mfe_kcal"], df["mfe_score"], df["structure"] = -5.0, 0.8, np.nan
    assert any("every prediction failed" in e for e in errors_for(df, "rnafold"))


def test_ranked_stage_requires_rank_permutation(shipped_config):
    ranked = pd.DataFrame({"rank": [1, 1]})
    errors = []
    df = candidates(2)
    for col in vo.OFFTARGET_COLS + ["blast_offtarget_hits"]:
        df[col] = 0
    df["mfe_kcal"], df["mfe_score"] = -5.0, 0.8
    df["gc_score"] = df["seed_score"] = df["ot_score"] = df["efficiency_score"] = 0.9
    df["verdict"] = "CANDIDATE"
    df["rank"] = ranked["rank"]
    vo.check_table(df, "ranked", vo.CSV_STAGES["ranked"], shipped_config, errors)
    assert any("permutation" in e for e in errors)


def test_text_stage_blast_requires_hits_and_field_count(tmp_path, shipped_config):
    errors = []
    empty = tmp_path / "b.tsv"
    empty.write_text("")
    vo.check_text(str(empty), "blast_results", shipped_config, errors)
    assert any("no hits" in e for e in errors)
    errors.clear()
    bad = tmp_path / "b2.tsv"
    bad.write_text("a\tb\tc\n")
    vo.check_text(str(bad), "blast_results", shipped_config, errors)
    assert any("12 tab-separated" in e for e in errors)


def test_text_stage_casoffinder_allows_empty(tmp_path, shipped_config):
    errors = []
    path = tmp_path / "c.txt"
    path.write_text("")
    vo.check_text(str(path), "casoffinder_results", shipped_config, errors)
    assert errors == []


def test_fasta_queries_stage(tmp_path, shipped_config):
    errors = []
    path = tmp_path / "q.fa"
    path.write_text(">a\n>b\nACGT\n")
    vo.check_text(str(path), "blast_queries", shipped_config, errors)
    assert any("empty sequence" in e for e in errors)
