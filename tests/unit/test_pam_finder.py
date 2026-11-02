"""Guide extraction behaviour (scripts/02_pam_finder.py).

The default ``guide.pam_regex`` must implement SpCas9's NGG PAM: a 20-nt
spacer followed by one arbitrary base and then ``GG``. (It used to match
``GGN``; see REFACTORING_REPORT.md, section 4.)
"""
import copy

import pytest

from conftest import load_script
from utils.seq import reverse_complement

pam = load_script("02_pam_finder")
GOOD_SPACER = "ACGTTGCAGCTAGCATCGAT"           # GC 0.5, no homopolymer


@pytest.fixture
def guide(shipped_config):
    return shipped_config["guide"]


def extract(seq, guide, strand="+"):
    return pam.extract_grnas(seq, strand, "x", guide)


def test_default_regex_is_ngg(guide, repo_root):
    assert guide["pam_regex"] == "(?=(.{20}).GG)"
    schema = (repo_root / "config" / "schema.yaml").read_text()
    assert "'(?=(.{20}).GG)'" in schema, "schema default must match the shipped config"


@pytest.mark.parametrize("n", list("ACGT"))
def test_spacer_and_ngg_pam_are_split_correctly(guide, n):
    hits = extract(GOOD_SPACER + n + "GG" + "ACGTA", guide)
    assert [(h["spacer"], h["pam"]) for h in hits] == [(GOOD_SPACER, n + "GG")]


def test_spacer_does_not_swallow_the_n(guide):
    """Regression for the old GGN behaviour: A*20 + T + GG."""
    hits = extract("A" * 20 + "TGG" + "C" * 5, guide)
    assert [(h["spacer"], h["pam"]) for h in hits] == [("A" * 20, "TGG")]


def test_pam_always_ends_in_gg(guide):
    seq = "ACGTTGCAGCTAGCATCGATAGGCTTAGCCGGATTACGTAGGCTAGCTAGGATCGATCGGCATCGAGG"
    hits = extract(seq, guide) + extract(reverse_complement(seq), guide, "-")
    assert hits and all(h["pam"].endswith("GG") and len(h["pam"]) == 3 for h in hits)
    assert all(len(h["spacer"]) == 20 for h in hits)


def test_gg_directly_after_spacer_without_n_is_not_enough(guide):
    """Only 2 nt after the spacer, so there is no room for N+GG."""
    assert extract(GOOD_SPACER + "GG", guide) == []


def test_old_ggn_regex_still_usable_via_config(guide):
    legacy = copy.deepcopy(guide)
    legacy["pam_regex"] = "(?=(.{20})GG)"
    hits = extract("A" * 20 + "TGG" + "C" * 5, legacy)
    assert [(h["spacer"], h["pam"]) for h in hits] == [("A" * 19 + "T", "GGC")]


def test_schema_row_has_expected_columns(guide):
    hits = extract(GOOD_SPACER + "AGG" + "ACGTA", guide)
    assert set(hits[0]) == {
        "grna_id", "spacer", "pam", "strand", "position_nt", "gc_content",
        "pass_gc", "pass_homopol", "pre_filter_pass", "seed_gc", "good_seed"}


def test_grna_id_encodes_record_strand_and_position(guide):
    hits = pam.extract_grnas("A" * 20 + "TGG", "+", "rec", guide)
    assert hits[0]["grna_id"] == "rec_+_20"


def test_no_gg_no_candidates(guide):
    assert extract("ACT" * 30, guide) == []


def test_gc_filter_and_homopolymer_filter(guide):
    low = extract("A" * 20 + "TGG", guide)[0]
    assert not low["pass_gc"] and not low["pass_homopol"] and not low["pre_filter_pass"]
    assert extract(GOOD_SPACER + "TGG", guide)[0]["pre_filter_pass"]
