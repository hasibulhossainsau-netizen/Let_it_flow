"""Tests for the RNAfold output parser, including the padded-MFE regression.

RNAfold prints the energy in a fixed-width field: "(-12.30)" but "( -8.30)" and
"(  0.00)". The old regex ``\\((-?\\d+\\.\\d+)\\)`` rejected every value in
(-10, 0], so those guides silently received MFE 0.0 (the *best* score).
"""
import os
import shutil
import stat

import pytest

from interfaces.rnafold import RNAfoldPredictor, parse_rnafold_output

STRUCT = ".....(.(((((.......((((((((.((((....))))...))))))))..))).)).)."


class TestParseRnafoldOutput:
    @pytest.mark.parametrize("field, expected", [
        ("(-12.30)", -12.30),
        ("(-20.00)", -20.00),
        ("( -8.30)", -8.30),      # regression: padded single-digit energy
        ("( -0.30)", -0.30),
        ("(  0.00)", 0.0),        # regression: zero is padded too
        ("(-100.50)", -100.50),
    ])
    def test_energy_field_formats(self, field, expected):
        assert parse_rnafold_output(f"{STRUCT} {field}\n") == (expected, STRUCT)

    def test_returns_structure_token(self):
        mfe, structure = parse_rnafold_output(f"((..))  (-1.20)")
        assert (mfe, structure) == (-1.20, "((..))")

    def test_full_stdout_with_header_and_sequence(self):
        out = f">seq\nGGGAAACCC\n{STRUCT} ( -8.30)\n"
        assert parse_rnafold_output(out) == (-8.30, STRUCT)

    def test_trailing_blank_lines_and_spaces(self):
        assert parse_rnafold_output(f"{STRUCT} ( -8.30)  \n\n\n")[0] == -8.30

    @pytest.mark.parametrize("bad", ["", "\n\n", "no energy here", ".... (abc)",
                                     "....  -8.30", ">seq\nACGU"])
    def test_unparseable_raises_value_error(self, bad):
        with pytest.raises(ValueError, match="unparseable"):
            parse_rnafold_output(bad)

    def test_energy_must_be_on_last_line(self):
        with pytest.raises(ValueError):
            parse_rnafold_output(f"{STRUCT} (-5.00)\nsomething else")


def _stub(tmp_path, body, exit_code=0):
    """Executable that stands in for RNAfold (args are ignored)."""
    script = tmp_path / "RNAfold"
    script.write_text(f"#!/bin/sh\ncat <<'EOF'\n{body}\nEOF\nexit {exit_code}\n")
    script.chmod(script.stat().st_mode | stat.S_IEXEC)
    return str(script)


class TestPredictorWithStubBinary:
    def test_padded_energy_is_parsed_end_to_end(self, tmp_path):
        exe = _stub(tmp_path, f">seq\nACGU\n{STRUCT} ( -8.30)")
        assert RNAfoldPredictor(executable=exe).predict("ACGU") == (-8.30, STRUCT)

    def test_zero_energy_is_a_real_zero_not_a_failure(self, tmp_path):
        exe = _stub(tmp_path, ">seq\nACGU\n.... (  0.00)")
        assert RNAfoldPredictor(executable=exe).predict("ACGU") == (0.0, "....")

    def test_nonzero_exit_raises(self, tmp_path):
        exe = _stub(tmp_path, "boom", exit_code=2)
        with pytest.raises(RuntimeError, match="exit 2"):
            RNAfoldPredictor(executable=exe).predict("ACGU")

    def test_garbage_output_raises(self, tmp_path):
        exe = _stub(tmp_path, "nothing useful")
        with pytest.raises(RuntimeError, match="exit 0.*unparseable"):
            RNAfoldPredictor(executable=exe).predict("ACGU")

    def test_missing_binary_raises_runtime_error(self, tmp_path):
        with pytest.raises(RuntimeError, match="could not be run"):
            RNAfoldPredictor(executable=str(tmp_path / "absent")).predict("ACGU")


@pytest.mark.skipif(shutil.which("RNAfold") is None, reason="RNAfold not installed")
class TestRealRnafold:
    def test_weakly_folded_spacer_is_not_scored_as_zero(self):
        # A known case from the example data whose MFE (-8.3) used to fail.
        spacer = "GGCACCACCAACGATATCGC".replace("T", "U")
        scaffold = "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCG".replace("T", "U")
        mfe, structure = RNAfoldPredictor().predict(spacer + scaffold)
        assert -30 < mfe < 0
        assert len(structure) == len(spacer + scaffold)

    def test_unstructured_sequence_returns_zero_mfe(self):
        mfe, structure = RNAfoldPredictor().predict("AAAAAAAAAAAAAAAAAAAA")
        assert mfe == 0.0 and set(structure) == {"."}
