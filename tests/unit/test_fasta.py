"""Unit tests for FASTA parsing (utils.fasta) and input validation (01_validate_inputs)."""
import io

import pytest

from conftest import load_script
from utils.fasta import parse_fasta, read_fasta, write_fasta
from utils.seq import gc_content, has_homopolymer, reverse_complement

validate = load_script("01_validate_inputs")


class TestParseFasta:
    def test_multi_record_and_wrapped_lines(self):
        text = ">a first record\nACGT\nACGT\n\n>b\nTTTT\n"
        assert list(parse_fasta(io.StringIO(text))) == [("a", "ACGTACGT"), ("b", "TTTT")]

    def test_id_is_first_header_token(self):
        assert list(parse_fasta(io.StringIO(">NC_1.1 desc words\nAC\n"))) == [("NC_1.1", "AC")]

    def test_sequence_before_header_raises(self):
        with pytest.raises(ValueError, match="before first"):
            list(parse_fasta(io.StringIO("ACGT\n>a\nAC\n")))

    def test_empty_input_yields_nothing(self):
        assert list(parse_fasta(io.StringIO(""))) == []

    def test_write_read_roundtrip(self, tmp_path):
        path = tmp_path / "x.fa"
        records = [("g1", "ACGT" * 10), ("g2", "TTGCA")]
        write_fasta(path, records)
        assert read_fasta(path) == records


class TestSeqHelpers:
    def test_reverse_complement(self):
        assert reverse_complement("AACGT") == "ACGTT"

    def test_reverse_complement_is_involution(self):
        assert reverse_complement(reverse_complement("ACGTNRY")) == "ACGTNRY"

    def test_gc_content(self):
        assert gc_content("GGCC") == 1.0
        assert gc_content("acgt") == 0.5

    def test_gc_content_empty_raises(self):
        with pytest.raises(ValueError):
            gc_content("")

    def test_homopolymer(self):
        assert has_homopolymer("ACGAAAAATC", 5)
        assert not has_homopolymer("ACGAAAATC", 5)


def _write(tmp_path, text, name="in.fa"):
    path = tmp_path / name
    path.write_text(text)
    return str(path)


def _gene_problems(tmp_path, text, min_len=20):
    _, problems = validate.validate_gene(_write(tmp_path, text), min_len)
    return problems


class TestValidateGene:
    def test_valid_passes(self, tmp_path):
        assert _gene_problems(tmp_path, ">g\n" + "ACGT" * 10 + "\n").total == 0

    def test_lowercase_is_valid(self, tmp_path):
        assert _gene_problems(tmp_path, ">g\n" + "acgt" * 10 + "\n").total == 0

    def test_too_short(self, tmp_path):
        problems = _gene_problems(tmp_path, ">g\nACGT\n")
        assert problems.total == 1 and "below the minimum" in problems.items[0]

    def test_no_header(self, tmp_path):
        problems = _gene_problems(tmp_path, "ACGTACGTACGTACGTACGTACGT\n")
        assert any("before the first '>' header" in m for m in problems.items)

    def test_invalid_characters_report_position_and_hint(self, tmp_path):
        problems = _gene_problems(tmp_path, ">g\n" + "ACGT" * 5 + "U1N" + "ACGT" * 5 + "\n")
        text = " | ".join(problems.items)
        assert "col 21" in text and "RNA 'U'" in text and "digits" in text

    def test_iupac_not_allowed_in_gene(self, tmp_path):
        assert _gene_problems(tmp_path, ">g\n" + "ACGT" * 5 + "R" + "ACGT" * 5 + "\n").total >= 1

    def test_duplicate_id(self, tmp_path):
        seq = "ACGT" * 10
        problems = _gene_problems(tmp_path, f">a\n{seq}\n>a\n{seq}\n")
        assert any("duplicate record ID 'a'" in m for m in problems.items)

    def test_empty_record(self, tmp_path):
        problems = _gene_problems(tmp_path, ">a\n>b\n" + "ACGT" * 10 + "\n")
        assert any("record 'a' has no sequence" in m for m in problems.items)

    def test_header_without_id(self, tmp_path):
        problems = _gene_problems(tmp_path, ">\n" + "ACGT" * 10 + "\n")
        assert any("no ID" in m for m in problems.items)

    def test_whitespace_inside_sequence(self, tmp_path):
        problems = _gene_problems(tmp_path, ">g\nACGT ACGT ACGT ACGT ACGT\n")
        assert any("whitespace inside" in m for m in problems.items)

    def test_binary_file(self, tmp_path):
        path = tmp_path / "bin.fa"
        path.write_bytes(bytes(range(128, 256)))
        _, problems = validate.validate_gene(str(path), 20)
        assert any("not a UTF-8 text file" in m for m in problems.items)

    def test_empty_file(self, tmp_path):
        assert _gene_problems(tmp_path, "").total == 1

    def test_missing_file(self, tmp_path):
        _, problems = validate.validate_gene(str(tmp_path / "nope.fa"), 20)
        assert any("cannot open file" in m for m in problems.items)

    def test_problem_list_is_capped_but_count_is_exact(self, tmp_path):
        text = ">g\n" + "".join("ACGT" * 5 + "1\n" for _ in range(30))
        problems = _gene_problems(tmp_path, text)
        assert len(problems.items) == validate.MAX_REPORTED and problems.total == 30


class TestValidateGenome:
    def test_iupac_allowed_below_threshold(self, tmp_path):
        seq = "ACGT" * 250 + "N"
        _, problems = validate.validate_genome(_write(tmp_path, f">c\n{seq}\n"), 1000, 0.01)
        assert problems.total == 0

    def test_too_many_ambiguous(self, tmp_path):
        seq = "N" * 100 + "ACGT" * 225
        _, problems = validate.validate_genome(_write(tmp_path, f">c\n{seq}\n"), 1000, 0.01)
        assert any("ambiguity" in m for m in problems.items)

    def test_too_short_genome(self, tmp_path):
        _, problems = validate.validate_genome(_write(tmp_path, ">c\nACGT\n"), 1000, 0.01)
        assert any("below the minimum" in m for m in problems.items)

    def test_total_length_across_records(self, tmp_path):
        text = ">c1\n" + "ACGT" * 150 + "\n>c2\n" + "ACGT" * 150 + "\n"
        _, problems = validate.validate_genome(_write(tmp_path, text), 1000, 0.01)
        assert problems.total == 0

    def test_non_iupac_char_rejected(self, tmp_path):
        seq = "ACGT" * 250 + "X"
        _, problems = validate.validate_genome(_write(tmp_path, f">c\n{seq}\n"), 1000, 0.01)
        assert any("invalid character 'X'" in m for m in problems.items)
