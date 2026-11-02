"""Unit tests for config loading/validation (utils.config + config/schema.yaml)."""
import copy

import pytest
import yaml

from utils.config import apply_schema_defaults, load_config, validate_config

MINIMAL = {"genome_fasta": "g.fa", "target_gene_fasta": "t.fa"}


@pytest.fixture
def schema_path(repo_root):
    return str(repo_root / "config" / "schema.yaml")


def check(cfg, schema_path):
    return validate_config(copy.deepcopy(cfg), schema_path)


class TestDefaults:
    def test_minimal_config_gets_every_section(self, schema_path):
        cfg = check(MINIMAL, schema_path)
        for section in ("tools", "validation", "guide", "blast", "structure",
                        "scoring", "export"):
            assert section in cfg
        assert cfg["max_mismatches"] == 3 and cfg["seed"] == 42
        assert cfg["checksums_file"] == "data/raw/checksums.md5"

    def test_partial_section_keeps_user_value_and_fills_rest(self, schema_path):
        cfg = check({**MINIMAL, "guide": {"gc_min": 0.3}}, schema_path)
        assert cfg["guide"]["gc_min"] == 0.3 and cfg["guide"]["gc_max"] == 0.75

    def test_defaults_are_not_shared_between_loads(self, schema_path):
        a, b = check(MINIMAL, schema_path), check(MINIMAL, schema_path)
        a["scoring"]["off_target_scores"].append({"max_hits": 99, "score": 0})
        assert len(b["scoring"]["off_target_scores"]) == 5

    def test_apply_schema_defaults_is_recursive(self):
        schema = {"properties": {"a": {"type": "object", "properties": {
            "b": {"default": 1}}, "default": {}}}}
        assert apply_schema_defaults(schema, {}) == {"a": {"b": 1}}


class TestShippedConfig:
    def test_shipped_config_is_valid(self, shipped_config):
        assert shipped_config["scoring"]["weights"].keys() == {"gc", "seed_gc", "mfe", "off_target"}

    def test_weights_sum_to_one(self, shipped_config):
        assert sum(shipped_config["scoring"]["weights"].values()) == pytest.approx(1.0)

    def test_toy_config_is_valid(self, repo_root):
        cfg = load_config(str(repo_root / "tests/data/toy/config.yaml"),
                          str(repo_root / "config/schema.yaml"))
        assert list(cfg["genes"]) == ["toy_gene"]


class TestRejections:
    def test_missing_required_key(self, schema_path):
        with pytest.raises(ValueError, match="genome_fasta"):
            check({"target_gene_fasta": "t.fa"}, schema_path)

    def test_unknown_top_level_key_is_a_typo_guard(self, schema_path):
        with pytest.raises(ValueError, match="Invalid config"):
            check({**MINIMAL, "max_mismatchez": 3}, schema_path)

    def test_error_names_location(self, schema_path):
        with pytest.raises(ValueError, match=r"guide\.gc_min"):
            check({**MINIMAL, "guide": {"gc_min": 2}}, schema_path)

    @pytest.mark.parametrize("override", [
        {"max_mismatches": -1}, {"max_mismatches": 7}, {"threads": 0},
        {"max_mismatches": "three"}, {"seed": -1},
    ])
    def test_out_of_range_values(self, schema_path, override):
        with pytest.raises(ValueError):
            check({**MINIMAL, **override}, schema_path)

    def test_weights_must_sum_to_one(self, schema_path):
        bad = {**MINIMAL, "scoring": {"weights": {"gc": 0.5, "seed_gc": 0.5, "mfe": 0.5, "off_target": 0.5}}}
        with pytest.raises(ValueError, match="sum to 1.0"):
            check(bad, schema_path)

    def test_gc_min_above_gc_max(self, schema_path):
        with pytest.raises(ValueError, match="gc_min"):
            check({**MINIMAL, "guide": {"gc_min": 0.8, "gc_max": 0.5}}, schema_path)

    def test_seed_longer_than_spacer(self, schema_path):
        with pytest.raises(ValueError, match="seed_length"):
            check({**MINIMAL, "guide": {"seed_length": 25}}, schema_path)

    def test_unknown_searcher(self, schema_path):
        with pytest.raises(ValueError):
            check({**MINIMAL, "tools": {"off_target_searchers": ["magic"]}}, schema_path)

    def test_bad_gene_name(self, schema_path):
        with pytest.raises(ValueError):
            check({**MINIMAL, "genes": {"bad name/": "x.fa"}}, schema_path)

    def test_valid_genes_block(self, schema_path):
        cfg = check({**MINIMAL, "genes": {"blaNDM-1": "ndm.fa"}}, schema_path)
        assert cfg["genes"] == {"blaNDM-1": "ndm.fa"}

    def test_bad_color(self, schema_path):
        with pytest.raises(ValueError):
            check({**MINIMAL, "export": {"colors": {"header": "blue"}}}, schema_path)


class TestLoadConfig:
    def test_missing_file(self, tmp_path, schema_path):
        with pytest.raises(FileNotFoundError):
            load_config(str(tmp_path / "nope.yaml"), schema_path)

    def test_loads_yaml_from_disk(self, tmp_path, schema_path):
        path = tmp_path / "c.yaml"
        path.write_text(yaml.safe_dump({**MINIMAL, "threads": 2}))
        assert load_config(str(path), schema_path)["threads"] == 2

    def test_empty_file_reports_missing_required(self, tmp_path, schema_path):
        path = tmp_path / "c.yaml"
        path.write_text("")
        with pytest.raises(ValueError, match="genome_fasta"):
            load_config(str(path), schema_path)
