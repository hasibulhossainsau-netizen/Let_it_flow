"""Configuration loading, default filling and schema validation."""

import copy
from pathlib import Path

import jsonschema
import yaml

DEFAULT_CONFIG = "config/config.yaml"
DEFAULT_SCHEMA = "config/schema.yaml"


def _read_yaml(path):
    with open(path) as handle:
        return yaml.safe_load(handle)


def apply_schema_defaults(schema: dict, data: dict) -> dict:
    """Fill missing keys of ``data`` with the ``default`` values in ``schema``.

    Nested objects are handled recursively, so a config that omits a whole
    section (e.g. an old 5-key config) still receives every default.
    The input is modified in place and also returned.
    """
    for key, sub in schema.get("properties", {}).items():
        if key not in data and "default" in sub:
            data[key] = copy.deepcopy(sub["default"])
        if sub.get("type") == "object" and isinstance(data.get(key), dict):
            apply_schema_defaults(sub, data[key])
    return data


def validate_config(config: dict, schema_path: str = DEFAULT_SCHEMA) -> dict:
    """Apply schema defaults, validate, and run cross-field checks.

    Raises ``ValueError`` with an actionable message on any problem.
    """
    schema = _read_yaml(schema_path)
    apply_schema_defaults(schema, config)
    try:
        jsonschema.validate(config, schema)
    except jsonschema.ValidationError as err:
        location = ".".join(str(p) for p in err.absolute_path) or "<root>"
        raise ValueError(f"Invalid config at '{location}': {err.message}") from err

    weights = sum(config["scoring"]["weights"].values())
    if abs(weights - 1.0) > 1e-9:
        raise ValueError(f"scoring.weights must sum to 1.0 (got {weights})")
    guide = config["guide"]
    if guide["gc_min"] > guide["gc_max"]:
        raise ValueError("guide.gc_min must not exceed guide.gc_max")
    if guide["seed_length"] > guide["spacer_length"]:
        raise ValueError("guide.seed_length must not exceed guide.spacer_length")
    return config


def load_config(path: str = DEFAULT_CONFIG,
                schema_path: str = DEFAULT_SCHEMA) -> dict:
    """Load a YAML config file, fill defaults and validate it."""
    if not Path(path).is_file():
        raise FileNotFoundError(f"Config file not found: {path}")
    return validate_config(_read_yaml(path) or {}, schema_path)
