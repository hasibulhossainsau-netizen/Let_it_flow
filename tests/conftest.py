"""Shared pytest setup: import paths and a loader for numbered scripts."""
import importlib.util
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"

# Scripts import siblings as top-level modules (``import scoring_metrics``,
# ``from utils.config import ...``), exactly as when run as
# ``python scripts/<name>.py``; ``utils`` also puts the repo root on sys.path
# so ``interfaces`` resolves.
for path in (SCRIPTS, ROOT):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))


def load_script(filename: str):
    """Import ``scripts/<filename>.py`` (works for names like 01_validate_inputs)."""
    spec = importlib.util.spec_from_file_location(
        f"script_{filename}", SCRIPTS / f"{filename}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="session")
def repo_root() -> Path:
    return ROOT


@pytest.fixture(scope="session")
def shipped_config() -> dict:
    """The repository's config/config.yaml, validated with defaults filled."""
    from utils.config import load_config
    return load_config(str(ROOT / "config" / "config.yaml"),
                       str(ROOT / "config" / "schema.yaml"))
