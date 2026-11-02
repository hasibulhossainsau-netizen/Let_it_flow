"""Shared runtime for every pipeline script: determinism + structured logging.

* ``seed_everything`` seeds Python's ``random`` and NumPy (default 42).
* Logging is JSON-lines on **stderr**. Snakemake rules run scripts through
  ``scripts/run_guarded.sh``, which appends stdout+stderr to
  ``logs/<rule>/<gene>.log``, so the script's structured records, tool output
  and tracebacks all end up in that one file.
* ``init_script`` logs, once per run, the execution parameters, the seed,
  library/tool versions and statistics of every input file.
* ``run_main`` wraps ``main()``: logs start/finish with elapsed time and turns
  any uncaught exception into a structured error record + exit code 1.
"""

import json
import logging
import os
import platform
import random
import sys
import time
from datetime import datetime, timezone
from importlib import metadata

DEFAULT_SEED = 42
_START = time.monotonic()


# ── Determinism ──────────────────────────────────────────────────────────────
def seed_everything(seed: int = DEFAULT_SEED) -> int:
    """Seed ``random`` and NumPy (if importable); returns the seed used."""
    random.seed(seed)
    try:
        import numpy as np
        np.random.seed(seed)
    except ImportError:  # pragma: no cover - numpy ships with pandas
        pass
    return seed


# ── Structured logging ───────────────────────────────────────────────────────
class JsonFormatter(logging.Formatter):
    """One JSON object per line: ts, level, script, msg + any extra fields."""

    def format(self, record):
        entry = {
            "ts": datetime.fromtimestamp(record.created, timezone.utc)
                          .isoformat(timespec="milliseconds"),
            "level": record.levelname,
            "script": record.name,
            "msg": record.getMessage(),
        }
        entry.update(getattr(record, "fields", {}))
        if record.exc_info:
            entry["traceback"] = self.formatException(record.exc_info)
        return json.dumps(entry, default=str)


def get_logger(name: str) -> logging.Logger:
    """Logger that writes JSON lines to stderr (configured once per name)."""
    logger = logging.getLogger(name)
    if not logger.handlers:
        handler = logging.StreamHandler(sys.stderr)
        handler.setFormatter(JsonFormatter())
        logger.addHandler(handler)
        logger.setLevel(os.environ.get("GRNA_LOG_LEVEL", "INFO"))
        logger.propagate = False
    return logger


def event(logger: logging.Logger, level: str, msg: str, **fields) -> None:
    """Log ``msg`` at ``level`` with structured key/value ``fields``."""
    logger.log(getattr(logging, level.upper()), msg, extra={"fields": fields})


def file_stats(path) -> dict:
    """Size, mtime and (for text files) line count of ``path``."""
    path = str(path)
    if not os.path.isfile(path):
        return {"path": path, "exists": False}
    st = os.stat(path)
    stats = {"path": path, "exists": True, "size_bytes": st.st_size,
             "mtime": datetime.fromtimestamp(st.st_mtime, timezone.utc)
                              .isoformat(timespec="seconds")}
    binary = path.lower().endswith((".xlsx", ".gz", ".nsq", ".nin", ".nhr"))
    if not binary and st.st_size < 256 * 1024 * 1024:
        with open(path, "rb") as handle:
            stats["lines"] = sum(chunk.count(b"\n") for chunk in
                                 iter(lambda: handle.read(1 << 20), b""))
    return stats


def _lib_versions() -> dict:
    versions = {"python": platform.python_version()}
    for pkg in ("pandas", "numpy", "biopython", "PyYAML", "jsonschema",
                "openpyxl"):
        try:
            versions[pkg] = metadata.version(pkg)
        except metadata.PackageNotFoundError:
            pass
    return versions


def init_script(logger, args=None, cfg=None, inputs=None, tools=None) -> int:
    """Seed RNGs and log parameters, versions and input stats; returns seed.

    ``args``   argparse namespace (logged as execution parameters)
    ``cfg``    loaded config; ``cfg['seed']`` overrides the default of 42
    ``inputs`` {label: path} of input files to describe
    ``tools``  {name: obj_with_version()} external tools to record versions of
    """
    seed = seed_everything((cfg or {}).get("seed", DEFAULT_SEED))
    tool_versions = {}
    for name, tool in (tools or {}).items():
        try:
            tool_versions[name] = tool.version()
        except Exception as err:  # a missing --version must never kill a run
            tool_versions[name] = f"unavailable ({err})"
    event(logger, "INFO", "run parameters",
          params=vars(args) if args is not None else {},
          seed=seed, pythonhashseed=os.environ.get("PYTHONHASHSEED"),
          versions=_lib_versions(), tool_versions=tool_versions,
          cwd=os.getcwd(), argv=sys.argv)
    for label, path in (inputs or {}).items():
        event(logger, "INFO", f"input file: {label}", **file_stats(path))
    return seed


def run_main(main, name: str) -> None:
    """Run ``main()`` with start/finish records and structured crash logging."""
    logger = get_logger(name)
    event(logger, "INFO", "script started")
    try:
        main()
    except SystemExit as exc:
        code = exc.code if isinstance(exc.code, int) else (0 if exc.code is None else 1)
        event(logger, "INFO" if code == 0 else "ERROR", "script finished",
              exit_code=code, elapsed_s=round(time.monotonic() - _START, 3))
        raise
    except Exception as exc:
        logger.error("unhandled exception: %s", exc, exc_info=True)
        event(logger, "ERROR", "script finished", exit_code=1,
              elapsed_s=round(time.monotonic() - _START, 3))
        sys.exit(1)
    event(logger, "INFO", "script finished", exit_code=0,
          elapsed_s=round(time.monotonic() - _START, 3))
