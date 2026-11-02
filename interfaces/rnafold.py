"""ViennaRNA RNAfold implementation of ``StructurePredictor``."""

import os
import re
import subprocess
import tempfile
from typing import Tuple

from .base import StructurePredictor

# RNAfold prints the energy right-aligned in a fixed-width field, e.g.
# "(-12.30)" but also "( -8.30)" and "(  0.00)". The ``\s*`` after "(" is
# essential: without it every MFE in (-10, 0] fails to parse.
_MFE_RE = re.compile(r"\(\s*(-?\d+\.\d+)\)\s*$")


def parse_rnafold_output(stdout: str) -> Tuple[float, str]:
    """Extract ``(mfe, dot_bracket)`` from RNAfold's standard output.

    Looks at the last non-empty line, ``<structure> (<mfe>)``. Raises
    ``ValueError`` if that line does not have the expected shape.
    """
    lines = [ln for ln in stdout.strip().split("\n") if ln.strip()]
    last = lines[-1].strip() if lines else ""
    match = _MFE_RE.search(last)
    if not match:
        raise ValueError(f"unparseable RNAfold output line: {last!r}")
    return float(match.group(1)), last.split()[0]


class RNAfoldPredictor(StructurePredictor):
    """Fold sequences with the ``RNAfold`` executable."""

    name = "rnafold"

    def __init__(self, timeout_s: int = 30, executable: str = "RNAfold"):
        self.timeout_s = timeout_s
        self.executable = executable

    def version(self) -> str:
        out = subprocess.run([self.executable, "--version"],
                             capture_output=True, text=True, timeout=self.timeout_s)
        return out.stdout.strip() or out.stderr.strip()

    def predict(self, rna_seq: str) -> Tuple[float, str]:
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fa",
                                         delete=False) as tmp:
            tmp.write(f">seq\n{rna_seq}\n")
            tmp_path = tmp.name
        try:
            result = subprocess.run(
                [self.executable, "--noPS", tmp_path],
                capture_output=True, text=True, timeout=self.timeout_s)
        except (OSError, subprocess.SubprocessError) as err:
            raise RuntimeError(f"RNAfold could not be run: {err}") from err
        finally:
            os.unlink(tmp_path)

        if result.returncode != 0:
            raise RuntimeError(
                f"RNAfold failed (exit {result.returncode}): "
                f"{result.stderr.strip() or 'no error output'}")
        try:
            return parse_rnafold_output(result.stdout)
        except ValueError as err:
            raise RuntimeError(f"RNAfold failed (exit 0): {err}") from err
