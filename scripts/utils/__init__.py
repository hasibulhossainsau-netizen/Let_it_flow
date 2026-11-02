"""Shared helpers for the gRNA design pipeline scripts.

Importing this package also puts the repository root on ``sys.path`` so that
scripts launched as ``python scripts/<name>.py`` can import ``interfaces``.
"""

import sys
from pathlib import Path

_ROOT = str(Path(__file__).resolve().parents[2])
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)
