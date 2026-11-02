"""Abstract base classes for swappable external-tool back ends."""

from abc import ABC, abstractmethod
from typing import Tuple

import pandas as pd

HIT_COLUMNS = ["grna_id", "spacer", "chromosome", "position", "strand",
               "mismatches"]


class OffTargetSearcher(ABC):
    """Finds genomic sites resembling each gRNA (candidate off-targets)."""

    #: Registry key used in ``tools.off_target_searchers`` of the config.
    name: str = ""

    @abstractmethod
    def version(self) -> str:
        """Return the version string of the underlying tool."""

    @abstractmethod
    def search(self, queries: pd.DataFrame, genome: str, workdir: str,
               max_mismatches: int) -> pd.DataFrame:
        """Search ``genome`` for sites similar to each query gRNA.

        ``queries`` has columns ``grna_id``, ``spacer`` and ``pam``.
        Returns one row per hit with the columns in ``HIT_COLUMNS``; an empty
        DataFrame (with those columns) means no hits. Intermediate files are
        written under ``workdir``. Raises ``RuntimeError`` if the tool fails.
        """


class StructurePredictor(ABC):
    """Predicts RNA secondary structure and its minimum free energy."""

    #: Registry key used in ``tools.structure_predictor`` of the config.
    name: str = ""

    @abstractmethod
    def version(self) -> str:
        """Return the version string of the underlying tool."""

    @abstractmethod
    def predict(self, rna_seq: str) -> Tuple[float, str]:
        """Return ``(mfe_kcal_per_mol, dot_bracket_structure)`` for ``rna_seq``.

        Raises ``RuntimeError`` if the tool fails or output cannot be parsed.
        """
