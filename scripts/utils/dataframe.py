"""DataFrame validation helpers."""

import pandas as pd


def require_columns(df: pd.DataFrame, columns, label: str = "table") -> None:
    """Raise ``ValueError`` if any of ``columns`` is missing from ``df``."""
    missing = [c for c in columns if c not in df.columns]
    if missing:
        raise ValueError(
            f"{label} is missing required column(s): {', '.join(missing)}")


def fill_missing_columns(df: pd.DataFrame, defaults: dict) -> pd.DataFrame:
    """Add any column in ``defaults`` that ``df`` lacks, with its default."""
    for col, value in defaults.items():
        if col not in df.columns:
            df[col] = value
    return df
