"""Table selection and Excel styling for the final gRNA report."""

import pandas as pd
from openpyxl import load_workbook
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter

VERDICT_COLOR_KEYS = {
    "TOP CANDIDATE": "top_candidate",
    "CANDIDATE": "candidate",
    "MARGINAL": "marginal",
}


def select_final_columns(df: pd.DataFrame, columns: list) -> pd.DataFrame:
    """Keep the configured columns (in order) that exist in ``df``."""
    return df[[c for c in columns if c in df.columns]].copy()


def _fill(hex_color: str) -> PatternFill:
    return PatternFill("solid", fgColor=hex_color)


def style_workbook(xlsx_path: str, sheet_name: str, colors: dict) -> None:
    """Apply header styling, column widths, verdict row colours and a frozen
    header to the workbook at ``xlsx_path`` (modified in place)."""
    wb = load_workbook(xlsx_path)
    ws = wb[sheet_name]

    for cell in ws[1]:
        cell.fill = _fill(colors["header"])
        cell.font = Font(bold=True, color="FFFFFF")
        cell.alignment = Alignment(horizontal="center")

    for col_idx in range(1, ws.max_column + 1):
        max_len = max(
            (len(str(ws.cell(row=r, column=col_idx).value or ""))
             for r in range(1, ws.max_row + 1)),
            default=10,
        )
        ws.column_dimensions[get_column_letter(col_idx)].width = min(max_len + 4, 30)

    headers = [c.value for c in ws[1]]
    if "verdict" in headers:
        verdict_idx = headers.index("verdict")
        for row in ws.iter_rows(min_row=2, max_row=ws.max_row):
            key = VERDICT_COLOR_KEYS.get(str(row[verdict_idx].value or ""), "reject")
            fill = _fill(colors[key])
            for cell in row:
                cell.fill = fill

    ws.freeze_panes = "A2"
    wb.save(xlsx_path)
