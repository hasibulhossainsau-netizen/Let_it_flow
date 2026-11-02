#!/usr/bin/env python3
"""Export the ranked gRNA table to CSV and a formatted Excel workbook.

Usage:
    python scripts/export_final.py --input ranked.csv --output results/prefix \
        [--config config/config.yaml]
"""

import argparse

import pandas as pd

from report_formatting import select_final_columns, style_workbook
from utils.config import load_config
from utils.runtime import event, get_logger, init_script, run_main

log = get_logger("export_final")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True,
                        help="Output prefix (no extension)")
    parser.add_argument("--config", default="config/config.yaml")
    args = parser.parse_args()

    cfg = load_config(args.config)
    export_cfg = cfg["export"]
    init_script(log, args, cfg, {"input": args.input})
    df_final = select_final_columns(pd.read_csv(args.input),
                                    export_cfg["columns"])

    csv_path, xlsx_path = args.output + ".csv", args.output + ".xlsx"
    df_final.to_csv(csv_path, index=False)
    event(log, "INFO", "CSV saved", path=csv_path, rows=len(df_final))

    sheet = export_cfg["sheet_name"]
    df_final.to_excel(xlsx_path, index=False, sheet_name=sheet)
    style_workbook(xlsx_path, sheet, export_cfg["colors"])
    event(log, "INFO", "Excel saved with formatting", path=xlsx_path)

    if "verdict" in df_final.columns:
        event(log, "INFO", "final verdict summary",
              verdicts={str(k): int(v) for k, v in
                        df_final["verdict"].value_counts().items()})


if __name__ == "__main__":
    run_main(main, "export_final")
