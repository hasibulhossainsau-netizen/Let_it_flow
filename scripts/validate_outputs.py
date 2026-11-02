#!/usr/bin/env python3
"""Post-step sanity checks for pipeline intermediates and results.

Each ``--stage`` knows what its file should look like:
  * the file exists and is non-empty (unless the stage legitimately may be),
  * schema alignment: required columns / field counts are present,
  * key columns contain zero NaN / empty values,
  * a few cheap invariants (unique IDs, scores within [0, 1], ...).

All problems are collected and reported together; exit code 1 on any failure.

Usage:
    python scripts/validate_outputs.py --stage ranked --input ranked.csv \
        [--config config/config.yaml]
"""
import argparse
import os
import sys

import pandas as pd

from utils.config import load_config
from utils.runtime import event, get_logger, init_script, run_main

log = get_logger("validate_outputs")

CANDIDATE_COLS = ["grna_id", "spacer", "pam", "strand", "position_nt",
                  "gc_content", "pass_gc", "pass_homopol", "pre_filter_pass",
                  "seed_gc", "good_seed"]
OFFTARGET_COLS = ["ot_0mm", "ot_1mm", "ot_2mm", "ot_3mm",
                  "total_offtargets", "ot_0mm_adjusted"]
UNIT_INTERVAL = ["gc_content", "seed_gc", "mfe_score", "gc_score",
                 "seed_score", "ot_score", "efficiency_score"]
NON_NEGATIVE = ["position_nt"] + OFFTARGET_COLS + ["blast_offtarget_hits"]

# stage -> (kind, required columns, extras)
CSV_STAGES = {
    "candidates_raw":      CANDIDATE_COLS,
    "candidates_filtered": CANDIDATE_COLS,
    "offtargets":          CANDIDATE_COLS + OFFTARGET_COLS,
    "blast_annotated":     CANDIDATE_COLS + OFFTARGET_COLS
                           + ["blast_offtarget_hits"],
    "rnafold":             CANDIDATE_COLS + OFFTARGET_COLS
                           + ["blast_offtarget_hits", "mfe_kcal", "mfe_score"],
    "ranked":              CANDIDATE_COLS + OFFTARGET_COLS
                           + ["blast_offtarget_hits", "mfe_kcal", "mfe_score",
                              "rank", "gc_score", "seed_score", "ot_score",
                              "efficiency_score", "verdict"],
    "final_csv":           None,     # taken from export.columns in the config
    "final_xlsx":          None,
}
TEXT_STAGES = ("casoffinder_results", "blast_results", "blast_queries")
ALL_STAGES = sorted(list(CSV_STAGES) + list(TEXT_STAGES))


def check_table(df, stage, required, cfg, errors):
    if df.empty:
        errors.append("table has a header but no rows")
        return
    missing = [c for c in required if c not in df.columns]
    if missing:
        errors.append(f"missing required column(s): {missing}")
    present = [c for c in required if c in df.columns]

    # Zero NaN / empty values in every key column.
    for col in present:
        n_nan = int(df[col].isna().sum())
        if n_nan:
            errors.append(f"column '{col}' has {n_nan} NaN/empty value(s)")

    if "grna_id" in df.columns and df["grna_id"].duplicated().any():
        errors.append("column 'grna_id' contains duplicates")
    if stage == "candidates_filtered" and "pre_filter_pass" in df.columns:
        if not (df["pre_filter_pass"] == True).all():  # noqa: E712
            errors.append("rows with pre_filter_pass != True survived filtering")

    for col in UNIT_INTERVAL:
        if col in present:
            vals = pd.to_numeric(df[col], errors="coerce")
            if vals.notna().all() and not vals.between(0, 1).all():
                errors.append(f"column '{col}' has values outside [0, 1]")
            elif vals.isna().any() and not df[col].isna().any():
                errors.append(f"column '{col}' is not numeric")
    for col in NON_NEGATIVE:
        if col in present:
            vals = pd.to_numeric(df[col], errors="coerce")
            if vals.notna().all() and (vals < 0).any():
                errors.append(f"column '{col}' has negative values")

    if stage in ("rnafold", "ranked") and "mfe_kcal" in df.columns:
        if (pd.to_numeric(df["mfe_kcal"], errors="coerce") > 0).any():
            errors.append("column 'mfe_kcal' has positive MFE values")
    if stage == "rnafold" and "structure" in df.columns:
        if df["structure"].isna().all():
            errors.append("RNAfold produced no structure for ANY spacer "
                          "(every prediction failed)")
    if stage == "ranked" and "rank" in df.columns:
        if sorted(df["rank"].tolist()) != list(range(1, len(df) + 1)):
            errors.append("column 'rank' is not a permutation of 1..N")


def check_text(path, stage, cfg, errors):
    with open(path) as handle:
        lines = [ln.rstrip("\n") for ln in handle if ln.strip()]

    if stage == "casoffinder_results":
        # Zero hits is a legitimate Cas-OFFinder outcome -> empty is allowed.
        if not lines:
            log.warning("Cas-OFFinder reported zero hits (empty result file)")
        for i, ln in enumerate(lines, 1):
            if len(ln.split("\t")) != 6:
                errors.append(f"line {i}: expected 6 tab-separated fields")
                break
    elif stage == "blast_results":
        if not lines:
            errors.append("BLAST produced no hits at all (is the gene in the "
                          "genome, and is the BLAST DB intact?)")
        n_fields = len(cfg["blast"]["outfmt"].split()) - 1
        for i, ln in enumerate(lines, 1):
            if len(ln.split("\t")) != n_fields:
                errors.append(f"line {i}: expected {n_fields} tab-separated "
                              f"fields (blast.outfmt)")
                break
    elif stage == "blast_queries":
        if not lines:
            errors.append("FASTA is empty")
        elif not lines[0].startswith(">"):
            errors.append("FASTA does not start with a '>' header")
        else:
            seq_lens, n = [], -1
            for ln in lines:
                if ln.startswith(">"):
                    seq_lens.append(0)
                    n += 1
                elif n >= 0:
                    seq_lens[n] += len(ln.strip())
            if any(length == 0 for length in seq_lens):
                errors.append("FASTA contains a record with an empty sequence")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", required=True, choices=ALL_STAGES)
    parser.add_argument("--input", required=True)
    parser.add_argument("--config", default="config/config.yaml")
    args = parser.parse_args()

    cfg = load_config(args.config)
    init_script(log, args, cfg, {"input": args.input})
    errors = []

    if not os.path.isfile(args.input):
        errors.append("file does not exist")
    elif os.path.getsize(args.input) == 0 and args.stage != "casoffinder_results":
        errors.append("file is empty (0 bytes)")
    elif args.stage in TEXT_STAGES:
        check_text(args.input, args.stage, cfg, errors)
    else:
        required = CSV_STAGES[args.stage] or cfg["export"]["columns"]
        try:
            if args.stage == "final_xlsx":
                df = pd.read_excel(args.input, sheet_name=cfg["export"]["sheet_name"])
            else:
                df = pd.read_csv(args.input)
        except Exception as err:  # unreadable / wrong sheet / bad CSV
            errors.append(f"cannot be read: {err}")
        else:
            check_table(df, args.stage, required, cfg, errors)
            if not errors:
                event(log, "INFO", "table checked", rows=len(df),
                      cols=len(df.columns))

    if errors:
        event(log, "ERROR", "validation failed", stage=args.stage,
              input=args.input, n_errors=len(errors))
        for msg in errors:
            log.error("%s: %s", args.input, msg)
        sys.exit(1)
    event(log, "INFO", "validation passed", stage=args.stage, input=args.input)


if __name__ == "__main__":
    run_main(main, "validate_outputs")
