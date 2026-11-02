#!/usr/bin/env python3
"""Write the pre-filter-passing gRNA spacers to a FASTA file for BLAST."""
import argparse

import pandas as pd

from utils.fasta import write_fasta
from utils.runtime import event, get_logger, init_script, run_main

log = get_logger("04b_generate_fasta")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    init_script(log, args, None, {"input": args.input})

    df = pd.read_csv(args.input)
    if "pre_filter_pass" in df.columns:
        df = df[df["pre_filter_pass"]]
    write_fasta(args.output, zip(df["grna_id"], df["spacer"]))
    event(log, "INFO", "BLAST query FASTA written", output=args.output,
          n_sequences=len(df))


if __name__ == "__main__":
    run_main(main, "04b_generate_fasta")
