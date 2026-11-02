#!/usr/bin/env python3
import argparse

import pandas as pd

from utils.config import load_config
from utils.runtime import event, get_logger, init_script, run_main

log = get_logger("05_parse_blast")

def parse_args():
    parser = argparse.ArgumentParser(description="Parse BLAST results for CRISPR gRNA off-targets.")
    parser.add_argument("--blast", required=True, help="Path to blast_results.tsv")
    parser.add_argument("--grnas", required=True, help="Path to grnas_with_offtargets.csv")
    parser.add_argument("--output", required=True, help="Path to output final csv")
    parser.add_argument("--config", default="config/config.yaml")
    return parser.parse_args()

def main():
    args = parse_args()
    cfg = load_config(args.config)
    blast_cfg = cfg["blast"]
    init_script(log, args, cfg, {"blast": args.blast, "grnas": args.grnas})

    # Load BLAST results (outfmt 6 defaults)
    blast_cols = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", 
                  "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
    
    try:
        blast_df = pd.read_csv(args.blast, sep="\t", names=blast_cols)
    except Exception as e:
        log.error("Cannot read BLAST file %s: %s. An empty file means BLAST "
                  "found no hits - check the BLAST DB and query FASTA.",
                  args.blast, e)
        raise SystemExit(1)

    # Strict filtering: sequence alignment length/identity thresholds from config
    # This captures true potential off-targets
    valid_hits = blast_df[(blast_df["pident"] >= blast_cfg["min_pident"])
                           & (blast_df["length"] >= blast_cfg["min_length"])]

    # Group by gRNA and count total hits
    # Instead of just subtracting 1, we drop hits that represent the actual gene copies
    # Real off-targets are usually lower bitscore than the absolute 100% matches, 
    # but since it's multicopy, we count how many true high-risk hits exist.
    blast_counts = {}
    for grna, group in valid_hits.groupby("qseqid"):
        # If a gRNA has length=20 and pident>=95, it's an on-target copy. 
        # We find how many such on-target copies exist in the genome.
        on_target_copies = len(group[(group["length"] == blast_cfg["on_target_length"])
                                    & (group["pident"] >= blast_cfg["on_target_min_pident"])])
        
        # If no perfect copy found (rare), assume 1
        if on_target_copies == 0:
            on_target_copies = 1
            
        total_hits = len(group)
        # Off-targets = Total valid hits minus the valid on-target copies
        real_off_targets = max(0, total_hits - on_target_copies)
        blast_counts[grna] = real_off_targets

    # Load Cas-OFFinder results
    grna_df = pd.read_csv(args.grnas)

    # Map the true BLAST off-target counts
    grna_df["blast_offtarget_hits"] = grna_df["grna_id"].map(blast_counts).fillna(0).astype(int)

    # Save the updated candidates
    grna_df.to_csv(args.output, index=False)
    event(log, "INFO", "BLAST data merged", output=args.output,
          blast_rows=len(blast_df), valid_hits=len(valid_hits),
          grnas=len(grna_df))

if __name__ == "__main__":
    run_main(main, "05_parse_blast")