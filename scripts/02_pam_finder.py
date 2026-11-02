#!/usr/bin/env python3
"""
PAM finder for SpCas9.
Scans forward and reverse-complement strands of the target gene and outputs a
CSV of candidate gRNA spacers. PAM regex, spacer length, GC bounds, homopolymer
limit and seed settings come from the ``guide`` section of config/config.yaml.

Usage:
    python scripts/02_pam_finder.py \
        --target data/raw/blaCTX-M-15_ref.fasta \
        --output data/processed/candidate_grnas_raw.csv
"""

import argparse
import csv
import re

from utils.config import load_config
from utils.fasta import read_fasta
from utils.runtime import event, get_logger, init_script, run_main
from utils.seq import gc_content, has_homopolymer, reverse_complement, seed_gc


log = get_logger("02_pam_finder")


def extract_grnas(seq: str, strand: str, seq_id: str, guide: dict) -> list:
    """Extract all candidate gRNAs from one strand sequence."""
    spacer_len = guide["spacer_length"]
    candidates = []
    for m in re.finditer(guide["pam_regex"], seq.upper()):
        spacer = m.group(1)
        pam_pos = m.start() + spacer_len
        pam = seq[pam_pos:pam_pos + guide["pam_length"]].upper()
        gc = gc_content(spacer)
        pass_gc = guide["gc_min"] <= gc <= guide["gc_max"]
        pass_hom = not has_homopolymer(spacer, guide["max_homopolymer"])
        sgc = seed_gc(spacer, guide["seed_length"])

        candidates.append({
            "grna_id": f"{seq_id}_{strand}_{pam_pos}",
            "spacer": spacer,
            "pam": pam,
            "strand": strand,
            "position_nt": pam_pos,
            "gc_content": round(gc, 3),
            "pass_gc": pass_gc,
            "pass_homopol": pass_hom,
            "pre_filter_pass": pass_gc and pass_hom,
            "seed_gc": round(sgc, 3),
            "good_seed": sgc >= guide["seed_min_gc"],
        })
    return candidates


def main():
    parser = argparse.ArgumentParser(description="PAM Finder — SpCas9")
    parser.add_argument("--target", required=True, help="Target gene FASTA")
    parser.add_argument("--output", required=True, help="Output CSV path")
    parser.add_argument("--config", default="config/config.yaml")
    args = parser.parse_args()

    cfg = load_config(args.config)
    guide = cfg["guide"]
    init_script(log, args, cfg, {"target": args.target})
    all_candidates = []

    for rec_id, fwd_seq in read_fasta(args.target):
        rev_seq = reverse_complement(fwd_seq)
        seq_id = rec_id.split(".")[0]

        fwd_hits = extract_grnas(fwd_seq, "+", seq_id, guide)
        rev_hits = extract_grnas(rev_seq, "-", seq_id, guide)
        event(log, "INFO", "scanned record", record=seq_id, length_bp=len(fwd_seq),
              forward_candidates=len(fwd_hits), reverse_candidates=len(rev_hits))
        all_candidates.extend(fwd_hits)
        all_candidates.extend(rev_hits)

    if not all_candidates:
        log.error("No gRNA candidates found. Check the input sequence and "
                  "guide.pam_regex in config.yaml.")
        raise SystemExit(1)

    with open(args.output, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(all_candidates[0].keys()))
        writer.writeheader()
        writer.writerows(all_candidates)

    passed = sum(1 for c in all_candidates if c["pre_filter_pass"])
    event(log, "INFO", "candidates written", total=len(all_candidates),
          passed_pre_filter=passed, output=args.output)


if __name__ == "__main__":
    run_main(main, "02_pam_finder")
