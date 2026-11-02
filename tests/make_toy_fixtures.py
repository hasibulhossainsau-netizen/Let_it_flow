#!/usr/bin/env python3
"""Regenerate the synthetic integration-test fixtures (deterministic).

Writes tests/data/toy/{genome.fasta,gene.fasta,checksums.md5}:
  * gene.fasta   - 100 bp random target (seeded)
  * genome.fasta - 1000 bp random "genome" with the gene embedded once, so the
                   target has exactly one on-target site and no repeats.

The fixtures are committed; run this only to change them:
    python tests/make_toy_fixtures.py [--seed N]
"""
import argparse
import hashlib
import random
from pathlib import Path

OUT = Path(__file__).resolve().parent / "data" / "toy"
GENE_LEN, GENOME_LEN, GENE_OFFSET = 100, 1000, 450
DEFAULT_SEED = 7


def wrap(seq, width=70):
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def build(seed):
    rng = random.Random(seed)
    gene = "".join(rng.choice("ACGT") for _ in range(GENE_LEN))
    genome = [rng.choice("ACGT") for _ in range(GENOME_LEN)]
    genome[GENE_OFFSET:GENE_OFFSET + GENE_LEN] = gene
    return gene, "".join(genome)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    args = parser.parse_args()
    gene, genome = build(args.seed)
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "gene.fasta").write_text(f">toy_gene synthetic 100 bp target\n{wrap(gene)}\n")
    (OUT / "genome.fasta").write_text(f">toy_chr synthetic 1 kb genome\n{wrap(genome)}\n")
    lines = []
    for name, dest in (("genome.fasta", "data/raw/toy_genome.fasta"),
                       ("gene.fasta", "data/raw/toy_gene.fasta")):
        digest = hashlib.md5((OUT / name).read_bytes()).hexdigest()
        lines.append(f"{digest}  {dest}")
    (OUT / "checksums.md5").write_text("\n".join(lines) + "\n")
    print(f"seed={args.seed}: wrote fixtures to {OUT}")


if __name__ == "__main__":
    main()
