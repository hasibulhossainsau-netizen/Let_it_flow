#!/usr/bin/env python3
"""Build a self-contained toy working directory for Snakemake.

Copies the workflow code and drops the synthetic 1 kb genome / 100 bp gene
fixtures in the layout the pipeline expects, so that plain
``snakemake --cores 1 --use-conda`` can be run inside it.

    python tests/toy_workdir.py /tmp/toy && cd /tmp/toy && snakemake -n
"""
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TOY = ROOT / "tests" / "data" / "toy"


def build_workdir(dest) -> Path:
    work = Path(dest)
    work.mkdir(parents=True, exist_ok=True)
    for name in ("workflow", "scripts", "interfaces", "envs"):
        shutil.copytree(ROOT / name, work / name, dirs_exist_ok=True,
                        ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
    (work / "config").mkdir(exist_ok=True)
    shutil.copy(ROOT / "config" / "schema.yaml", work / "config" / "schema.yaml")
    shutil.copy(TOY / "config.yaml", work / "config" / "config.yaml")
    raw = work / "data" / "raw"
    raw.mkdir(parents=True, exist_ok=True)
    shutil.copy(TOY / "genome.fasta", raw / "toy_genome.fasta")
    shutil.copy(TOY / "gene.fasta", raw / "toy_gene.fasta")
    shutil.copy(TOY / "checksums.md5", raw / "checksums.md5")
    return work


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    print(build_workdir(sys.argv[1]))
