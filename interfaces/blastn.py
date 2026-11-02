"""BLASTn (blastn-short) implementation of ``OffTargetSearcher``."""

import os
import subprocess

import pandas as pd

from .base import HIT_COLUMNS, OffTargetSearcher

_BLAST_COLUMNS = ["qseqid", "sseqid", "pident", "length", "mismatch",
                  "gapopen", "qstart", "qend", "sstart", "send", "evalue",
                  "bitscore"]


class BlastnSearcher(OffTargetSearcher):
    """Homology screen of spacers against a genome with ``blastn``.

    Reported ``mismatches`` count aligned mismatches plus unaligned query bases.
    """

    name = "blastn"

    def __init__(self, task: str = "blastn-short", word_size: int = 7,
                 perc_identity: float = 80, threads: int = 1):
        self.task = task
        self.word_size = word_size
        self.perc_identity = perc_identity
        self.threads = threads

    def version(self) -> str:
        out = subprocess.run(["blastn", "-version"], capture_output=True,
                             text=True)
        return out.stdout.splitlines()[0] if out.stdout else out.stderr.strip()

    @staticmethod
    def _run(cmd):
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(
                f"{cmd[0]} failed (exit {result.returncode}): "
                f"{result.stderr.strip()}")

    def search(self, queries, genome, workdir, max_mismatches):
        os.makedirs(workdir, exist_ok=True)
        db = os.path.join(workdir, "blast_db")
        fasta = os.path.join(workdir, "blast_queries.fasta")
        out_path = os.path.join(workdir, "blast_results.tsv")

        with open(fasta, "w") as handle:
            for gid, spacer in zip(queries["grna_id"], queries["spacer"]):
                handle.write(f">{gid}\n{spacer}\n")

        self._run(["makeblastdb", "-in", genome, "-dbtype", "nucl", "-out", db])
        self._run(["blastn", "-query", fasta, "-db", db, "-out", out_path,
                   "-outfmt", "6 " + " ".join(_BLAST_COLUMNS),
                   "-word_size", str(self.word_size),
                   "-perc_identity", str(self.perc_identity),
                   "-num_threads", str(self.threads), "-task", self.task])

        try:
            hits = pd.read_csv(out_path, sep="\t", header=None,
                               names=_BLAST_COLUMNS)
        except pd.errors.EmptyDataError:
            return pd.DataFrame(columns=HIT_COLUMNS)

        spacers = queries.set_index("grna_id")["spacer"]
        hits["mismatches"] = (hits["mismatch"] + hits["gapopen"]
                              + hits["qseqid"].map(spacers.str.len())
                              - hits["length"])
        hits = hits[hits["mismatches"] <= max_mismatches]
        hits["grna_id"] = hits["qseqid"]
        hits["spacer"] = hits["qseqid"].map(spacers)
        hits["chromosome"] = hits["sseqid"]
        hits["position"] = hits[["sstart", "send"]].min(axis=1)
        hits["strand"] = (hits["sstart"] <= hits["send"]).map({True: "+", False: "-"})
        return hits[HIT_COLUMNS].reset_index(drop=True)
