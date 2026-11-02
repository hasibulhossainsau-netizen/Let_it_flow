"""Cas-OFFinder implementation of ``OffTargetSearcher``."""

import os
import subprocess

import pandas as pd

from .base import HIT_COLUMNS, OffTargetSearcher

_RESULT_COLUMNS = ["query", "chromosome", "position", "strand",
                   "off_target_seq", "mismatches"]


class CasOFFinderSearcher(OffTargetSearcher):
    """Run Cas-OFFinder on spacer+PAM queries."""

    name = "cas-offinder"

    def __init__(self, device: str = "C", pam: str = "NGG",
                 spacer_length: int = 20, executable: str = "cas-offinder"):
        self.device = device
        self.pam = pam
        self.spacer_length = spacer_length
        self.executable = executable

    def version(self) -> str:
        out = subprocess.run([self.executable, "--version"],
                             capture_output=True, text=True)
        return (out.stdout or out.stderr).strip()

    def search(self, queries, genome, workdir, max_mismatches):
        os.makedirs(workdir, exist_ok=True)
        in_path = os.path.join(workdir, "casoffinder_input.txt")
        out_path = os.path.join(workdir, "casoffinder_results.txt")

        full = (queries["spacer"].str.upper() + queries["pam"].str.upper())
        with open(in_path, "w") as handle:
            handle.write(os.path.abspath(genome) + "\n")
            handle.write("N" * self.spacer_length + self.pam + "\n")
            for target in full:
                handle.write(f"{target}\t{max_mismatches}\n")

        result = subprocess.run([self.executable, in_path, self.device, out_path],
                                capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(
                f"cas-offinder failed (exit {result.returncode}): "
                f"{result.stderr.strip()}")

        try:
            hits = pd.read_csv(out_path, sep="\t", header=None,
                               names=_RESULT_COLUMNS)
        except pd.errors.EmptyDataError:
            return pd.DataFrame(columns=HIT_COLUMNS)

        ids = pd.DataFrame({"grna_id": queries["grna_id"].values,
                            "spacer": queries["spacer"].values,
                            "query": full.values})
        merged = hits.merge(ids, on="query", how="inner")
        return merged[HIT_COLUMNS].reset_index(drop=True)
