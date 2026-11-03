"""End-to-end toy run: 1 kb genome + 100 bp gene through the real workflow.

Runs exactly what a user runs, in a throw-away copy of the repository layout:

    snakemake --cores 1 --use-conda

Environment knobs:
  GRNA_IT_USE_CONDA=0      skip --use-conda (tools already on PATH, e.g. inside
                           the Docker image or an activated crispr_pipeline env)
  GRNA_CONDA_PREFIX=<dir>  where Snakemake caches the built rule environment
                           (default ~/.cache/crispr_pipeline_tests/conda), so
                           the ~450 MB env is built once, not per test run.
"""
import hashlib
import json
import os
import shutil
import subprocess
from pathlib import Path

import pandas as pd
import pytest

from toy_workdir import build_workdir

pytestmark = pytest.mark.integration

ROOT = Path(__file__).resolve().parents[2]
GENE = "toy_gene"
USE_CONDA = os.environ.get("GRNA_IT_USE_CONDA", "1") != "0"
CONDA_PREFIX = os.environ.get(
    "GRNA_CONDA_PREFIX", str(Path.home() / ".cache" / "crispr_pipeline_tests" / "conda"))

EXPECTED_COLUMNS = ["rank", "grna_id", "spacer", "pam", "strand", "position_nt",
                    "gc_content", "seed_gc", "mfe_kcal", "mfe_score",
                    "ot_0mm_adjusted", "ot_1mm", "ot_2mm", "ot_3mm",
                    "total_offtargets", "blast_offtarget_hits",
                    "efficiency_score", "verdict"]


def _skip_reason():
    if shutil.which("snakemake") is None:
        return "snakemake not on PATH"
    if USE_CONDA and shutil.which("conda") is None:
        return "conda not on PATH (needed for --use-conda; set GRNA_IT_USE_CONDA=0 to use PATH tools)"
    return None


def write_checksums(work):
    """Write the MD5s of the synthetic inputs into the workdir's checksums.md5.

    Computed from the files actually present, so the closed checksum_verify
    rule can never fail on a stale or line-ending-altered fixture checksum.
    """
    raw = work / "data" / "raw"
    lines = [f"{hashlib.md5((raw / name).read_bytes()).hexdigest()}  data/raw/{name}"
             for name in ("toy_genome.fasta", "toy_gene.fasta")]
    (raw / "checksums.md5").write_text("\n".join(lines) + "\n")


@pytest.fixture(scope="module")
def workdir(tmp_path_factory):
    reason = _skip_reason()
    if reason:
        pytest.skip(reason)
    work = build_workdir(tmp_path_factory.mktemp("toy_run"))
    write_checksums(work)
    return work


def run_snakemake(work, *extra):
    cmd = ["snakemake", "--cores", "1", *extra]
    if USE_CONDA:
        # Always use the conda frontend: Snakemake 8.5 defaults to `mamba env create
        # --quiet ...`, which fails against the mamba 2.x shipped by recent
        # Miniforge/setup-miniconda. conda (libmamba solver) is just as fast.
        cmd += ["--use-conda", "--conda-prefix", CONDA_PREFIX,
                "--conda-frontend", "conda"]
    env = {**os.environ, "PYTHONDONTWRITEBYTECODE": "1"}
    return subprocess.run(cmd, cwd=work, env=env, capture_output=True, text=True,
                          timeout=3000)


@pytest.fixture(scope="module")
def first_run(workdir):
    proc = run_snakemake(workdir)
    if proc.returncode != 0:
        # Full, untruncated output (run_snakemake uses text=True, so these are str).
        pytest.fail(
            f"snakemake failed (exit {proc.returncode})\n"
            f"command: {' '.join(proc.args)}\n"
            f"--- stdout ---\n{proc.stdout}\n--- stderr ---\n{proc.stderr}",
            pytrace=False)
    return proc


@pytest.fixture(scope="module")
def results(workdir, first_run):
    return pd.read_csv(workdir / "results" / f"{GENE}_final.csv")


def test_workflow_succeeds_and_produces_reports(workdir, first_run):
    assert (workdir / "results" / f"{GENE}_final.csv").stat().st_size > 0
    assert (workdir / "results" / f"{GENE}_final.xlsx").stat().st_size > 0


def test_final_table_schema_and_content(results):
    assert list(results.columns) == EXPECTED_COLUMNS
    assert 1 <= len(results) <= 30
    assert not results.isna().any().any(), "zero NaN anywhere in the final table"
    assert results["rank"].tolist() == list(range(1, len(results) + 1))
    assert results["grna_id"].is_unique
    assert results["spacer"].str.len().eq(20).all()
    assert results["efficiency_score"].between(0, 1).all()
    assert results["efficiency_score"].is_monotonic_decreasing


def test_rnafold_parser_fix_end_to_end(results):
    """No guide may carry the old silent-failure MFE of exactly 0.0 unless it is
    genuinely unstructured; padded energies like ( -8.30) must be parsed."""
    assert (results["mfe_kcal"] <= 0).all()
    assert (results["mfe_kcal"] < 0).any()


def test_gene_embedded_once_has_no_exact_offtargets(results):
    assert (results["ot_0mm_adjusted"] == 0).all()
    assert not results["verdict"].str.contains("exact off-target").any()


def test_every_guide_finds_its_own_site_with_a_true_ngg_pam(workdir, results):
    """Regression for the PAM-regex fix: spacer + NGG must round-trip through
    Cas-OFFinder, so each guide gets exactly one 0-mismatch hit (its own site)."""
    assert results["pam"].str.fullmatch("[ACGT]GG").all()
    hits = pd.read_csv(workdir / "data" / "processed" / GENE / "casoffinder_results.txt",
                       sep="\t", header=None)
    assert (hits[5] == 0).sum() == len(results)
    assert (results["total_offtargets"] == 1).all()
    assert (results["ot_0mm_adjusted"] == 0).all()


def test_all_sanity_flags_written(workdir):
    flags = workdir / "logs" / "flags" / GENE / "checks"
    stages = {p.stem for p in flags.glob("*.ok")}
    assert {"candidates_raw", "candidates_filtered", "casoffinder_results",
            "offtargets", "blast_queries", "blast_results", "blast_annotated",
            "rnafold", "ranked", "final_csv", "final_xlsx"} <= stages
    assert (workdir / "logs" / "flags" / "checksum_verify.ok").exists()


def test_logs_are_structured_and_record_tool_versions(workdir):
    log_dir = workdir / "logs"
    expected = ["checksum_verify/all.log", "validate_genome/all.log",
                f"validate_gene/{GENE}.log", f"pam_finder/{GENE}.log",
                f"casoffinder_run/{GENE}.log", f"blast_run/{GENE}.log",
                f"rnafold_score/{GENE}.log", f"rank_grnas/{GENE}.log",
                f"export_final/{GENE}.log"]
    for rel in expected:
        lines = [ln for ln in (log_dir / rel).read_text().splitlines() if ln.strip()]
        records = []
        for ln in lines:
            try:
                records.append(json.loads(ln))
            except json.JSONDecodeError:
                pass  # raw tool output (e.g. Cas-OFFinder progress) is allowed
        assert records, f"{rel}: no JSON records"
        assert records[0]["script"] == "run_guarded" and records[0]["msg"] == "started"
        assert records[-1]["msg"] == "finished: OK", rel
    started = json.loads((log_dir / f"blast_run/{GENE}.log").read_text().splitlines()[0])
    assert "blastn" in started["tool_version"]
    started = json.loads((log_dir / f"rnafold_score/{GENE}.log").read_text().splitlines()[0])
    assert "RNAfold" in started["tool_version"]


def test_second_run_is_a_noop(workdir, first_run):
    proc = run_snakemake(workdir)
    assert proc.returncode == 0
    assert "Nothing to be done" in proc.stdout + proc.stderr


def test_corrupted_input_is_caught_by_checksum_rule(workdir, first_run):
    gene = workdir / "data" / "raw" / "toy_gene.fasta"
    original = gene.read_text()
    gene.write_text(original + "ACGT\n")
    try:
        proc = run_snakemake(workdir, "--forceall", "--keep-going")
        assert proc.returncode != 0
        assert "MD5 mismatch" in proc.stdout + proc.stderr
    finally:
        gene.write_text(original)
