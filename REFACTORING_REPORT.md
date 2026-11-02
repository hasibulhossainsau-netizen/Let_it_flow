# Refactoring Report

This report covers the five-phase refactoring of the gRNA design pipeline: what
changed architecturally, what breaks for existing users and how to migrate, and
two **bug fixes that change scientific results**: RNAfold MFE parsing (§3) and
the PAM regex, which matched `GGN` instead of `NGG` (§4).

**Contents**

1. [Architectural changes](#1-architectural-changes)
2. [Breaking changes & migration notes](#2-breaking-changes--migration-notes)
3. [Bug fix: RNAfold MFE parsing](#3-bug-fix-rnafold-mfe-parsing)
4. [Bug fix: the PAM regex matched `GGN` instead of `NGG`](#4-bug-fix-the-pam-regex-matched-ggn-instead-of-ngg)
5. [Other known limitations & observations](#5-other-known-limitations--observations)
6. [Verification](#6-verification)

---

## 1. Architectural changes

### Phase 1: configuration, modularity, plugins
* **Single source of truth.** Every tunable (GC window, homopolymer limit, seed
  definition, BLAST thresholds, sgRNA scaffold, MFE floor, score weights/curves,
  off-target tiers, verdict thresholds, export columns and colours, tool timeouts)
  moved out of the scripts into `config/config.yaml`.
* **Schema-validated config.** `config/schema.yaml` (JSON Schema) validates the
  config, fills defaults recursively, and rejects unknown keys. Cross-field
  checks (weights sum to 1, `gc_min ≤ gc_max`, `seed_length ≤ spacer_length`) live
  in `scripts/utils/config.py`. An old 5-key config still loads: every omitted
  section takes its default, and the defaults reproduce the previous behaviour.
* **Shared library.** `scripts/utils/` (`config`, `fasta`, `seq`, `dataframe`)
  replaces duplicated helper code; scoring maths became pure functions in
  `scripts/scoring_metrics.py`.
* **Script split.** `07_score_and_rank.py` and `08_export_final.py` were replaced
  by `rank_candidates.py`, `export_final.py`, `scoring_metrics.py` and
  `report_formatting.py`.
* **Plugin interfaces.** `interfaces/` defines `OffTargetSearcher` and
  `StructurePredictor` abstract bases with Cas-OFFinder, BLASTn and RNAfold
  implementations selected by name from `config.tools`. *(See
  [§5](#5-other-known-limitations--observations): only the structure predictor is
  wired into the workflow.)*

### Phase 2: Snakefile decomposition
* **Atomic rules.** Composite rules were split into single-responsibility steps:
  `casoffinder_prep → casoffinder_run → casoffinder_parse`,
  `blast_db`, `blast_fasta → blast_run → blast_parse`, plus `checksum_verify`,
  `validate_genome`, `validate_gene`, `pam_finder`, `filter_grnas`,
  `rnafold_score`, `rank_grnas`, `export_final`. `rnafold_score` stays one rule
  because its script drives RNAfold per spacer itself.
* **`{gene}` wildcard.** An optional `genes: {name: fasta}` block routes any
  number of targets through the same rules (`data/processed/{gene}/`,
  `results/{gene}_final.*`). Genome-level work (checksums, genome validation,
  BLAST DB) runs once.
* **Declarations on every rule:** `log:`, `threads:`, `resources:` (`mem_mb`,
  `runtime`) and `conda:`. External-tool rules also have `retries: 3`, with
  memory/runtime scaled by attempt.
* **Fail-fast wrapper.** `scripts/run_guarded.sh` appends to the log, records the
  exit code, dumps the log tail to the console on failure, and propagates the
  exit code so Snakemake can retry.
* **Integrity gates.** `checksum_verify` checks genome and genes against
  `data/raw/checksums.md5` (closed check: an unlisted file fails). A generic
  `validate_stage` rule (`scripts/validate_outputs.py`) checks, for each
  intermediate, non-empty output, schema alignment, zero NaN in key columns and a
  few invariants (unique IDs, scores in [0, 1], ranks are 1..N, RNAfold produced
  at least one structure). Downstream rules take the `.ok` flag as an input, so a
  bad file stops the DAG.

### Phase 3: script robustness
* **Input validation.** `01_validate_inputs.py` checks FASTA syntax, ID
  uniqueness, empty records, the ACGT alphabet (IUPAC tolerated in the genome up
  to a configurable fraction) and minimum lengths, with line/column-accurate,
  actionable messages.
* **Determinism.** Python and NumPy are seeded (`seed: 42`), `PYTHONHASHSEED=42`
  is set by the wrapper and images, and ranking sorts by
  `(efficiency_score desc, grna_id asc)` with a stable sort.
* **Structured logging.** Every script emits JSON lines (parameters, seed,
  library and tool versions, input-file statistics, start/finish with elapsed
  time, structured tracebacks) to stderr; the wrapper captures it together with
  raw tool output in `logs/<rule>/<gene>.log`. Sentinel flags moved to
  `logs/flags/`.
* **No more silent failures.** Zero candidates, or an unreadable/empty BLAST
  file, now exit with code 1 and an explanatory message.

### Phase 4: reproducibility
* **One pinned environment** `envs/crispr_pipeline.yaml` (exact versions for every
  direct dependency; `conda-forge` + `bioconda`, no `defaults`) and a fully
  resolved `envs/crispr_pipeline.lock.yaml` (`conda env export --no-builds`).
* **Containers.** `Dockerfile` (micromamba base) and `Singularity.def`
  (Apptainer) install exactly that environment.
* **Docs.** `docs/ENVIRONMENT.md`.

### Phase 5: tests, CI, docs
* **Unit tests** (`tests/unit`): scoring functions, FASTA parsing and input
  validation, config/schema validation, the RNAfold parser, ranking determinism,
  output sanity checks, PAM extraction.
* **Integration test** (`tests/integration`): a synthetic 1 kb genome with an
  embedded 100 bp gene run through `snakemake --cores 1 --use-conda`; it also
  checks outputs, structured logs, tool versions, idempotence and checksum
  failure handling.
* **CI** (`.github/workflows/ci.yml`): unit tests + config validation + dry-run,
  integration test, and the full suite inside the Docker image.
* **Docs.** Rewritten `README.md` (setup, configuration, testing, "How to
  Extend"), this report, and a regenerated workflow DAG
  (`example_output/pipeline_dag.png`, previously showing the pre-refactor rules).

---

## 2. Breaking changes & migration notes

### Breaking changes

| Area | Before | Now |
| :--- | :--- | :--- |
| **Config** | free-form | strictly validated: **unknown keys are errors**. New optional keys: `genes`, `checksums_file`, `seed`, `validation.max_genome_ambiguous_fraction`, and the sections listed in Phase 1. |
| **Required input** | none | `data/raw/checksums.md5` must list the genome **and every gene FASTA**, or the run stops at `checksum_verify`. |
| **Rule names** | `validate_inputs`, `casoffinder`, `parse_offtargets`, `blast_screen`, … | `checksum_verify`, `validate_genome`, `validate_gene`, `casoffinder_prep/run/parse`, `blast_fasta/run/parse`, `validate_stage`, … Update any `snakemake <rule>` targets, cluster profiles or `--until` flags. |
| **Intermediate paths** | `data/processed/*.csv` | `data/processed/<gene>/*.csv` (BLAST DB stays in `data/processed/genome_db/`). |
| **Result names** | `results/<prefix>_final.*` | `results/<gene>_final.*`; in single-gene mode `<gene>` = `output_prefix`, so **names are unchanged** unless you adopt `genes:`. |
| **Logs** | `logs/NN_name.log`, `logs/*.flag` | `logs/<rule>/<gene>.log` (JSON lines + raw tool output), flags in `logs/flags/`. |
| **Script stdout** | human-readable `print` output | JSON lines on **stderr**; nothing is printed to stdout. Anything that parsed script stdout must read the log instead. |
| **Scripts** | `07_score_and_rank.py`, `08_export_final.py` | removed; use `rank_candidates.py`, `export_final.py`. `01_validate_inputs.py` gained `--only/--genome/--gene`. |
| **Validation strictness** | short sequences only warned; no alphabet check | below-minimum length, non-ACGT characters in a gene, bad FASTA syntax and duplicate IDs are **errors**. |
| **Failure behaviour** | zero candidates / unreadable BLAST output returned quietly | exit code 1 with a message. |
| **Ranking** | unstable sort on `efficiency_score` | `(efficiency_score desc, grna_id asc)`: **tied guides can appear in a different order** than in older outputs. |
| **Guide candidates (PAM)** | `guide.pam_regex` `(?=(.{20})GG)` (GGN; spacers shifted by one base) | `(?=(.{20}).GG)` (true `NGG`): **different spacers, positions, PAMs and off-target counts; the candidate set changes** ([§4](#4-bug-fix-the-pam-regex-matched-ggn-instead-of-ngg)). Pin the old regex in your config to reproduce old candidates. |
| **RNAfold scores** | MFE ∈ (−10, 0] silently became 0.0 | parsed correctly; **scores, ranks and one verdict change** ([§3](#3-bug-fix-rnafold-mfe-parsing)). |
| **Environment** | `envs/grna_env.yaml` (`defaults` channel; jupyter, seaborn, entrez-direct, regex) | `envs/crispr_pipeline.yaml` (env name `crispr_pipeline`); `defaults` dropped; unused packages removed; `pyyaml`, `jsonschema`, `numpy`, `cas-offinder`, `pocl`, `pytest` added. |
| **Snakemake 8 + conda** | n/a | `--use-conda` needs `mamba` on `PATH` or `--conda-frontend conda`. |

### Migration steps

1. **Environment:** `conda env create -f envs/crispr_pipeline.yaml` (or build a
   container); retire `grna_pipeline` / `envs/grna_env.yaml`.
2. **Checksums:** `md5sum data/raw/<genome>.fasta data/raw/<gene>.fasta >> data/raw/checksums.md5`.
3. **Config:** an old 5-key config works unchanged. Remove any keys the schema
   rejects; the error names the offending path (e.g. `Invalid config at 'guide.gc_min'`).
   To process several genes, add a `genes:` block.
4. **Re-run from scratch.** Old `data/processed/*.csv` are not picked up by the new
   layout. Remove `data/processed/` and `logs/`, then run `snakemake`.
5. **Compare against old results with care:** expect a different candidate set
   (PAM fix), different MFE-derived scores/ranks (RNAfold fix) and a different
   order among tied guides. Results from before these fixes should not be reused.
6. **Update tooling** that referenced old rule names, log file names or script
   stdout.

---

## 3. Bug fix: RNAfold MFE parsing

**Status: fixed** (`interfaces/rnafold.py`).

### The bug
RNAfold prints the energy right-aligned in a fixed-width field: `(-12.30)` but
also `( -8.30)` and `(  0.00)`. The parser used the regex `\((-?\d+\.\d+)\)`,
which requires the number to follow `(` immediately, so **every MFE in (−10, 0]
failed to parse**. `06_rnafold_score.py` caught the resulting error, logged a
warning, and substituted `0.0`, which `mfe_score` maps to the **best possible
score (1.0)**. Weakly structured, entirely ordinary guides were therefore
silently promoted.

### The fix
```diff
-_MFE_RE = r"\((-?\d+\.\d+)\)"
+_MFE_RE = r"\(\s*(-?\d+\.\d+)\)\s*$"      # allow RNAfold's padding
```
Parsing was moved into a pure function, `parse_rnafold_output()`, so it can be
unit-tested without the binary. A real zero (`(  0.00)`) is now a genuine 0.0,
not a failure.

### Impact (measured on the original example, before the PAM fix)
Evidence that the bug was real: with the old parser, the pipeline reported 18
"unparseable output (exit 0)" failures on a 92-guide run, and in the old
`example_output/` **every one of those 18 guides carried `mfe_kcal = 0.0`,
including the old rank-1 and rank-2 guides**. Re-scoring that same 92-guide set
with the fixed parser (RNAfold 2.6.4):

| | Before | After |
| :--- | :---: | :---: |
| Guides with a bogus `mfe_kcal = 0.0` | 18 of 92 | 0 |
| MFE range | 0.0 to −20.0 | −7.5 to −20.0 |
| TOP CANDIDATE / CANDIDATE / REJECT | 74 / 6 / 12 | 73 / 7 / 12 |
| Rank-1 guide | `DQ915955_-_828` | `DQ915955_+_749` |

Exactly the 18 previously-failing guides changed (`mfe_kcal`, `mfe_score`,
`efficiency_score`, therefore ranks); one verdict moved. Spacers, GC, seed GC
and every off-target column were unchanged. (This intermediate comparison used
only the RNAfold-dependent columns and carried the original off-target columns
over, because the original genome is unavailable. The final `example_output/` was
subsequently regenerated end to end; see [§4](#regenerated-example).)

### Regression protection
`tests/unit/test_rnafold_parser.py` covers `(-12.30)`, `( -8.30)`, `(  0.00)`,
`( -0.30)`, multi-line output, garbage input, a stub-binary end-to-end path, and
(when installed) the real RNAfold. Reintroducing the old regex makes 9 of its
tests fail. The integration test also asserts that all toy MFEs are negative.

---

## 4. Bug fix: the PAM regex matched `GGN` instead of `NGG`

**Status: fixed** (`config/config.yaml` and the default in `config/schema.yaml`).

### The bug
The original `guide.pam_regex` was `(?=(.{20})GG)`, which requires `GG` to start
**immediately after** the 20-nt spacer. SpCas9's PAM is `NGG`: the spacer is
followed by one arbitrary base *and then* `GG`.

```
sequence:  AAAAAAAAAAAAAAAAAAAA T GG CCCCC
old      :  spacer = AAAAAAAAAAAAAAAAAAAT   pam = GGC     (spacer swallows the N)
fixed    :  spacer = AAAAAAAAAAAAAAAAAAAA   pam = TGG
```

1. **Spacers were shifted by one base**: the reported "spacer" was the last 19 nt
   of the true protospacer plus the PAM's `N`, and the reported "PAM" was `GG`
   plus a base that is not part of the PAM.
2. **Cas-OFFinder under-reported.** Queries are `spacer + PAM` against the
   pattern `20N + NGG`. With the old PAM (`GGx`) a guide's own site matched only
   when `x = G`, so off-target counts depended on an arbitrary PAM base.

### The fix
```diff
-  pam_regex: '(?=(.{20})GG)'
+  pam_regex: '(?=(.{20}).GG)'    # 20-nt spacer, then N, then GG; PAM = NGG
```
`pam_length: 3` is unchanged: the PAM is the three bases after the spacer. The
`{20}` must equal `guide.spacer_length`. Anyone who needs the old behaviour can
pin the old regex in their config (a unit test keeps that path working).

### Evidence
* **Toy integration run (1 kb genome, 100 bp gene, 5 guides):** Cas-OFFinder
  returned 2 hits before the fix (only the two `GGG` guides found their own site)
  and **5 hits after: every guide finds exactly its own site** (`total_offtargets
  = 1`, `ot_0mm_adjusted = 0`). The integration test now asserts this.
* **Real gene (blaCTX-M-15, 876 bp):** 99 raw candidates either way, but the
  spacers differ and the pre-filter pass count moves **92 → 90**. Against the
  public HS11286 genome Cas-OFFinder went from **0 hits to 25**. Every reported
  PAM now matches `[ACGT]GG`.
* Unit tests (`tests/unit/test_pam_finder.py`) cover all four `N` values, the
  "spacer must not swallow N" regression, both strands, and the legacy regex.

<a id="regenerated-example"></a>
### Regenerated example output
Because the guides themselves changed, `example_output/` could no longer be
patched and was **regenerated end to end** with the current pipeline (all fixes
applied) against the public *K. pneumoniae* HS11286 genome: RefSeq `NC_016845.1`
(chromosome) followed by plasmids `NC_016838.1`, `NC_016846.1`, `NC_016839.1`,
`NC_016840.1`, `NC_016847.1`, `NC_016841.1` (7 records, 5,764,167 bytes of FASTA,
MD5 `e2f05042654ce0a7a703c68c129a026e`), with the gene file
`data/raw/blaCTX-M-15_ref.fasta` and `config/config.yaml` as shipped.

| | Original example | Regenerated |
| :--- | :---: | :---: |
| Guides screened | 92 | 90 |
| TOP CANDIDATE / CANDIDATE / REJECT | 74 / 6 / 12 | 89 / 1 / 0 |
| Genome | not reproducible | public, accessions above |

The verdict mix is different mainly because the genome is different: the 12
rejections in the original example came from exact off-target hits that this
public genome does not produce. The example is now self-consistent and anyone can
reproduce it. `data/raw/checksums.md5` must contain the MD5 of whichever genome
file you use.

---

## 5. Other known limitations & observations

* **Off-target searcher plugins are not wired into the workflow.** The
  `OffTargetSearcher` classes exist, but the Snakefile runs the Cas-OFFinder and
  BLAST rule chains directly, so `tools.off_target_searchers` currently has no
  effect. Only `tools.structure_predictor` is live.
* **`total_offtargets` includes the guide's own on-target hit** whenever the
  gene is present in the genome (while `ot_0mm_adjusted` subtracts one). Since the
  PAM fix every guide finds its own site, so the on-target contributes a uniform +1
  to those guides and no longer biases the ranking, but absolute off-target counts
  still include it. (Pre-existing; surfaced by the toy run.)
* **The public HS11286 genome does not contain a near-exact copy of the target**
  (the regenerated example has no exact off-targets and its Cas-OFFinder output
  has 25 near-matches for 90 guides), so it exercises the specificity logic only
  weakly. Use a genome that carries the target to see rejections.
* **Biopython is pinned but no longer imported** (the new validator has its own
  parser); it was kept because it was requested.
* **Untested here:** `Singularity.def` (no Apptainer available) and the GitHub
  Actions workflow (not run on GitHub; its commands were rehearsed locally, see
  §6). `snakemake --lint` reports roughly 75 stylistic findings; CI runs it as
  informational only.
* **Platform:** pins were solved and tested on `linux-64` only.

---

## 6. Verification

| Check | Result |
| :--- | :--- |
| Unit tests (`pytest tests/unit tests/test_toy_fixtures.py`) | 162 passed |
| Integration test with the exact spec command, real `conda` building the env from `envs/crispr_pipeline.yaml` (`snakemake --cores 1 --use-conda`) | 8 passed (~93 s; before the PAM-fix test was added) |
| Integration test with tools on `PATH` (`GRNA_IT_USE_CONDA=0`), incl. the new "every guide finds its own site" test | 9 passed |
| Whole suite, tools on `PATH` | 171 passed (162 unit + 9 integration) |
| Full suite inside the Docker image (CI `container` job rehearsal, before the PAM fix) | 162 passed |
| Mutation check: old RNAfold regex restored | 9 parser tests fail, as intended |
| Full 50-job run, 2 genes, real BLAST+ / Cas-OFFinder / RNAfold | all jobs succeeded; Docker and micromamba results byte-identical |
| Full 27-job run on the public HS11286 genome with all fixes (regenerated `example_output/`) | succeeded; 90 guides, all PAMs `[ACGT]GG`, no `mfe_kcal = 0.0` |
| Environment solve + create from `envs/crispr_pipeline.yaml`; lock file re-solve (dry run) | OK (242 packages in the lock) |
| Determinism: 3 shuffled input orders, random `PYTHONHASHSEED` | byte-identical ranked output (89 of 92 rows tied) |

Not exercised: the workflow on GitHub-hosted runners, `apptainer build`.
