# Reproducible environments

This pipeline needs five external toolchains (Python data stack, BLAST+,
ViennaRNA, Cas-OFFinder, Snakemake). Everything is pinned in **one** file,
[`envs/crispr_pipeline.yaml`](../envs/crispr_pipeline.yaml), which is the single
source of truth for local installs, `snakemake --use-conda`, Docker and
Singularity/Apptainer.

| Package      | Pin      | Provides                                   |
|--------------|----------|--------------------------------------------|
| python       | 3.11.9   | interpreter                                |
| pandas       | 2.2.1    | tables                                     |
| numpy        | 1.26.4   | scoring math                               |
| openpyxl     | 3.1.2    | Excel export                               |
| pyyaml       | 6.0.1    | config loading                             |
| jsonschema   | 4.21.1   | `config/schema.yaml` validation            |
| biopython    | 1.83     | pinned tool dependency (not imported by the current scripts) |
| blast        | 2.15.0   | `blastn`, `makeblastdb`                    |
| viennarna    | 2.6.4    | `RNAfold`                                  |
| cas-offinder | 2.4.1    | off-target search                          |
| pocl         | 5.0      | CPU OpenCL runtime for Cas-OFFinder (`tools.cas_offinder_device: C`) |
| snakemake    | 8.5.3    | workflow engine                            |
| pytest       | 8.2.0    | test runner (`tests/`)                     |

Channels are `conda-forge` then `bioconda`; the `defaults` channel is
intentionally not used. Platform: `linux-64` (WSL2 works).

## 1. Create the environment (conda / mamba / micromamba)

```bash
# any one of:
conda  env create -f envs/crispr_pipeline.yaml
mamba  env create -f envs/crispr_pipeline.yaml
micromamba create -f envs/crispr_pipeline.yaml

conda activate crispr_pipeline
```

Sanity check:

```bash
python -c "import pandas, numpy, openpyxl, yaml, jsonschema"
blastn -version | head -1      # blastn: 2.15.0+
RNAfold --version              # RNAfold 2.6.4
cas-offinder 2>&1 | head -1    # Cas-OFFinder v2.4.1 ...
snakemake --version            # 8.5.3
```

Run the pipeline from the repository root:

```bash
snakemake --cores 4                 # tools already on PATH (activated env)
snakemake --cores 4 --use-conda     # let Snakemake build the env per rule
```

Before the first run, `data/raw/checksums.md5` must list the genome and every
gene FASTA (`md5sum data/raw/<genome>.fasta >> data/raw/checksums.md5`).

## 2. Export / lock the environment (`--no-builds`)

The YAML above pins *direct* dependencies. A fully resolved snapshot of every
transitive package is committed as
[`envs/crispr_pipeline.lock.yaml`](../envs/crispr_pipeline.lock.yaml) (242
packages, linux-64); recreate it exactly with
`conda env create -f envs/crispr_pipeline.lock.yaml`.

It was produced by exporting the solved environment **without build strings**,
which keeps the file portable across minor rebuilds of the same package version.
Regenerate it whenever `crispr_pipeline.yaml` changes:

```bash
conda activate crispr_pipeline
conda env export --no-builds | grep -v '^prefix:' > envs/crispr_pipeline.lock.yaml
# micromamba equivalent:
micromamba env export --no-builds -n crispr_pipeline | grep -v '^prefix:' > envs/crispr_pipeline.lock.yaml
```

Notes:

* Drop the `prefix:` line (done by `grep -v` above); it is a local path.
* Check the `channels:` order in the exported file; keep `conda-forge` first.
* Rebuild from the lock file with
  `conda env create -n crispr_pipeline_locked -f envs/crispr_pipeline.lock.yaml`.
* `--no-builds` is portable but not bit-for-bit. For an exact, checksummed lock
  use `conda-lock` (`conda-lock -f envs/crispr_pipeline.yaml -p linux-64`) and commit the
  generated `conda-lock.yml`.
* `--from-history` instead exports only what you explicitly installed, which is
  roughly the hand-written file again.

### Updating a pin

1. Edit the version in `envs/crispr_pipeline.yaml`.
2. Confirm it solves: `micromamba create -n tmp -f envs/crispr_pipeline.yaml --dry-run`.
3. Re-run the pipeline and compare results with `example_output/`.
4. Rebuild the images (below) and regenerate the lock file if you keep one.

## 3. Docker

```bash
docker build -t crispr-pipeline:latest .

# bind-mount the project so data/, results/ and logs/ stay on the host
docker run --rm -v "$PWD":/pipeline -w /pipeline \
    --user "$(id -u):$(id -g)" crispr-pipeline:latest \
    snakemake --cores 4
```

* The image is `mambaorg/micromamba:1.5.8` + `envs/crispr_pipeline.yaml`.
* Do **not** pass `--use-conda` inside the container; the tools are already installed.
* `--user "$(id -u):$(id -g)"` keeps files in `results/` and `logs/` owned by you.
* The genome FASTA is not in the repository; place it under `data/raw/` on the host first.

## 4. Singularity / Apptainer (HPC)

```bash
# build from the repository root (the definition copies envs/crispr_pipeline.yaml)
apptainer build crispr-pipeline.sif Singularity.def
#   without root:            apptainer build --fakeroot crispr-pipeline.sif Singularity.def
#   or reuse the Docker image: apptainer build crispr-pipeline.sif docker-daemon://crispr-pipeline:latest

apptainer exec --bind "$PWD":/pipeline --pwd /pipeline crispr-pipeline.sif \
    snakemake --cores 4
```

On a scheduler, wrap the same `apptainer exec` line in your job script (or use
Snakemake's Slurm executor and pass `--software-deployment-method apptainer`
with `--apptainer-args "--bind ..."`). `singularity` can be substituted for
`apptainer` in every command.

## Determinism checklist

* tool/library versions: pinned in `envs/crispr_pipeline.yaml`;
* random seeds: `seed: 42` in `config/config.yaml` (Python + NumPy), `PYTHONHASHSEED=42` set by `scripts/run_guarded.sh` and the images;
* ranking ties: broken by `grna_id`;
* inputs: verified against `data/raw/checksums.md5` before anything runs;
* every run logs versions and parameters to `logs/<rule>/<gene>.log`.
