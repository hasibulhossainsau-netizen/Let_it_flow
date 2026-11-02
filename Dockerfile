# Dockerfile — reproducible gRNA pipeline image
#
# Build : docker build -t crispr-pipeline:latest .
# Run   : docker run --rm -v "$PWD":/pipeline -w /pipeline \
#             --user "$(id -u):$(id -g)" crispr-pipeline:latest \
#             snakemake --cores 4
#
# The base image tag is pinned; the toolchain comes from envs/crispr_pipeline.yaml
# (exact version pins), so rebuilding yields the same tool versions.
FROM mambaorg/micromamba:1.5.8

LABEL org.opencontainers.image.title="crispr-pipeline" \
      org.opencontainers.image.description="Cas9 gRNA design pipeline (Snakemake, BLAST+, Cas-OFFinder, ViennaRNA)" \
      org.opencontainers.image.licenses="see LICENSE"

# Install the pinned environment into the image's base env (no per-rule envs
# are needed inside the container: run snakemake WITHOUT --use-conda).
COPY --chown=$MAMBA_USER:$MAMBA_USER envs/crispr_pipeline.yaml /tmp/crispr_pipeline.yaml
RUN micromamba install --yes --name base --file /tmp/crispr_pipeline.yaml \
    && micromamba clean --all --yes \
    && rm /tmp/crispr_pipeline.yaml

# Makes `RUN`/CMD see the activated env (micromamba images activate via entrypoint).
ARG MAMBA_DOCKERFILE_ACTIVATE=1

# Bundle the pipeline code so the image also works stand-alone. Normally you
# bind-mount your working directory over /pipeline (see "Run" above) so that
# data/, results/ and logs/ live on the host.
WORKDIR /pipeline
COPY --chown=$MAMBA_USER:$MAMBA_USER . /pipeline

# Fixed seeds/ordering: matches scripts/run_guarded.sh
ENV PYTHONHASHSEED=42

ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]
CMD ["snakemake", "--cores", "4", "--snakefile", "workflow/Snakefile"]
