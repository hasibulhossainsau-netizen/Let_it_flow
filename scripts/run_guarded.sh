#!/usr/bin/env bash
# Run a command, append its stdout/stderr to a log, and fail fast.
#
#   [VERSION_CMD='tool --version'] run_guarded.sh <logfile> <command> [args...]
#
# * <logfile> follows logs/<rule>/<gene>.log; rule and gene are recorded in
#   the structured records below (taken from that path).
# * Every attempt appends JSON-lines records (same shape as the Python
#   scripts' structured logging: ts/level/script/msg + fields): "started"
#   with the command, working dir, host and the probed tool version, and
#   "finished" with exit code and elapsed seconds. Snakemake `retries:`
#   therefore leave a complete per-attempt history in one file.
# * The command's stdout AND stderr are appended verbatim, so the Python
#   scripts' own structured records and any tool stderr land in the same log.
# * Non-zero exit: the exit code is reported, the tail of the log is dumped
#   to the console, and the same code is propagated (Snakemake marks the job
#   failed / retries).
# * PYTHONHASHSEED is pinned so hash-dependent ordering is reproducible.

set -uo pipefail

if [ "$#" -lt 2 ]; then
    echo "usage: $0 <logfile> <command> [args...]" >&2
    exit 64
fi

log="$1"; shift
mkdir -p "$(dirname "$log")"
export PYTHONHASHSEED=42

rule="$(basename "$(dirname "$log")")"
gene="$(basename "$log" .log)"

json_escape() {
    local s="$1"
    s="${s//\\/\\\\}"; s="${s//\"/\\\"}"; s="${s//$'\n'/\\n}"; s="${s//$'\t'/\\t}"
    printf '%s' "$s"
}

# record <LEVEL> <msg> [extra-json-fields-without-braces]
record() {
    printf '{"ts": "%s", "level": "%s", "script": "run_guarded", "msg": "%s", "rule": "%s", "gene": "%s"%s}\n' \
        "$(date -u +%Y-%m-%dT%H:%M:%S.000+00:00)" "$1" "$(json_escape "$2")" \
        "$(json_escape "$rule")" "$(json_escape "$gene")" "${3:+, $3}" >> "$log"
}

tool_version=""
if [ -n "${VERSION_CMD:-}" ]; then
    # Best-effort probe; a missing/odd --version must never fail the job.
    tool_version="$(timeout 20 bash -c "$VERSION_CMD" 2>&1 | head -n 3 || true)"
fi

start_s=$SECONDS
record INFO "started" \
    "\"command\": \"$(json_escape "$*")\", \"cwd\": \"$(json_escape "$PWD")\", \"host\": \"$(json_escape "$(hostname)")\", \"tool_version_cmd\": \"$(json_escape "${VERSION_CMD:-}")\", \"tool_version\": \"$(json_escape "$tool_version")\", \"pythonhashseed\": 42"

"$@" >> "$log" 2>&1
rc=$?
elapsed=$((SECONDS - start_s))

if [ "$rc" -ne 0 ]; then
    record ERROR "finished: FAILED" "\"exit_code\": $rc, \"elapsed_s\": $elapsed"
    echo "[run_guarded] ${rule}/${gene} FAILED (exit code $rc): $*" >&2
    echo "--- last 40 lines of $log ---" >&2
    tail -n 40 "$log" >&2
    exit "$rc"
fi

record INFO "finished: OK" "\"exit_code\": 0, \"elapsed_s\": $elapsed"
