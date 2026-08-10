#!/bin/bash -l

set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 CONFIG_PATH MANIFEST_PATH" >&2
    exit 2
fi

CONFIG_PATH="$(realpath "$1")"
MANIFEST_PATH="$(realpath "$2")"

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
REPOSITORY_ROOT="$(cd "$SCRIPT_DIR/../.." >/dev/null 2>&1 && pwd)"

if [[ ! -f "$CONFIG_PATH" ]]; then
    echo "Configuration does not exist: $CONFIG_PATH" >&2
    exit 2
fi

if [[ ! -f "$MANIFEST_PATH" ]]; then
    echo "Manifest does not exist: $MANIFEST_PATH" >&2
    exit 2
fi

cd "$REPOSITORY_ROOT"

METADATA="$(
python - "$CONFIG_PATH" "$MANIFEST_PATH" <<'PY'
import sys

import pyarrow.parquet as pq

from pdbclean.config import load_config
from pdbclean.manifest import (
    manifest_partition_count,
    resolve_manifest_snapshot,
    validate_manifest_table,
)

config_path = sys.argv[1]
manifest_path = sys.argv[2]

loaded = load_config(config_path)
config = loaded.data

snapshot_config = config["snapshot"]
execution_config = config["execution"]

manifest = pq.read_table(manifest_path)

snapshot = resolve_manifest_snapshot(
    manifest,
    snapshot_config,
)

expected_count = None
expected_total_bytes = None

if snapshot_config["mode"] == "fixed":
    expected_count = snapshot_config.get(
        "expected_mmcif_count"
    )
    expected_total_bytes = snapshot_config.get(
        "expected_total_bytes"
    )

summary = validate_manifest_table(
    manifest,
    expected_snapshot=snapshot,
    expected_count=expected_count,
    expected_total_bytes=expected_total_bytes,
)

batch_size = execution_config["batch_size"]

task_count = manifest_partition_count(
    summary.row_count,
    batch_size,
)

print(
    f"{snapshot}\t"
    f"{summary.row_count}\t"
    f"{batch_size}\t"
    f"{task_count}"
)
PY
)"

IFS=$'\t' read -r SNAPSHOT MANIFEST_ROWS BATCH_SIZE TASK_COUNT <<< "$METADATA"

if [[ ! "$TASK_COUNT" =~ ^[0-9]+$ ]] || (( TASK_COUNT < 1 )); then
    echo "Invalid derived task count: $TASK_COUNT" >&2
    exit 2
fi

LAST_TASK=$((TASK_COUNT - 1))

# Production entrypoints require a clean repository. Checking here prevents
# submitting an array whose workers would all fail provenance validation.
GIT_COMMIT="$(
python - <<'PY'
from pathlib import Path

from pdbclean.provenance import resolve_clean_git_commit

print(resolve_clean_git_commit(Path.cwd()))
PY
)"

LOG_ROOT="$HOME/fastscratch/pdbclean_logs"
mkdir -p "$LOG_ROOT"

QUALITY_SCRIPT="$SCRIPT_DIR/run_quality_array.sbatch"
MERGE_SCRIPT="$SCRIPT_DIR/merge_quality_outputs.sbatch"

echo "========================================"
echo "PDBClean quality pipeline submission"
echo "Snapshot:       $SNAPSHOT"
echo "Manifest rows:  $MANIFEST_ROWS"
echo "Batch size:     $BATCH_SIZE"
echo "Task count:     $TASK_COUNT"
echo "Array range:    0-$LAST_TASK"
echo "Git commit:     $GIT_COMMIT"
echo "Config:         $CONFIG_PATH"
echo "Manifest:       $MANIFEST_PATH"
echo "Logs:           $LOG_ROOT"
echo "========================================"

ARRAY_RESULT="$(sbatch --parsable --array="0-${LAST_TASK}" "$QUALITY_SCRIPT" "$CONFIG_PATH" "$MANIFEST_PATH")"
ARRAY_JOB_ID="${ARRAY_RESULT%%;*}"

if [[ ! "$ARRAY_JOB_ID" =~ ^[0-9]+$ ]]; then
    echo "Unexpected sbatch array response: $ARRAY_RESULT" >&2
    exit 1
fi

echo "Quality array submitted: $ARRAY_JOB_ID"

MERGE_RESULT="$(sbatch --parsable --dependency="afterok:${ARRAY_JOB_ID}" "$MERGE_SCRIPT" "$CONFIG_PATH" "$MANIFEST_PATH")"
MERGE_JOB_ID="${MERGE_RESULT%%;*}"

if [[ ! "$MERGE_JOB_ID" =~ ^[0-9]+$ ]]; then
    echo "Unexpected sbatch merge response: $MERGE_RESULT" >&2
    exit 1
fi

echo "Quality merge submitted: $MERGE_JOB_ID"
echo "Dependency: afterok:$ARRAY_JOB_ID"
echo
echo "Submission complete."
