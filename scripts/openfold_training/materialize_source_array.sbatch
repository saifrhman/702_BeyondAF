#!/bin/bash -l
#SBATCH --job-name=of26_mmcif
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=/users/%u/fastscratch/COMP702_pdbclean_pipeline/outputs/openfold_training/20260101/pdbclean-dedup-v1/source_materialisation/logs/%A_%a.out
#SBATCH --error=/users/%u/fastscratch/COMP702_pdbclean_pipeline/outputs/openfold_training/20260101/pdbclean-dedup-v1/source_materialisation/logs/%A_%a.err

set -euo pipefail

REPOSITORY_ROOT="$HOME/fastscratch/COMP702_pdbclean_pipeline"

TRAIN="$REPOSITORY_ROOT/outputs/openfold_training/20260101/pdbclean-dedup-v1"

WORKLIST="$TRAIN/source_materialisation/source_worklist.parquet"

STATUS_DIR="$TRAIN/source_materialisation/status"

PY="$HOME/fastscratch/envs/bri_env_1.2.2/bin/python"

EXPECTED_WORKLIST_SHA256="01f20d56079fd1e8c02f99dbe8d1fb0abf59fbd00596f3d5e278ab7d7f579766"

mkdir -p \
    "$TRAIN/source_materialisation/logs" \
    "$STATUS_DIR"

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    echo "SLURM_ARRAY_TASK_ID is not set" >&2
    exit 2
fi

BATCH_OFFSET="${BATCH_OFFSET:-0}"
BATCH_ID=$((SLURM_ARRAY_TASK_ID + BATCH_OFFSET))

if (( BATCH_ID < 0 || BATCH_ID > 1181 )); then
    echo "Invalid logical batch ID: $BATCH_ID" >&2
    exit 2
fi

cd "$REPOSITORY_ROOT"

export PYTHONPATH="$REPOSITORY_ROOT/src"
export PYTHONUNBUFFERED=1
export PYTHONDONTWRITEBYTECODE=1
export XDG_CACHE_HOME="$HOME/fastscratch/.cache"

echo "========================================"
echo "OpenFold 20260101 source materialisation"
echo "Job ID:       ${SLURM_JOB_ID}"
echo "Array task:   ${SLURM_ARRAY_TASK_ID}"
echo "Batch offset: ${BATCH_OFFSET}"
echo "Logical batch:${BATCH_ID}"
echo "Node:         $(hostname)"
echo "Started:      $(date --iso-8601=seconds)"
echo "Repository:   ${REPOSITORY_ROOT}"
echo "Worklist:     ${WORKLIST}"
echo "========================================"

"$PY" \
    scripts/openfold_training/materialize_source_batch.py \
    --worklist "$WORKLIST" \
    --batch-id "$BATCH_ID" \
    --status-dir "$STATUS_DIR" \
    --expected-worklist-sha256 \
        "$EXPECTED_WORKLIST_SHA256"

echo
echo "Completed: $(date --iso-8601=seconds)"
