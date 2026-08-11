#!/bin/bash -l
#SBATCH --job-name=pdbclean_quality
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=/users/%u/fastscratch/pdbclean_logs/quality_%A_%a.out
#SBATCH --error=/users/%u/fastscratch/pdbclean_logs/quality_%A_%a.err

set -euo pipefail

if [[ $# -ne 5 ]]; then
    echo "Usage: $0 CONFIG_PATH MANIFEST_PATH LOGICAL_TASK_COUNT PHYSICAL_WORKER_COUNT REPOSITORY_ROOT" >&2
    exit 2
fi

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    echo "SLURM_ARRAY_TASK_ID is not set" >&2
    exit 2
fi

CONFIG_PATH="$(realpath "$1")"
MANIFEST_PATH="$(realpath "$2")"
LOGICAL_TASK_COUNT="$3"
PHYSICAL_WORKER_COUNT="$4"
REPOSITORY_ROOT="$(realpath "$5")"

if [[ ! "$SLURM_ARRAY_TASK_ID" =~ ^[0-9]+$ ]]; then
    echo "Invalid SLURM_ARRAY_TASK_ID: $SLURM_ARRAY_TASK_ID" >&2
    exit 2
fi

if [[ ! "$LOGICAL_TASK_COUNT" =~ ^[0-9]+$ ]] || (( LOGICAL_TASK_COUNT < 1 )); then
    echo "Invalid logical task count: $LOGICAL_TASK_COUNT" >&2
    exit 2
fi

if [[ ! "$PHYSICAL_WORKER_COUNT" =~ ^[0-9]+$ ]] || (( PHYSICAL_WORKER_COUNT < 1 )); then
    echo "Invalid physical worker count: $PHYSICAL_WORKER_COUNT" >&2
    exit 2
fi

if (( SLURM_ARRAY_TASK_ID >= PHYSICAL_WORKER_COUNT )); then
    echo "Physical worker ID is outside worker range: $SLURM_ARRAY_TASK_ID" >&2
    exit 2
fi

CONFIG_PATH="$(realpath "$CONFIG_PATH")"
MANIFEST_PATH="$(realpath "$MANIFEST_PATH")"

if [[ ! -d "$REPOSITORY_ROOT" ]]; then
    echo "Repository root does not exist: $REPOSITORY_ROOT" >&2
    exit 2
fi

if [[ ! -f "$REPOSITORY_ROOT/scripts/pdbclean/run_quality_task.py" ]]; then
    echo "PDBClean task runner not found under repository root: $REPOSITORY_ROOT" >&2
    exit 2
fi

if [[ ! -f "$CONFIG_PATH" ]]; then
    echo "Configuration does not exist: $CONFIG_PATH" >&2
    exit 2
fi

if [[ ! -f "$MANIFEST_PATH" ]]; then
    echo "Manifest does not exist: $MANIFEST_PATH" >&2
    exit 2
fi

module load miniforge3/25.3.0-python3.12.10
eval "$(conda shell.bash hook)"
conda activate "$HOME/fastscratch/envs/bri_env_1.2.2"

cd "$REPOSITORY_ROOT"

export PYTHONUNBUFFERED=1

echo "========================================"
echo "PDBClean quality physical worker"
echo "Job ID:               ${SLURM_JOB_ID}"
echo "Array job ID:         ${SLURM_ARRAY_JOB_ID:-}"
echo "Physical worker ID:   ${SLURM_ARRAY_TASK_ID}"
echo "Physical worker count:${PHYSICAL_WORKER_COUNT}"
echo "Logical task count:   ${LOGICAL_TASK_COUNT}"
echo "Node:                 $(hostname)"
echo "Started:              $(date --iso-8601=seconds)"
echo "Repository:           ${REPOSITORY_ROOT}"
echo "Config:               ${CONFIG_PATH}"
echo "Manifest:             ${MANIFEST_PATH}"
echo "========================================"

PROCESSED_COUNT=0

for (( LOGICAL_TASK_ID=SLURM_ARRAY_TASK_ID; \
       LOGICAL_TASK_ID<LOGICAL_TASK_COUNT; \
       LOGICAL_TASK_ID+=PHYSICAL_WORKER_COUNT )); do

    echo
    echo "----------------------------------------"
    echo "Physical worker ${SLURM_ARRAY_TASK_ID}: logical task ${LOGICAL_TASK_ID}"
    echo "Started logical task: $(date --iso-8601=seconds)"
    echo "----------------------------------------"

    python scripts/pdbclean/run_quality_task.py \
        --config "$CONFIG_PATH" \
        --manifest "$MANIFEST_PATH" \
        --task-id "$LOGICAL_TASK_ID"

    PROCESSED_COUNT=$((PROCESSED_COUNT + 1))

    echo "Completed logical task ${LOGICAL_TASK_ID}: $(date --iso-8601=seconds)"
done

echo
echo "========================================"
echo "Physical worker completed"
echo "Logical tasks processed: $PROCESSED_COUNT"
echo "Completed: $(date --iso-8601=seconds)"
echo "========================================"
