#!/bin/bash -l
#SBATCH --job-name=pdbclean_quality
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --output=/users/%u/fastscratch/pdbclean_logs/quality_%A_%a.out
#SBATCH --error=/users/%u/fastscratch/pdbclean_logs/quality_%A_%a.err

set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 CONFIG_PATH MANIFEST_PATH" >&2
    exit 2
fi

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    echo "SLURM_ARRAY_TASK_ID is not set" >&2
    exit 2
fi

CONFIG_PATH="$(realpath "$1")"
MANIFEST_PATH="$(realpath "$2")"

SCRIPT_DIR="$(
    cd -- "$(dirname -- "${BASH_SOURCE[0]}")" >/dev/null 2>&1
    pwd
)"
REPOSITORY_ROOT="$(
    cd "$SCRIPT_DIR/../.." >/dev/null 2>&1
    pwd
)"

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
echo "PDBClean quality task"
echo "Job ID:       ${SLURM_JOB_ID}"
echo "Array job ID: ${SLURM_ARRAY_JOB_ID:-}"
echo "Array task:   ${SLURM_ARRAY_TASK_ID}"
echo "Node:         $(hostname)"
echo "Started:      $(date --iso-8601=seconds)"
echo "Repository:   ${REPOSITORY_ROOT}"
echo "Config:       ${CONFIG_PATH}"
echo "Manifest:     ${MANIFEST_PATH}"
echo "========================================"

python scripts/pdbclean/run_quality_task.py \
    --config "$CONFIG_PATH" \
    --manifest "$MANIFEST_PATH" \
    --task-id "$SLURM_ARRAY_TASK_ID"

echo "========================================"
echo "Completed: $(date --iso-8601=seconds)"
echo "========================================"
