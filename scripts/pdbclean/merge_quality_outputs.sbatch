#!/bin/bash -l
#SBATCH --job-name=pdbclean_qmerge
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --output=/users/%u/fastscratch/pdbclean_logs/merge_%j.out
#SBATCH --error=/users/%u/fastscratch/pdbclean_logs/merge_%j.err

set -euo pipefail

if [[ $# -ne 3 ]]; then
    echo "Usage: $0 CONFIG_PATH MANIFEST_PATH REPOSITORY_ROOT" >&2
    exit 2
fi

CONFIG_PATH="$(realpath "$1")"
MANIFEST_PATH="$(realpath "$2")"
REPOSITORY_ROOT="$(realpath "$3")"

if [[ ! -d "$REPOSITORY_ROOT" ]]; then
    echo "Repository root does not exist: $REPOSITORY_ROOT" >&2
    exit 2
fi

if [[ ! -f "$REPOSITORY_ROOT/scripts/pdbclean/merge_quality_outputs.py" ]]; then
    echo "PDBClean merge runner not found under repository root: $REPOSITORY_ROOT" >&2
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
conda activate "${PDBCLEAN_CONDA_ENV:-$HOME/fastscratch/envs/bri_env_1.2.2}"

cd "$REPOSITORY_ROOT"

export PYTHONUNBUFFERED=1

echo "========================================"
echo "PDBClean quality merge"
echo "Job ID:     ${SLURM_JOB_ID}"
echo "Node:       $(hostname)"
echo "Started:    $(date --iso-8601=seconds)"
echo "Repository: ${REPOSITORY_ROOT}"
echo "Config:     ${CONFIG_PATH}"
echo "Manifest:   ${MANIFEST_PATH}"
echo "========================================"

python scripts/pdbclean/merge_quality_outputs.py \
    --config "$CONFIG_PATH" \
    --manifest "$MANIFEST_PATH"

echo "========================================"
echo "Completed: $(date --iso-8601=seconds)"
echo "========================================"
