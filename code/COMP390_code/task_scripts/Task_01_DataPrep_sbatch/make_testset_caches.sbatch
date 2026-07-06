#!/bin/bash -l
#SBATCH -J make_testset_caches
#SBATCH -p long
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=1-00:00:00
#SBATCH --output=/users/%u/COMP702_BeyondAF/logs/slurm/make_testset_caches_%j.log
#SBATCH --error=/users/%u/COMP702_BeyondAF/logs/slurm/make_testset_caches_%j.err

set -euo pipefail

source "$HOME/COMP702_BeyondAF/config/comp702_paths.sh"

# Activate the shared OpenFold environment.
module load miniforge3/25.3.0-python3.12.10
eval "$(conda shell.bash hook)"
conda activate "$OPENFOLD_ENV"

# Centralised paths.
MMCIF_DIR="$TESTSET_MMCIF_DIR"
CACHE_DIR="$TESTSET_CACHE_DIR"
CLUSTER_FILE="$TESTSET_CLUSTER_FILE"

MMCIF_CACHE="$TESTSET_MMCIF_CACHE"
CHAIN_CACHE="$TESTSET_CHAIN_CACHE"

MMCIF_TMP="${MMCIF_CACHE}.tmp"
CHAIN_TMP="${CHAIN_CACHE}.tmp"

mkdir -p "$CACHE_DIR"

echo "Job $SLURM_JOB_ID started on $(hostname) at $(date)"
echo "MMCIF_DIR=$MMCIF_DIR"
echo "CACHE_DIR=$CACHE_DIR"
echo "CLUSTER_FILE=$CLUSTER_FILE"

if [ ! -d "$MMCIF_DIR" ]; then
  echo "ERROR: mmCIF directory not found: $MMCIF_DIR"
  exit 1
fi

if [ ! -f "$CLUSTER_FILE" ]; then
  echo "ERROR: cluster file not found: $CLUSTER_FILE"
  exit 1
fi

echo "=== Step 1: generate mmCIF cache ==="
rm -f "$MMCIF_TMP"
python "$OPENFOLD_CODE_DIR/scripts/generate_mmcif_cache.py" \
  "$MMCIF_DIR" \
  "$MMCIF_TMP" \
  --no_workers "$SLURM_CPUS_PER_TASK"

mv "$MMCIF_TMP" "$MMCIF_CACHE"
ls -lh "$MMCIF_CACHE"

echo "=== Step 2: generate chain data cache ==="
rm -f "$CHAIN_TMP"
python "$OPENFOLD_CODE_DIR/scripts/generate_chain_data_cache.py" \
  "$MMCIF_DIR" \
  "$CHAIN_TMP" \
  --cluster_file "$CLUSTER_FILE" \
  --no_workers "$SLURM_CPUS_PER_TASK"

mv "$CHAIN_TMP" "$CHAIN_CACHE"
ls -lh "$CHAIN_CACHE"

echo "All done at $(date)"
