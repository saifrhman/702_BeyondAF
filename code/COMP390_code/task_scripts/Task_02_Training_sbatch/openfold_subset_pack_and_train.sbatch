#!/bin/bash -l
#SBATCH -J openfold_subset
#SBATCH -p gpu-a-lowsmall
#SBATCH -N 1
#SBATCH --gres=gpu:2
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=1-00:00:00
#SBATCH --output=/users/%u/COMP702_BeyondAF/logs/slurm/openfold_subset_%j.log
#SBATCH --error=/users/%u/COMP702_BeyondAF/logs/slurm/openfold_subset_%j.err

set -euo pipefail

source "$HOME/COMP702_BeyondAF/config/comp702_paths.sh"

echo "======================================"
echo "Job ID: ${SLURM_JOB_ID:-manual}"
echo "Node: $(hostname)"
echo "Start time: $(date)"
echo "CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-unset}"
echo "======================================"

# Activate the shared OpenFold environment.
module load miniforge3/25.3.0-python3.12.10
eval "$(conda shell.bash hook)"
conda activate "$OPENFOLD_ENV"

N_WORKERS="${SLURM_CPUS_PER_TASK:-16}"
GPUS_PER_NODE="${SLURM_GPUS_ON_NODE:-2}"

export OMP_NUM_THREADS="$N_WORKERS"
export MKL_NUM_THREADS="$N_WORKERS"
export PYTHONPATH="$OPENFOLD_CODE_DIR:${PYTHONPATH:-}"

# Source data on shared storage.
SRC_DATA="$OPENFOLD_DATA_ROOT"

# GPU-local staged data.
GPU_LOCAL_BASE="$OPENFOLD_GPU_LOCAL_BASE"
GPU_DATA="$OPENFOLD_GPU_DATA_ROOT"
TARFILE="$OPENFOLD_GPU_DATA_TAR"
OUT_ROOT="${OPENFOLD_GPU_SUBSET_OUT_ROOT:-$OPENFOLD_GPU_RUNS_ROOT/run_subset_a100}"

mkdir -p "$GPU_LOCAL_BASE"
mkdir -p "$OUT_ROOT"

echo "[CHECK] SRC_DATA=$SRC_DATA"
test -d "$SRC_DATA" || { echo "Missing: $SRC_DATA"; exit 1; }

# Stage the data onto GPU-local storage if it is not already present.
if [ ! -d "$GPU_DATA" ]; then
  echo "Packing OpenFold data from shared storage..."

  TAR_COMPRESS="gzip"
  command -v pigz >/dev/null && TAR_COMPRESS="pigz"

  tar -C "$(dirname "$SRC_DATA")" -cf - "$(basename "$SRC_DATA")" \
    | "$TAR_COMPRESS" > "$TARFILE"

  echo "Extracting OpenFold data to GPU-local storage..."
  tar -xf "$TARFILE" -C "$GPU_LOCAL_BASE"
else
  echo "GPU-local OpenFold data already exists, skipping transfer."
fi

# Use GPU-local data for training.
DATA_ROOT="$GPU_DATA"

MMCIF_DIR="$DATA_ROOT/pdb_data/mmcif/mmcif_files"
ALIGN_DB_DIR="$DATA_ROOT/alignment_subset/alignment_dbs"
ALIGN_INDEX="$ALIGN_DB_DIR/alignment_db.index"
CHAIN_CACHE="$DATA_ROOT/pdb_data/data_caches/chain_data_cache.json"
MMCIF_CACHE="$DATA_ROOT/pdb_data/data_caches/mmcif_cache.json"

echo "[CHECK] MMCIF_DIR=$MMCIF_DIR"
test -d "$MMCIF_DIR" || { echo "Missing: $MMCIF_DIR"; exit 1; }

echo "[CHECK] ALIGN_DB_DIR=$ALIGN_DB_DIR"
test -d "$ALIGN_DB_DIR" || { echo "Missing: $ALIGN_DB_DIR"; exit 1; }

echo "[CHECK] ALIGN_INDEX=$ALIGN_INDEX"
test -f "$ALIGN_INDEX" || { echo "Missing: $ALIGN_INDEX"; exit 1; }

echo "[CHECK] CHAIN_CACHE=$CHAIN_CACHE"
test -f "$CHAIN_CACHE" || { echo "Missing: $CHAIN_CACHE"; exit 1; }

echo "[CHECK] MMCIF_CACHE=$MMCIF_CACHE"
test -f "$MMCIF_CACHE" || { echo "Missing: $MMCIF_CACHE"; exit 1; }

echo "Starting OpenFold subset training..."
python -u "$OPENFOLD_CODE_DIR/train_openfold.py" \
  "$MMCIF_DIR" \
  "$ALIGN_DB_DIR" \
  "$MMCIF_DIR" \
  "$OUT_ROOT" \
  --alignment_index_path "$ALIGN_INDEX" \
  --train_chain_data_cache_path "$CHAIN_CACHE" \
  --template_release_dates_cache_path "$MMCIF_CACHE" \
  --max_template_date 2021-10-10 \
  --config_preset initial_training \
  --seed 42 \
  --num_nodes 1 \
  --gpus "$GPUS_PER_NODE" \
  --num_workers "$N_WORKERS" \
  --precision bf16-mixed

echo "======================================"
echo "Finished at: $(date)"
echo "======================================"
