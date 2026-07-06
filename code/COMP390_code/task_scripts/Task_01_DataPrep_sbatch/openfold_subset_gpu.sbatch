#!/bin/bash -l
#SBATCH -J of_subset
#SBATCH -p gpu-a-lowsmall
#SBATCH -N 1
#SBATCH --gres=gpu:a100:2
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=1-00:00:00
#SBATCH --output=/users/%u/COMP702_BeyondAF/logs/slurm/of_subset_%j.log
#SBATCH --error=/users/%u/COMP702_BeyondAF/logs/slurm/of_subset_%j.err

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

# Centralised paths.
MMCIF_DIR="$SUBSET_MMCIF_DIR"
ALIGN_DB_DIR="$SUBSET_ALIGN_DB_DIR"
ALIGN_INDEX="$SUBSET_ALIGN_INDEX"
CHAIN_CACHE="$SUBSET_CHAIN_CACHE"
MMCIF_CACHE="$SUBSET_MMCIF_CACHE"
OBSOLETE="$SUBSET_OBSOLETE_FILE"
OUT_DIR="${OPENFOLD_SUBSET_OUT_DIR:-$OPENFOLD_RUNS_ROOT/run_subset_gpu_a100}"

mkdir -p "$OUT_DIR"

# Check required input paths before starting a long training job.
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

EXTRA_OBSOLETE=()
if [ -f "$OBSOLETE" ]; then
  echo "[CHECK] OBSOLETE found: $OBSOLETE"
  EXTRA_OBSOLETE+=( --obsolete_pdbs_file_path "$OBSOLETE" )
else
  echo "[CHECK] OBSOLETE not found, skipping"
fi

echo "Starting OpenFold training..."
python -u "$OPENFOLD_CODE_DIR/train_openfold.py" \
  "$MMCIF_DIR" \
  "$ALIGN_DB_DIR" \
  "$MMCIF_DIR" \
  "$OUT_DIR" \
  --alignment_index_path "$ALIGN_INDEX" \
  --train_chain_data_cache_path "$CHAIN_CACHE" \
  --template_release_dates_cache_path "$MMCIF_CACHE" \
  --max_template_date 2021-10-10 \
  --config_preset initial_training \
  --seed 42 \
  --num_nodes 1 \
  --gpus "$GPUS_PER_NODE" \
  --num_workers "$N_WORKERS" \
  --precision bf16-mixed \
  "${EXTRA_OBSOLETE[@]}"

echo "======================================"
echo "Finished at: $(date)"
echo "======================================"
