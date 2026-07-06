#!/bin/bash -l
#SBATCH -J of_subset_h100_ep10
#SBATCH -p gpu-h100
#SBATCH -N 1
#SBATCH --gres=gpu:h100:1
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=3-00:00:00
#SBATCH --output=/users/%u/COMP702_BeyondAF/logs/slurm/of_subset_h100_ep10_%j.log
#SBATCH --error=/users/%u/COMP702_BeyondAF/logs/slurm/of_subset_h100_ep10_%j.err

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
GPUS_PER_NODE="${SLURM_GPUS_ON_NODE:-1}"

# Keep Triton cache away from networked home directories.
mkdir -p "$TRITON_CACHE_DIR"

export PYTHONPATH="$OPENFOLD_CODE_DIR:${PYTHONPATH:-}"
export OMP_NUM_THREADS="$N_WORKERS"
export MKL_NUM_THREADS="$N_WORKERS"
export TORCH_CUDA_ARCH_LIST="9.0"

# Full PDB data layout.
PDB_FULL="$PDB_DATA_FULL_ROOT"
MMCIF_DIR="$PDB_FULL/mmcif_files"
CHAIN_CACHE="$PDB_FULL/chain_data_cache.json"
MMCIF_CACHE="$PDB_FULL/mmcif_cache.json"
OBSOLETE="$PDB_FULL/obsolete.dat"

# Alignment subset layout.
ALIGN_DB_DIR="$SUBSET_ALIGN_DB_DIR"
ALIGN_INDEX="$SUBSET_ALIGN_INDEX"

OUT_DIR="${OPENFOLD_SUBSET_H100_EP10_OUT_DIR:-$OPENFOLD_RUNS_ROOT/run_subset_h100_ep10}"

mkdir -p "$OUT_DIR"

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

echo "Starting OpenFold H100 training up to epoch 10..."
python -u "$OPENFOLD_CODE_DIR/train_openfold.py" \
  "$MMCIF_DIR" \
  "$ALIGN_DB_DIR" \
  "$MMCIF_DIR" \
  "$OUT_DIR" \
  2021-10-10 \
  --alignment_index_path "$ALIGN_INDEX" \
  --train_chain_data_cache_path "$CHAIN_CACHE" \
  --template_release_dates_cache_path "$MMCIF_CACHE" \
  --config_preset initial_training \
  --seed 42 \
  --num_nodes 1 \
  --gpus "$GPUS_PER_NODE" \
  --precision bf16 \
  --max_epochs 10 \
  --log_every_n_steps 25 \
  --checkpoint_every_epoch \
  "${EXTRA_OBSOLETE[@]}"

echo "======================================"
echo "Finished at: $(date)"
echo "======================================"
