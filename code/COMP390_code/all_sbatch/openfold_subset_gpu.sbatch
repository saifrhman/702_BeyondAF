#!/bin/bash -l
#SBATCH -J of_subset
#SBATCH -p gpu-a-lowsmall
#SBATCH -N 1
#SBATCH --gres=gpu:a100:2
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=1-00:00:00
#SBATCH --output=/users/%u/fastscratch/of_subset_%j.log
#SBATCH --error=/users/%u/fastscratch/of_subset_%j.err

set -euo pipefail

echo "======================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Start time: $(date)"
echo "CUDA_VISIBLE_DEVICES=$CUDA_VISIBLE_DEVICES"
echo "======================================"

# ---- 环境（不要 source ~/.bashrc，避免 CONDA_PREFIX 未定义触发 -u）----
source /mnt/data2/users/sgmwu14/tools/miniforge3/etc/profile.d/conda.sh
conda activate openfold_env

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}

# ---- 路径 ----
DATA_DIR=/users/sgmwu14/localscratch
OUT_DIR=/users/sgmwu14/scratch/openfold_runs/run_subset_gpu13_a100

MMCIF_DIR=$DATA_DIR/pdb_data/mmcif/mmcif_files
ALIGN_DB_DIR=$DATA_DIR/alignment_subset/alignment_dbs
ALIGN_INDEX=$ALIGN_DB_DIR/alignment_db.index

CHAIN_CACHE=$DATA_DIR/pdb_data/data_caches/chain_data_cache.json
MMCIF_CACHE=$DATA_DIR/pdb_data/data_caches/mmcif_cache.json

OBSOLETE=$DATA_DIR/pdb_data/obsolete.dat

mkdir -p "$OUT_DIR"

# ---- 路径存在性检查（避免跑半天才报错）----
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
  echo "[CHECK] OBSOLETE not found, skip"
fi

# ---- 开始训练（stdout 会写入 log，方便 tail -f）----
echo "🚀 Starting OpenFold training..."
python -u /users/sgmwu14/openfold/train_openfold.py \
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
  --gpus 2 \
  --num_workers 8 \
  --precision bf16-mixed \
  "${EXTRA_OBSOLETE[@]}"

echo "======================================"
echo "Finished at: $(date)"
echo "======================================"