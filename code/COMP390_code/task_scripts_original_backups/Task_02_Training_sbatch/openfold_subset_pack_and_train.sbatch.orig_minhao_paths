#!/bin/bash -l
#SBATCH -J openfold_subset
#SBATCH -p gpu-a-lowsmall
#SBATCH -N 1
#SBATCH --gres=gpu:2
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=1-00:00:00
#SBATCH -w gpu15
#SBATCH --output=/users/%u/fastscratch/openfold_subset_%j.log
#SBATCH --error=/users/%u/fastscratch/openfold_subset_%j.err

set -euo pipefail

echo "======================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Start time: $(date)"
echo "======================================"

# =============================
# 1️⃣ 环境
# =============================
source ~/.bashrc
conda activate openfold_env

export OMP_NUM_THREADS=16
export MKL_NUM_THREADS=16

# =============================
# 2️⃣ 路径定义
# =============================

# node037 上的数据源（共享文件系统可读）
SRC_DATA=/users/sgmwu14/localscratch

# GPU 节点本地 scratch（真正训练用）
GPU_LOCAL_BASE=/localscratch/$USER
GPU_DATA=$GPU_LOCAL_BASE/openfold_data

# 打包文件
TARFILE=$GPU_LOCAL_BASE/openfold_data.tar.gz

# 输出
OUT_ROOT=$GPU_LOCAL_BASE/openfold_runs/run_subset_a100

mkdir -p "$GPU_LOCAL_BASE"
mkdir -p "$OUT_ROOT"

# =============================
# 3️⃣ 打包并传输（只在 GPU 节点做）
# =============================

if [ ! -d "$GPU_DATA" ]; then
  echo "📦 Packing localscratch from node037..."

  # 使用 pigz（多核 gzip），如果没有就自动退化为 gzip
  TAR_COMPRESS="gzip"
  command -v pigz >/dev/null && TAR_COMPRESS="pigz"

  tar -C "$(dirname $SRC_DATA)" -cf - "$(basename $SRC_DATA)" \
    | $TAR_COMPRESS > "$TARFILE"

  echo "📦 Extracting to GPU local scratch..."
  mkdir -p "$GPU_DATA"
  tar -xf "$TARFILE" -C "$GPU_LOCAL_BASE"

else
  echo "✅ GPU local data already exists, skip transfer"
fi

# =============================
# 4️⃣ 训练路径（⚠️全部使用 GPU 本地）
# =============================

DATA_ROOT=$GPU_DATA

MMCIF_DIR=$DATA_ROOT/pdb_data/mmcif/mmcif_files
ALIGN_DB_DIR=$DATA_ROOT/alignment_subset/alignment_dbs
ALIGN_INDEX=$ALIGN_DB_DIR/alignment_db.index

# ✅ 使用 RODA 官方 cache
CHAIN_CACHE=$DATA_ROOT/pdb_data/data_caches/chain_data_cache.json
MMCIF_CACHE=$DATA_ROOT/pdb_data/data_caches/mmcif_cache.json

# =============================
# 5️⃣ 启动 OpenFold 训练
# =============================

echo "🚀 Starting OpenFold training..."

python /users/sgmwu14/openfold/train_openfold.py \
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
  --gpus 2 \
  --num_workers 8 \
  --precision bf16-mixed

echo "======================================"
echo "Finished at: $(date)"
echo "======================================"