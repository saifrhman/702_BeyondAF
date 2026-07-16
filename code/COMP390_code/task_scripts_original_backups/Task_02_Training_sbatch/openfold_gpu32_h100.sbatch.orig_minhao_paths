#!/bin/bash -l
#SBATCH -J of_h100
#SBATCH -p gpu-h100
#SBATCH -N 1
#SBATCH -w gpu32
#SBATCH --gres=gpu:h100:1
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=1-00:00:00
#SBATCH --output=/users/%u/fastscratch/of_h100_%j.log
#SBATCH --error=/users/%u/fastscratch/of_h100_%j.err

set -eo pipefail

echo "======================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Start time: $(date)"
echo "CUDA_VISIBLE_DEVICES=$CUDA_VISIBLE_DEVICES"
echo "======================================"

# ---------- Conda 环境 ----------
set +u
source /mnt/data2/users/sgmwu14/tools/miniforge3/etc/profile.d/conda.sh
conda activate openfold_env
set -u

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}

# ---------- 数据路径 ----------
DATA_DIR=/users/sgmwu14/localscratch
OUT_DIR=/users/sgmwu14/scratch/openfold_runs/run_subset_gpu32_h100

MMCIF_DIR=$DATA_DIR/pdb_data/mmcif/mmcif_files
ALIGN_DB_DIR=$DATA_DIR/alignment_subset/alignment_dbs
ALIGN_INDEX=$ALIGN_DB_DIR/alignment_db.index
CHAIN_CACHE=$DATA_DIR/pdb_data/data_caches/chain_data_cache.json
MMCIF_CACHE=$DATA_DIR/pdb_data/data_caches/mmcif_cache.json

mkdir -p "$OUT_DIR"

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
  --gpus 1 \
  --num_workers 8 \
  --precision bf16-mixed

echo "======================================"
echo "Finished at: $(date)"
echo "======================================"