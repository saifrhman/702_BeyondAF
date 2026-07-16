#!/bin/bash -l
#SBATCH -J of_subset_h100
#SBATCH -p gpu-h100
#SBATCH -N 1
#SBATCH -w gpu32
#SBATCH --gres=gpu:h100:1
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=3-00:00:00
#SBATCH --output=/users/%u/fastscratch/of_subset_h100_%j.log
#SBATCH --error=/users/%u/fastscratch/of_subset_h100_%j.err

set -euo pipefail

set +u
source /mnt/data2/users/sgmwu14/tools/miniforge3/etc/profile.d/conda.sh
conda activate openfold_env
set -u

# Avoid Triton cache on NFS
export TRITON_CACHE_DIR=/users/sgmwu14/fastscratch/triton_cache
mkdir -p "$TRITON_CACHE_DIR"

export PYTHONPATH=/users/sgmwu14/openfold:$PYTHONPATH
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}

DATA_DIR=/users/sgmwu14/localscratch
PDB_FULL=$DATA_DIR/pdb_data_full
OUT_DIR=/users/sgmwu14/scratch/openfold_runs/run_subset_gpu32_h100_ep10

# NEW full mmcif + caches paths
MMCIF_DIR=$PDB_FULL/mmcif_files
CHAIN_CACHE=$PDB_FULL/chain_data_cache.json
MMCIF_CACHE=$PDB_FULL/mmcif_cache.json
OBSOLETE=$PDB_FULL/obsolete.dat

# Alignment subset (unchanged)
ALIGN_DB_DIR=$DATA_DIR/alignment_subset/alignment_dbs
ALIGN_INDEX=$ALIGN_DB_DIR/alignment_db.index

mkdir -p "$OUT_DIR"

EXTRA_OBSOLETE=()
if [ -f "$OBSOLETE" ]; then
  EXTRA_OBSOLETE+=( --obsolete_pdbs_file_path "$OBSOLETE" )
fi

# Optional: compile only H100 arch
export TORCH_CUDA_ARCH_LIST="9.0"

python -u /users/sgmwu14/openfold/train_openfold.py \
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
  --gpus 1 \
  --precision bf16 \
  --max_epochs 10 \
  --log_every_n_steps 25 \
  --checkpoint_every_epoch \
  "${EXTRA_OBSOLETE[@]}"