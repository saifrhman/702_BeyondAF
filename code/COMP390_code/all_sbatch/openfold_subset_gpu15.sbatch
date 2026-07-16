#!/bin/bash -l
#SBATCH -J of_subset
#SBATCH -p gpu-a-lowsmall
#SBATCH -N 1
#SBATCH -w gpu15
#SBATCH --gres=gpu:a100:2
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=1-00:00:00
#SBATCH --output=/users/%u/fastscratch/of_subset_%j.log
#SBATCH --error=/users/%u/fastscratch/of_subset_%j.err

set -euo pipefail

set +u
source /mnt/data2/users/sgmwu14/tools/miniforge3/etc/profile.d/conda.sh
conda activate openfold_env
set -u

export PYTHONPATH=/users/sgmwu14/openfold:$PYTHONPATH
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}

DATA_DIR=/users/sgmwu14/localscratch
OUT_DIR=/users/sgmwu14/scratch/openfold_runs/run_subset_gpu15_a100

MMCIF_DIR=$DATA_DIR/pdb_data/mmcif/mmcif_files
ALIGN_DB_DIR=$DATA_DIR/alignment_subset/alignment_dbs
ALIGN_INDEX=$ALIGN_DB_DIR/alignment_db.index
CHAIN_CACHE=$DATA_DIR/pdb_data/data_caches/chain_data_cache.json
MMCIF_CACHE=$DATA_DIR/pdb_data/data_caches/mmcif_cache.json
OBSOLETE=$DATA_DIR/pdb_data/obsolete.dat

mkdir -p "$OUT_DIR"

EXTRA_OBSOLETE=()
if [ -f "$OBSOLETE" ]; then
  EXTRA_OBSOLETE+=( --obsolete_pdbs_file_path "$OBSOLETE" )
fi

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
  --gpus 2 \
  --precision bf16 \
  "${EXTRA_OBSOLETE[@]}"