#!/bin/bash -l
#SBATCH -J make_testset_caches
#SBATCH -p long
#SBATCH -N 1
#SBATCH -w node037
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=1-00:00:00
#SBATCH --output=/users/%u/fastscratch/make_testset_caches_%j.log
#SBATCH --error=/users/%u/fastscratch/make_testset_caches_%j.err

set -euo pipefail

# ---- conda ----
set +u
source /mnt/data2/users/sgmwu14/tools/miniforge3/etc/profile.d/conda.sh
conda activate openfold_env
set -u

# ---- paths ----
MMCIF_DIR=/users/sgmwu14/localscratch/pdb_data/mmcif/mmcif_files
CACHE_DIR=/users/sgmwu14/localscratch/pdb_data/data_caches
CLUSTER_FILE=/users/sgmwu14/localscratch/alignment_testset/all-seqs_clusters-40.txt

MMCIF_CACHE=${CACHE_DIR}/mmcif_cache_testset.json
CHAIN_CACHE=${CACHE_DIR}/chain_data_cache_testset.json

MMCIF_TMP=${MMCIF_CACHE}.tmp
CHAIN_TMP=${CHAIN_CACHE}.tmp

mkdir -p "$CACHE_DIR"

echo "=== Step 1: generate mmcif cache ==="
rm -f "$MMCIF_TMP"
python /users/sgmwu14/openfold/scripts/generate_mmcif_cache.py \
  "$MMCIF_DIR" \
  "$MMCIF_TMP" \
  --no_workers 16
mv "$MMCIF_TMP" "$MMCIF_CACHE"
ls -lh "$MMCIF_CACHE"

echo "=== Step 2: generate chain data cache ==="
rm -f "$CHAIN_TMP"
python /users/sgmwu14/openfold/scripts/generate_chain_data_cache.py \
  "$MMCIF_DIR" \
  "$CHAIN_TMP" \
  --cluster_file "$CLUSTER_FILE" \
  --no_workers 16
mv "$CHAIN_TMP" "$CHAIN_CACHE"
ls -lh "$CHAIN_CACHE"

echo "All done at $(date)"