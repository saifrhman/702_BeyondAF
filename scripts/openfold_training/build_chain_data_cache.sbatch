#!/bin/bash -l
#SBATCH --job-name=of_chain_cache
#SBATCH --partition=nodes
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=08:00:00

set -euo pipefail

SHARD_COUNT="${1:?shard count required}"
SHARD_ID="${SLURM_ARRAY_TASK_ID:-0}"

REPO=/mnt/fastscratch/users/sgsrehm1/COMP702_pdbclean_pipeline
BASE="$REPO/outputs/openfold_training/20260101/pdbclean-dedup-v1"
OUT="$BASE/training_view/chain_cache_shards"

ENV=/mnt/fastscratch/users/sgsrehm1/envs/openfold_train
. "$ENV/etc/conda/activate.d/10-cache-dirs.sh"

export PYTHONPATH="$REPO/src${PYTHONPATH:+:$PYTHONPATH}"

mkdir -p "$OUT"

echo "host    : $(hostname)"
echo "shard   : $SHARD_ID of $SHARD_COUNT"
echo "started : $(date -Is)"

"$ENV/bin/python" "$REPO/scripts/openfold_training/build_chain_data_cache.py" \
    --training-chains "$BASE/openfold_training_chains.parquet" \
    --mmcif-dir /users/sgsrehm1/scratch/COMP702_openfold_20260101/pdbclean-dedup-v1/mmcif \
    --shard-id "$SHARD_ID" \
    --shard-count "$SHARD_COUNT" \
    --output "$OUT/shard_$(printf '%05d' "$SHARD_ID").json"

echo "finished: $(date -Is)"
