#!/bin/bash -l
#SBATCH --job-name=of_chain_cache
#SBATCH --partition=nodes
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=08:00:00

# Portability: every path below comes from here, derived from $USER.
# Override any of them in the environment; nothing names an individual.
source "${PDBCLEAN_ENV_FILE:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/../../config/pdbclean/pipeline_env.sh}"

set -euo pipefail

SHARD_COUNT="${1:?shard count required}"
SHARD_ID="${SLURM_ARRAY_TASK_ID:-0}"

REPO=${PDBCLEAN_REPO_ROOT}
BASE="$REPO/outputs/openfold_training/20260101/pdbclean-dedup-v1"
OUT="$BASE/training_view/chain_cache_shards"

ENV=${OPENFOLD_TRAIN_ENV}
. "$ENV/etc/conda/activate.d/10-cache-dirs.sh"

export PYTHONPATH="$REPO/src${PYTHONPATH:+:$PYTHONPATH}"

mkdir -p "$OUT"

echo "host    : $(hostname)"
echo "shard   : $SHARD_ID of $SHARD_COUNT"
echo "started : $(date -Is)"

"$ENV/bin/python" "$REPO/scripts/openfold_training/build_chain_data_cache.py" \
    --training-chains "$BASE/openfold_training_chains.parquet" \
    --mmcif-dir ${OPENFOLD_MMCIF_DIR} \
    --shard-id "$SHARD_ID" \
    --shard-count "$SHARD_COUNT" \
    --output "$OUT/shard_$(printf '%05d' "$SHARD_ID").json"

echo "finished: $(date -Is)"
