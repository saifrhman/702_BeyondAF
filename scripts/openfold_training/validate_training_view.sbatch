#!/bin/bash -l
#SBATCH --job-name=of_view_val
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
FEATURE_CHECKS="${2:-25}"

SHARD_ID="${SLURM_ARRAY_TASK_ID:-0}"

REPO=${PDBCLEAN_REPO_ROOT}
BASE="$REPO/outputs/openfold_training/20260101/pdbclean-dedup-v1"
VIEW="$BASE/training_view"
REPORTS="$VIEW/validation/${SHARD_COUNT}"

ENV=${OPENFOLD_TRAIN_ENV}
. "$ENV/etc/conda/activate.d/10-cache-dirs.sh"

mkdir -p "$REPORTS"

echo "host       : $(hostname)"
echo "shard      : $SHARD_ID of $SHARD_COUNT"
echo "started    : $(date -Is)"

"$ENV/bin/python" "$REPO/scripts/openfold_training/validate_training_view.py" \
    --training-chains "$BASE/openfold_training_chains.parquet" \
    --projection-index "$VIEW/projection_index.json" \
    --mmcif-dir ${OPENFOLD_MMCIF_DIR} \
    --alignment-index "$BASE/alignment_index/pdbclean_dedup_v1.index" \
    --msa-store ${OPENFOLD_MSA_STORE} \
    --shard-id "$SHARD_ID" \
    --shard-count "$SHARD_COUNT" \
    --feature-checks "$FEATURE_CHECKS" \
    --report "$REPORTS/shard_$(printf '%05d' "$SHARD_ID").json"

echo "finished   : $(date -Is)"
