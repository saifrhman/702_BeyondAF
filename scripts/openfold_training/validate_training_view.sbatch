#!/bin/bash -l
#SBATCH --job-name=of_view_val
#SBATCH --partition=nodes
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=08:00:00

set -euo pipefail

SHARD_COUNT="${1:?shard count required}"
FEATURE_CHECKS="${2:-25}"

SHARD_ID="${SLURM_ARRAY_TASK_ID:-0}"

REPO=/mnt/fastscratch/users/sgsrehm1/COMP702_pdbclean_pipeline
BASE="$REPO/outputs/openfold_training/20260101/pdbclean-dedup-v1"
VIEW="$BASE/training_view"
REPORTS="$VIEW/validation/${SHARD_COUNT}"

ENV=/mnt/fastscratch/users/sgsrehm1/envs/openfold_train
. "$ENV/etc/conda/activate.d/10-cache-dirs.sh"

mkdir -p "$REPORTS"

echo "host       : $(hostname)"
echo "shard      : $SHARD_ID of $SHARD_COUNT"
echo "started    : $(date -Is)"

"$ENV/bin/python" "$REPO/scripts/openfold_training/validate_training_view.py" \
    --training-chains "$BASE/openfold_training_chains.parquet" \
    --projection-index "$VIEW/projection_index.json" \
    --mmcif-dir /users/sgsrehm1/scratch/COMP702_openfold_20260101/pdbclean-dedup-v1/mmcif \
    --alignment-index "$BASE/alignment_index/pdbclean_dedup_v1.index" \
    --msa-store /mnt/fastscratch/users/sgsrehm1/COMP702_openfold_msa/20260101/pdbclean-dedup-v1/msas \
    --shard-id "$SHARD_ID" \
    --shard-count "$SHARD_COUNT" \
    --feature-checks "$FEATURE_CHECKS" \
    --report "$REPORTS/shard_$(printf '%05d' "$SHARD_ID").json"

echo "finished   : $(date -Is)"
