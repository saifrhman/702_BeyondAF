#!/bin/bash -l
#SBATCH --job-name=of_view_val
#SBATCH --partition=nodes
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=08:00:00

# Portability: every path below comes from here, derived from $USER.
# Override any of them in the environment; nothing names an individual.
# Slurm copies this script into its spool directory before running it, so
# ${BASH_SOURCE[0]} does not point into the repository under sbatch. Try the
# explicit override, the repo root the submitter exported, the directory the
# job was submitted from, and only then this file's own location.
for _pdbclean_env in "${PDBCLEAN_ENV_FILE:-}" \
                     "${PDBCLEAN_REPO_ROOT:-}/config/pdbclean/pipeline_env.sh" \
                     "${SLURM_SUBMIT_DIR:-}/config/pdbclean/pipeline_env.sh" \
                     "$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd)/../../config/pdbclean/pipeline_env.sh"; do
    [[ -f "$_pdbclean_env" ]] && { source "$_pdbclean_env"; break; }
done
unset _pdbclean_env

if [[ -z "${PDBCLEAN_REPO_ROOT:-}" ]]; then
    echo "cannot locate config/pdbclean/pipeline_env.sh; set PDBCLEAN_ENV_FILE" >&2
    exit 1
fi

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
