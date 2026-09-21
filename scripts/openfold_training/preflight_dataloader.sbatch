#!/bin/bash -l
#SBATCH --job-name=of_preflight
#SBATCH --partition=nodes
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00:00

# Portability: every path below comes from here, derived from $USER.
# Override any of them in the environment; nothing names an individual.
source "${PDBCLEAN_ENV_FILE:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/../../config/pdbclean/pipeline_env.sh}"

set -euo pipefail

REPO=${PDBCLEAN_REPO_ROOT}
BASE="$REPO/outputs/openfold_training/20260101/pdbclean-dedup-v1"
VIEW="$BASE/training_view"
ENV=${OPENFOLD_TRAIN_ENV}

. "$ENV/etc/conda/activate.d/10-cache-dirs.sh"
export PYTHONPATH="$REPO/src${PYTHONPATH:+:$PYTHONPATH}"

echo "host    : $(hostname)"
echo "started : $(date -Is)"

"$ENV/bin/python" "$REPO/scripts/openfold_training/preflight_dataloader.py" \
    --mmcif-dir ${OPENFOLD_MMCIF_DIR} \
    --msa-store ${OPENFOLD_MSA_STORE} \
    --alignment-index "$BASE/alignment_index/pdbclean_dedup_v1.index" \
    --projection-index "$VIEW/projection_index.json" \
    --train-filter "$VIEW/train_filter_view.txt" \
    --samples 12 \
    --report "$VIEW/dataloader_preflight.json"

echo "finished: $(date -Is)"
