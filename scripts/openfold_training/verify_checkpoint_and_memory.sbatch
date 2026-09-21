#!/bin/bash -l
#SBATCH --job-name=of_ckpt_verify
#SBATCH --partition=gpu-l40s-low
#SBATCH --gres=gpu:l40s:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=01:00:00

# Portability: every path below comes from here, derived from $USER.
# Override any of them in the environment; nothing names an individual.
source "${PDBCLEAN_ENV_FILE:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/../../config/pdbclean/pipeline_env.sh}"

set -euo pipefail

CKPT="${1:?checkpoint required}"

REPO=${PDBCLEAN_REPO_ROOT}
BASE="$REPO/outputs/openfold_training/20260101/pdbclean-dedup-v1"
VIEW="$BASE/training_view"
ENV=${OPENFOLD_TRAIN_ENV}

. "$ENV/etc/conda/activate.d/10-cache-dirs.sh"
export PYTHONPATH="$REPO/src${PYTHONPATH:+:$PYTHONPATH}"

echo "host   : $(hostname)"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader

cd ${OPENFOLD_SRC}

"$ENV/bin/python" "$REPO/scripts/openfold_training/verify_checkpoint_and_memory.py" \
    --checkpoint "$CKPT" \
    --mmcif-dir ${OPENFOLD_MMCIF_DIR} \
    --msa-store ${OPENFOLD_MSA_STORE} \
    --alignment-index "$BASE/alignment_index/pdbclean_dedup_v1.index" \
    --projection-index "$VIEW/projection_index.json" \
    --train-filter "$VIEW/train_filter_view.txt" \
    --steps 2 \
    --report "$VIEW/checkpoint_reload_and_memory.json"
