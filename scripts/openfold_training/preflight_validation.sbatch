#!/bin/bash -l
#SBATCH --job-name=of_val_preflight
#SBATCH --partition=gpu-l40s
#SBATCH --gres=gpu:l40s:1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00:00
#SBATCH --output=val_preflight_%j.out
#SBATCH --error=val_preflight_%j.err

# Walk every validation chain through the real eval pipeline.
#
# Cheap insurance. Five separate validation-path faults in this project were
# each found by losing about two hours of H100 time; all five were data or
# wiring faults a dataloader pass would have caught. Run this after any change
# to the validation set, the projection index, or the cache -- it asks for one
# L40S rather than the H100 the run is queued for, and it trains nothing.

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

REPO=${PDBCLEAN_REPO_ROOT}
BASE="${OF_TRAIN_BASE:-$REPO/outputs/openfold_training/20260101/pdbclean-final-v1}"
ENV=${OPENFOLD_TRAIN_ENV}
OFSRC=${OPENFOLD_SRC}

. "$ENV/etc/conda/activate.d/10-cache-dirs.sh"
export PYTHONPATH="$OFSRC:$REPO/src${PYTHONPATH:+:$PYTHONPATH}"
export PATH="$ENV/bin:$PATH"

cd "$OFSRC"

echo "population   : $BASE"
echo "cache        : $BASE/val_mmcif_data_cache.json"
echo "projection   : $BASE/val_projection_index.json"
echo

"$ENV/bin/python" "$REPO/scripts/openfold_training/preflight_validation_dataloader.py" \
    --val-data-dir        "$BASE/val_data_dir" \
    --val-alignment-dir   "$BASE/val_alignment_dir" \
    --val-cache           "$BASE/val_mmcif_data_cache.json" \
    --val-projection-index "$BASE/val_projection_index.json" \
    --template-mmcif-dir  "${OPENFOLD_MMCIF_DIR}" \
    --report              "$BASE/val_preflight_report.json" \
    "${@}"
