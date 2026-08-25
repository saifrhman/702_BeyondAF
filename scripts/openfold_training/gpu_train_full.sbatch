#!/bin/bash -l
#SBATCH --job-name=of_train_full
#SBATCH --partition=gpu-l40s
#SBATCH --gres=gpu:l40s:2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=3-00:00:00
#SBATCH --signal=SIGUSR1@180
#SBATCH --requeue

# Full OpenFold retraining on the PDBClean 2026-01-01 geometry-cleaned dataset.
#
# The model is initialised RANDOMLY. No pretrained AlphaFold or OpenFold weights
# are ever loaded: --resume_from_ckpt is passed only to continue *this* run from
# its own latest checkpoint across Slurm wall-time boundaries. The scientific
# goal is a new model trained on the cleaned dataset, so initialising from
# released weights would defeat the experiment entirely.
#
# Barkla caps this partition at 3 days, which is shorter than the training
# budget, so the run continues across jobs. Continuation uses Lightning's own
# SLURM handling: --signal=SIGUSR1@180 gives it three minutes' warning before
# the wall clock expires, at which point it checkpoints and requeues this same
# job id. That path restores the global step, the optimizer and the LR
# scheduler from the checkpoint, so the step counter is never reset and the
# schedule never restarts. The explicit --resume_from_ckpt below is the
# fallback for a *fresh* submission after a hard failure, where there is no
# requeue to inherit from.

set -euo pipefail

RUN_NAME="${1:?run name required}"
TOTAL_EPOCHS="${2:?total epochs required}"
EPOCH_LEN="${3:?train_epoch_len required}"
GPUS="${4:-2}"
ACCUM="${5:-4}"

# Global batch = GPUS * per-GPU batch(1) * ACCUM. Keep that product fixed when
# changing GPU count, or the run is not the same experiment.
echo "global batch    : $(( GPUS * ACCUM ))  (gpus=$GPUS x batch 1 x accum=$ACCUM)"

REPO=/mnt/fastscratch/users/sgsrehm1/COMP702_pdbclean_pipeline
BASE="$REPO/outputs/openfold_training/20260101/pdbclean-dedup-v1"
VIEW="$BASE/training_view"
OFSRC=/mnt/fastscratch/users/sgsrehm1/openfold_src

MMCIF=/users/sgsrehm1/scratch/COMP702_openfold_20260101/pdbclean-dedup-v1/mmcif
MSAS=/mnt/fastscratch/users/sgsrehm1/COMP702_openfold_msa/20260101/pdbclean-dedup-v1/msas

RUN_ROOT="/mnt/fastscratch/users/sgsrehm1/openfold_runs/$RUN_NAME"
LOGDIR="$RUN_ROOT/slurm_logs"

ENV=/mnt/fastscratch/users/sgsrehm1/envs/openfold_train
. "$ENV/etc/conda/activate.d/10-cache-dirs.sh"

export PYTHONPATH="$REPO/src${PYTHONPATH:+:$PYTHONPATH}"

# Each checkpoint is ~1.5 GB, and upstream keeps every epoch's -- which would
# exceed the fastscratch quota long before the step budget is reached. Keep one
# rolling resume point, plus a sparse milestone series so the run still leaves a
# training history to analyse (400 epochs / 50 = 8 milestones, about 12 GB).
export OPENFOLD_SAVE_TOP_K=1
export OPENFOLD_MILESTONE_EVERY_N_EPOCHS=50

mkdir -p "$RUN_ROOT" "$LOGDIR"

echo "=== identity ==="
echo "host            : $(hostname)"
echo "slurm job       : ${SLURM_JOB_ID}"
echo "run root        : $RUN_ROOT"
echo "started         : $(date -Is)"
echo "pdbclean commit : $(git -C "$REPO" rev-parse HEAD)"
echo "openfold source : $OFSRC (from da89cd28446abcd7be95459b7b349dedaee666c0)"

nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader

# -- resume from this run's own newest checkpoint, never from released weights
RESUME=""
LATEST="$(find "$RUN_ROOT" -name '*.ckpt' -printf '%T@ %p\n' 2>/dev/null \
          | sort -rn | head -1 | cut -d' ' -f2- || true)"

if [[ -n "$LATEST" ]]; then
    echo "resuming from   : $LATEST"
    RESUME="--resume_from_ckpt $LATEST"
else
    echo "resuming from   : (none -- fresh random initialisation)"
fi

echo
echo "=== dataset identity ==="
sha256sum \
    "$BASE/alignment_index/pdbclean_dedup_v1.index" \
    "$VIEW/projection_index.json" \
    "$VIEW/train_chain_data_cache.json" \
    "$VIEW/train_filter_view.txt"

echo
echo "=== training ==="

cd "$OFSRC"

set -x
"$ENV/bin/python" train_openfold.py \
    "$MMCIF" \
    "$MSAS" \
    "$MMCIF" \
    "$RUN_ROOT" \
    2026-01-01 \
    --config_preset initial_training \
    --alignment_index_path "$BASE/alignment_index/pdbclean_dedup_v1.index" \
    --train_projection_index_path "$VIEW/projection_index.json" \
    --train_chain_data_cache_path "$VIEW/train_chain_data_cache.json" \
    --train_filter_path "$VIEW/train_filter_view.txt" \
    --experiment_config_json "$REPO/scripts/openfold_training/gpu_full_config.json" \
    --train_epoch_len "$EPOCH_LEN" \
    --max_epochs "$TOTAL_EPOCHS" \
    --gpus "$GPUS" \
    --num_nodes 1 \
    --precision bf16 \
    --accumulate_grad_batches "$ACCUM" \
    --seed 20260101 \
    --log_every_n_steps 50 \
    --log_lr \
    --num_sanity_val_steps 0 \
    --checkpoint_every_epoch \
    $RESUME
STATUS=$?
set +x

echo "train_openfold exit: $STATUS"
echo "finished        : $(date -Is)"

exit $STATUS
