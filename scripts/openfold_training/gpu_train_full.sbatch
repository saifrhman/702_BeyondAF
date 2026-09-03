#!/bin/bash -l
#SBATCH --job-name=of_train_full
#SBATCH --partition=gpu-l40s
#SBATCH --gres=gpu:l40s:2
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
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
SUPERSEDE_JOB="${6:-}"

# Global batch = GPUS * per-GPU batch(1) * ACCUM. Keep that product fixed when
# changing GPU count, or the run is not the same experiment.
echo "global batch    : $(( GPUS * ACCUM ))  (gpus=$GPUS x batch 1 x accum=$ACCUM)"

# Under SLURM, Lightning takes its world size from SLURM_NTASKS rather than from
# --gpus, and the training command must be launched by srun so that one task
# exists per GPU. Getting this wrong does not fail loudly: Lightning reports
# "MEMBER: 1/1", trains on a single GPU while holding the others idle, and the
# effective global batch silently becomes GPUS times smaller than intended.
# Refuse to start rather than train a different experiment than the one named.
NTASKS="${SLURM_NTASKS:-1}"

if [[ "$NTASKS" != "$GPUS" ]]; then
    echo "FATAL: SLURM_NTASKS=$NTASKS but GPUS=$GPUS." >&2
    echo "       Submit with --ntasks-per-node=$GPUS so Lightning starts one" >&2
    echo "       rank per GPU, or the run would use one GPU at global batch" >&2
    echo "       $(( GPUS * ACCUM / GPUS )) instead of $(( GPUS * ACCUM ))." >&2
    exit 2
fi

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
# torch JIT-compiles CUDA extensions by shelling out to ninja, so the env bin
# must be on PATH -- sourcing the cache-dirs hook alone does not do that, and
# the failure surfaces only as "Ninja is required to load C++ extensions".
export PATH="$ENV/bin:$PATH"

# Host C++ compiler for any JIT-built CUDA extension. Prepending $ENV/bin is
# required for ninja, but it also shadows /usr/bin/gcc with conda's -- and this
# env has gcc_impl_linux-64 without gxx_impl_linux-64, so $ENV/libexec/gcc ships
# cc1, cc1obj and cc1objplus but no cc1plus. C++ sources then die with
# "gcc: fatal error: cannot execute 'cc1plus'". The node's own gcc is complete.
export CC=/usr/bin/gcc
export CXX=/usr/bin/g++
export NVCC_PREPEND_FLAGS="-ccbin ${CXX}"

# One kernel cache per GPU architecture. torch keys its JIT extension cache by
# extension name alone, so every GPU this run touches shares one directory. The
# arch shows up in the compile flags (-gencode ... sm_80 vs sm_90, -DGPU_ARCH),
# so torch does notice and rebuild rather than loading a wrong-architecture .so
# -- correctness was never at risk. What it costs is about four minutes of
# rebuild at every handover, and the tier ratchet hands over deliberately
# (L40S -> A100 -> H100). Worse, a job building for one architecture while
# another still trains on a second writes into the same directory.
#
# Split by architecture and a handover pays the build once per GPU type ever,
# instead of once per handover, with no chance of two types colliding.
ARCH="$(nvidia-smi --query-gpu=compute_cap --format=csv,noheader 2>/dev/null | head -1 | tr -d '. ')"
if [[ ! "$ARCH" =~ ^[0-9]+$ ]]; then
    # Older drivers have no compute_cap field. The GPU name separates the
    # architectures just as well; it is only ever a cache key.
    ARCH="$(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -1 | tr -cd '[:alnum:]')"
fi
export TORCH_EXTENSIONS_DIR="${TORCH_EXTENSIONS_DIR}/sm_${ARCH:-unknown}"
mkdir -p "$TORCH_EXTENSIONS_DIR"
echo "kernel cache    : $TORCH_EXTENSIONS_DIR"

# Attention kernel. This OpenFold enables DeepSpeed's fused Evoformer attention
# only under long_sequence_inference, which asserts not-train, so training runs
# the unfused path by default. The fused kernel is the same arithmetic with a
# different reduction order, but that is verified rather than assumed: the
# benchmark writes a verdict and this job reads it at start-up. Deferring the
# decision to run time means a job can sit in the queue while the evidence is
# still being gathered, instead of forfeiting its queue position.
EXP_CONFIG="$REPO/scripts/openfold_training/gpu_full_config.json"
VERDICT="$VIEW/attention_kernel_benchmark.json"

if [[ -f "$VERDICT" ]] && grep -q '"verdict": "adopt"' "$VERDICT"; then
    export CUTLASS_PATH="$OFSRC/cutlass"

    # The kernel is JIT-built per GPU architecture -- the benchmark compiled it
    # for sm_89 on an L40S, and an H100 needs a fresh sm_90 build. Force that
    # build here, before training starts, so a compile failure costs a few
    # minutes at start-up instead of killing a run that is hours in. Falling
    # back to the unfused path is always correct: it is the same arithmetic.
    KBUILD="$LOGDIR/kernel_build_${SLURM_JOB_ID}.log"

    if "$ENV/bin/python" -c "
from deepspeed.ops.op_builder import EvoformerAttnBuilder
builder = EvoformerAttnBuilder()
assert builder.NAME == 'evoformer_attn', f'wrong builder: {builder.NAME}'
builder.load()
" >"$KBUILD" 2>&1; then
        EXP_CONFIG="$REPO/scripts/openfold_training/gpu_full_config_fused.json"
        echo "attention       : fused DeepSpeed Evoformer kernel (verdict adopt, built for this GPU)"
    else
        echo "attention       : unfused -- verdict is adopt but the kernel would not"
        echo "                  build on $(hostname); see $KBUILD"
    fi
else
    echo "attention       : unfused (no adopt verdict on file)"
fi

# Each checkpoint is ~1.5 GB, and upstream keeps every epoch's -- which would
# exceed the fastscratch quota long before the step budget is reached. Keep one
# rolling resume point, plus a sparse milestone series so the run still leaves a
# training history to analyse (400 epochs / 50 = 8 milestones, about 12 GB).
export OPENFOLD_SAVE_TOP_K=1
export OPENFOLD_MILESTONE_EVERY_N_EPOCHS=50

mkdir -p "$RUN_ROOT" "$LOGDIR"

# Handoff. Exactly one job may write this run root: two would race on the
# checkpoint file and interleave the logger. The claim is taken only now, once
# Slurm has actually granted the GPUs, so a job already training keeps making
# progress for the whole time this one sits in the queue. At most its in-flight
# epoch is lost, because the resume point below is its last checkpoint.
#
# The holder is read from the run root rather than passed in, because a job
# queued behind two successive predecessors cannot know at submission time which
# one it will actually be superseding.
CLAIM="$RUN_ROOT/.active_job"
CLAIM_TIER="$RUN_ROOT/.active_gpu_tier"

# Rank the GPU this job was actually given. Evoformer is memory-bandwidth
# bound, so the ordering follows HBM bandwidth rather than peak FLOPS:
# H100 3.35 TB/s > A100 2.04 > L40S 0.86 > A40 0.70.
gpu_tier() {
    local name
    name="$(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -1)"

    case "$name" in
        *H100*) echo 4 ;;
        *A100*) echo 3 ;;
        *L40S*) echo 2 ;;
        *A40*)  echo 1 ;;
        *)      echo 0 ;;
    esac
}

MY_TIER="$(gpu_tier)"
MY_GPU="$(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -1)"
echo "this gpu        : $MY_GPU (tier $MY_TIER)"

HOLDER=""
[[ -f "$CLAIM" ]] && HOLDER="$(cat "$CLAIM" 2>/dev/null || true)"
[[ -n "$SUPERSEDE_JOB" ]] && HOLDER="${HOLDER:-$SUPERSEDE_JOB}"

HOLDER_TIER=0
[[ -f "$CLAIM_TIER" ]] && HOLDER_TIER="$(cat "$CLAIM_TIER" 2>/dev/null || echo 0)"

# Only ever move the run to faster hardware. Several jobs are queued on
# different partitions at once so that whichever frees first can take over;
# without this check the last one to start would win regardless of speed, and
# an A100 landing after an H100 would quietly slow the run down.
if [[ -n "$HOLDER" && "$HOLDER" != "${SLURM_JOB_ID}" ]] \
   && squeue -j "$HOLDER" -h -o "%T" 2>/dev/null | grep -q . \
   && (( MY_TIER <= HOLDER_TIER )); then
    echo "STANDING DOWN: job $HOLDER already trains on tier $HOLDER_TIER hardware;" \
         "this job is tier $MY_TIER ($MY_GPU). Not superseding a faster or equal GPU."
    exit 0
fi

if [[ -n "$HOLDER" && "$HOLDER" != "${SLURM_JOB_ID}" ]]; then
    if squeue -j "$HOLDER" -h -o "%T" 2>/dev/null | grep -q .; then
        echo "superseding job : $HOLDER"
        scancel "$HOLDER" || true

        for _ in $(seq 1 60); do
            squeue -j "$HOLDER" -h -o "%T" 2>/dev/null | grep -q . || break
            sleep 5
        done

        # Let the dying ranks finish flushing before anything is read or written.
        sleep 20
        echo "superseded      : $HOLDER is gone"
    else
        echo "previous holder : $HOLDER (already finished)"
    fi
fi

echo "${SLURM_JOB_ID}" > "$CLAIM"
echo "$MY_TIER" > "$CLAIM_TIER"

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
srun --kill-on-bad-exit=1 "$ENV/bin/python" train_openfold.py \
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
    --experiment_config_json "$EXP_CONFIG" \
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
