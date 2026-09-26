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

REPO=${PDBCLEAN_REPO_ROOT}
# Which prepared population this run trains on. Defaults to the geometric
# (PDBClean-dedup-v1) view, so an unqualified submission is the run that has
# always been submitted here; the sequence-reduced Stage-14d view is selected by
# exporting OF_TRAIN_BASE and OF_ALIGN_INDEX, which sbatch carries into the job.
BASE="${OF_TRAIN_BASE:-$REPO/outputs/openfold_training/20260101/pdbclean-dedup-v1}"
VIEW="$BASE/training_view"
ALIGN_INDEX="${OF_ALIGN_INDEX:-$BASE/alignment_index/pdbclean_dedup_v1.index}"

# Held-out validation. Absent these, the run behaves exactly as it always has:
# no validation loader, no held-out metric, checkpoint selection by recency.
#
# The three paths are set together or not at all. Defaulting the other two off
# OF_VAL_DATA_DIR would guess at a layout; a bare `set -u` abort would say only
# "unbound variable" an hour into a GPU allocation. So say which one is missing.
VAL_ARGS=()
if [[ -n "${OF_VAL_DATA_DIR:-}" ]]; then
    for _v in OF_VAL_ALIGNMENT_DIR OF_VAL_CACHE OF_VAL_PROJECTION; do
        if [[ -z "${!_v:-}" ]]; then
            echo "OF_VAL_DATA_DIR is set but $_v is not; validation needs all three" >&2
            exit 2
        fi
    done
    unset _v

    VAL_ARGS=(
        --val_data_dir "$OF_VAL_DATA_DIR"
        --val_alignment_dir "$OF_VAL_ALIGNMENT_DIR"
        --val_mmcif_data_cache_path "$OF_VAL_CACHE"
        # Validation must read the same residues training reads. Without this
        # the featuriser falls back to the deposited polymer and peak memory
        # stops being bounded by anything the split controls.
        --val_projection_index_path "$OF_VAL_PROJECTION"
        --check_val_every_n_epoch "${OF_VAL_EVERY_N_EPOCHS:-5}"
    )
    echo "validation      : $OF_VAL_DATA_DIR every ${OF_VAL_EVERY_N_EPOCHS:-5} epochs"
else
    echo "validation      : none configured"
fi
OFSRC=${OPENFOLD_SRC}

MMCIF=${OPENFOLD_MMCIF_DIR}
MSAS=${OPENFOLD_MSA_STORE}

RUN_ROOT="${OPENFOLD_RUNS_ROOT}/$RUN_NAME"
LOGDIR="$RUN_ROOT/slurm_logs"

ENV=${OPENFOLD_TRAIN_ENV}
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

# Learning-rate schedule, rescaled to this run's budget.
#
# Upstream starts decaying at step 50,000 and steps again only every 50,000
# after that. Against a 51,200-step budget that is one 5% reduction across the
# final 1,200 steps: the run would sit at a flat 1e-3 from warmup to the end.
# AlphaFold pairs that 1e-3 with a global batch of 128; this run uses 8, so the
# plateau is noisier here than the value was ever chosen for, and the epoch
# metrics showed it -- ten-epoch means of lddt_ca swinging about +/-0.025 with
# no sustained trend from roughly epoch 100 onward, some five times the
# sampling noise of the metric itself.
#
# So anneal on this run's own scale, keeping AlphaFold's staircase shape and
# its 0.95 factor and changing only the spacing: decay begins at step 25,600 --
# epoch 200, the exact halfway point of 400 -- and steps every 1,200 steps
# (9.4 epochs) thereafter, reaching 3.24e-4, about a third of maximum, at the
# final step. That leaves 59% of the learning-rate mass still to spend, so the
# back half of training still learns rather than merely freezing in place.
#
# Two things this deliberately does NOT do. It does not touch max_lr, which
# would change the first 200 epochs retroactively in interpretation if not in
# fact. And it does not add the gradient-norm clipping at 0.1 that AlphaFold
# specifies and this fork omits entirely: that is a real and separate deviation,
# but changing both at once would make the effect of either unattributable, and
# the per-step loss shows no spike pathology to justify the urgency (median 35.0,
# max 71.3, nothing above three times the median across 190 logged steps).
#
# The schedule is NOT set here. It lives in $RUN_ROOT/lr_schedule.json, because
# it is a property of the experiment rather than of the job that happens to be
# running it -- and because exports added to this file cannot reach a lineage
# already in flight: Slurm requeues the batch script as it was stored at
# submission, so this run, submitted in September and requeueing every three
# days, would never see them. OPENFOLD_LR_* still override the file if set, but
# setting them here would let the job script silently overrule the run's own
# recorded schedule, which is the wrong way round.
#
# Absent that file and those variables, the scheduler is upstream's exactly.

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

# Tier floor: the best hardware this run has ever trained on, and from then on
# the worst it will accept. The existing ratchet only stops a slower GPU
# displacing a *live* holder -- it says nothing about the gap after one ends,
# where a slower job could pick the run up unopposed and quietly keep it there.
#
# Pinning closes that. It is a deliberate trade of throughput for consistency:
# a pinned run waits for its tier rather than training slower, so if the tier is
# busy for a week, the run waits a week. That is the intended behaviour, not an
# oversight.
#
# To unpin, lower or delete this file:  rm $RUN_ROOT/.min_gpu_tier
FLOOR_FILE="$RUN_ROOT/.min_gpu_tier"
FLOOR=0
[[ -f "$FLOOR_FILE" ]] && FLOOR="$(cat "$FLOOR_FILE" 2>/dev/null || echo 0)"
[[ "$FLOOR" =~ ^[0-9]+$ ]] || FLOOR=0

if (( MY_TIER < FLOOR )); then
    echo "STANDING DOWN: this run is pinned to tier $FLOOR and above, and this job"
    echo "               is tier $MY_TIER ($MY_GPU). Releasing the node now rather"
    echo "               than holding it idle. To unpin: rm $FLOOR_FILE"
    exit 0
fi

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

# Raise the floor only after the claim is held, never before. A job that stood
# down or died on the way here has not trained on this tier, and promising a
# tier the run has not actually reached would strand it: nothing slower would be
# allowed to pick it up, and nothing faster ever arrived.
if (( MY_TIER > FLOOR )); then
    echo "$MY_TIER" > "$FLOOR_FILE"
    echo "tier floor      : raised to $MY_TIER -- this run will not accept slower hardware again"
else
    echo "tier floor      : $FLOOR"
fi

# -- re-arm the watchdog ------------------------------------------------------
# The watchdog chain is serial: each heartbeat schedules the next, so any single
# failure to reach sbatch ends the watching for good. That happened on
# 2026-09-05 -- a Slurm controller outage outlasted the retry window, the chain
# stopped, and this run trained unobserved for 29 hours. It survived, but only
# because nothing went wrong during them.
#
# So the training job re-arms it. This runs at every start and every requeue,
# which bounds any gap in coverage to one requeue cycle instead of forever. The
# two mechanisms fail independently: the watchdog rescues dead training, and
# live training revives a dead watchdog.
#
# Only ever one chain. The watchdog's own duplicate guard would stop a second
# from scheduling successors, but it would still burn a slot and log a
# confusing line, so check here too.
WATCHDOG="$REPO/scripts/openfold_training/train_watchdog.sbatch"
if squeue -u "$USER" -h -o '%i' >/dev/null 2>&1; then
    if squeue -u "$USER" --name=of_train_watchdog -h -o '%i' 2>/dev/null | grep -q .; then
        echo "watchdog        : already armed"
    elif [[ -x "$WATCHDOG" ]]; then
        # The watchdog declares a stall when no new epoch has been recorded
        # within this many seconds, so the interval must exceed the time one
        # epoch actually takes. It was 1800 s against an epoch of roughly
        # 98 min (128 optimiser steps at ~46 s), which made a false stall
        # certain: every 30 minutes it submitted a rescue, the rescue claimed
        # the run root, the running job stood down and the in-flight epoch was
        # lost. The seqid100 lineage sat at step 199 for two and a half hours
        # that way, with eight jobs queued.
        #
        # EPOCH_LEN samples at ACCUM per step is EPOCH_LEN/(GPUS*ACCUM) steps;
        # at a conservative 60 s per step that is the expected epoch, and the
        # watchdog waits three of them before calling a stall.
        STEPS_PER_EPOCH=$(( EPOCH_LEN / (GPUS * ACCUM) ))
        STALL_SECONDS=$(( STEPS_PER_EPOCH * 60 * 3 ))
        (( STALL_SECONDS < 1800 )) && STALL_SECONDS=1800
        echo "watchdog        : stall threshold ${STALL_SECONDS}s (~3 epochs of ${STEPS_PER_EPOCH} steps)"
        WD="$(sbatch --parsable "$WATCHDOG" "$RUN_ROOT" "$STALL_SECONDS" "$TOTAL_EPOCHS" 2>&1)"
        if [[ "$WD" =~ ^[0-9]+$ ]]; then
            echo "watchdog        : re-armed as job $WD (chain was not running)"
        else
            echo "watchdog        : WARNING could not re-arm: $WD"
        fi
    fi
else
    echo "watchdog        : WARNING Slurm controller unreachable; not re-armed"
fi

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
    "$ALIGN_INDEX" \
    "$VIEW/projection_index.json" \
    "$VIEW/train_chain_data_cache.json" \
    "$VIEW/train_filter_view.txt"

echo
echo "=== training ==="

cd "$OFSRC"

# The #SBATCH header names one partition and the submitter overrides it with
# another, so the job inherits whichever memory default each one carries:
# gpu-l40s sets DefMemPerNode=UNLIMITED, gpu-h100 sets DefMemPerCPU=21150. With
# both in the environment srun refuses to start --
#
#   srun: fatal: SLURM_MEM_PER_CPU, SLURM_MEM_PER_GPU, and SLURM_MEM_PER_NODE
#         are mutually exclusive.
#
# and the job dies in about a second, which is long enough to look like a
# scheduling hiccup and short enough to leave no useful log. Clearing all three
# makes srun inherit the allocation sbatch already holds, which is what is
# wanted in every case: the memory is decided by the allocation, not by which
# partition's default happened to be read last.
unset SLURM_MEM_PER_CPU SLURM_MEM_PER_GPU SLURM_MEM_PER_NODE

set -x
srun --kill-on-bad-exit=1 "$ENV/bin/python" train_openfold.py \
    "$MMCIF" \
    "$MSAS" \
    "$MMCIF" \
    "$RUN_ROOT" \
    2026-01-01 \
    --config_preset initial_training \
    --alignment_index_path "$ALIGN_INDEX" \
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
    "${VAL_ARGS[@]}" \
    $RESUME
STATUS=$?
set +x

echo "train_openfold exit: $STATUS"
echo "finished        : $(date -Is)"

exit $STATUS
