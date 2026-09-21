#!/bin/bash -l
#
# Submit the OpenFold retraining lineage.
#
# One run, several queues. The same job is submitted to every GPU partition the
# account can reach, and whichever starts first begins training; when a faster
# GPU frees up, that job supersedes the slower one through the claim files in
# the run root. gpu_train_full.sbatch only ever ratchets upward in HBM
# bandwidth, so a late A100 stands down rather than displacing a running H100.
#
# This script is also what the watchdog replays to rescue a stalled run, which
# is why it is a file in the repo and not a command typed once into a terminal:
# a rescue must reproduce the submission exactly, including the arguments that
# define the experiment.
#
# Usage:  submit_train_full.sh [run_name] [total_epochs] [epoch_len] [gpus] [accum]

# Portability: every path below comes from here, derived from $USER.
# Override any of them in the environment; nothing names an individual.
source "${PDBCLEAN_ENV_FILE:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/../../config/pdbclean/pipeline_env.sh}"

set -uo pipefail

RUN_NAME="${1:-pdbclean_dedup_v1_scratch}"
TOTAL_EPOCHS="${2:-400}"
EPOCH_LEN="${3:-1024}"
GPUS="${4:-1}"
ACCUM="${5:-8}"

REPO=${PDBCLEAN_REPO_ROOT}
SBATCH="$REPO/scripts/openfold_training/gpu_train_full.sbatch"
LOGS=${OPENFOLD_LOG_ROOT}
RUN_ROOT="${OPENFOLD_RUNS_ROOT}/$RUN_NAME"

mkdir -p "$LOGS"

# Already covered? Submitting on top of a live lineage would only make the
# claim logic kill one of the two, losing an in-flight epoch for nothing.
# A second lineage (a different prepared population) is a different
# experiment and must not be mistaken for this one, so the guard and the claim
# logic key on a job name the caller can set.
JOB_NAME="${OF_JOB_NAME:-of_train_full}"
EXISTING="$(squeue -u "$USER" --name="$JOB_NAME" -h -o '%i %T %P' 2>/dev/null)"
if [[ -n "$EXISTING" && "${FORCE:-0}" != "1" ]]; then
    echo "training already queued or running:"
    echo "$EXISTING" | sed 's/^/  /'
    echo "nothing submitted. Set FORCE=1 to submit anyway."
    exit 0
fi

echo "run          : $RUN_NAME"
echo "experiment   : $TOTAL_EPOCHS epochs x $EPOCH_LEN samples, ${GPUS} gpu x batch 1 x accum $ACCUM"
echo "global batch : $((GPUS * ACCUM))"
echo "steps/epoch  : $((EPOCH_LEN / (GPUS * ACCUM)))"
echo "total steps  : $((TOTAL_EPOCHS * EPOCH_LEN / (GPUS * ACCUM)))"
echo

# Once the run is pinned to a tier there is no point queueing below it: those
# jobs would allocate a GPU, read the floor, stand down and release it. Cheap,
# but it burns queue priority and puts noise in the logs for no possible gain.
FLOOR=0
[[ -f "$RUN_ROOT/.min_gpu_tier" ]] && FLOOR="$(cat "$RUN_ROOT/.min_gpu_tier" 2>/dev/null || echo 0)"
[[ "$FLOOR" =~ ^[0-9]+$ ]] || FLOOR=0
(( FLOOR > 0 )) && echo "tier floor   : $FLOOR (submitting only to partitions at or above it)"

# tag  partition          gres          walltime   tier  -- ordered fastest first
TARGETS=(
    "h100 gpu-h100        gpu:h100:$GPUS 3-00:00:00 4"
    "a100 gpu-a100-lowbig gpu:a100:$GPUS 1-00:00:00 3"
    "l40s gpu-l40s        gpu:l40s:$GPUS 3-00:00:00 2"
)

SUBMITTED=()

for target in "${TARGETS[@]}"; do
    read -r tag part gres walltime tier <<<"$target"

    if (( tier < FLOOR )); then
        echo "skipped      : $part -- tier $tier is below this run's floor of $FLOOR"
        continue
    fi

    jid="$(sbatch --parsable \
              --partition="$part" \
              --gres="$gres" \
              --time="$walltime" \
              --ntasks-per-node="$GPUS" \
              --job-name="$JOB_NAME" \
              --output="$LOGS/${tag}_%j.out" \
              --error="$LOGS/${tag}_%j.err" \
              "$SBATCH" "$RUN_NAME" "$TOTAL_EPOCHS" "$EPOCH_LEN" "$GPUS" "$ACCUM" 2>&1)"

    if [[ "$jid" =~ ^[0-9]+$ ]]; then
        echo "submitted    : $jid  $part ($tag, $walltime)"
        SUBMITTED+=("$jid")
    else
        # A partition the account cannot reach is not a failure of the run --
        # the other queues still cover it. Say so and keep going.
        echo "skipped      : $part -- $jid"
    fi
done

if (( ${#SUBMITTED[@]} == 0 )); then
    echo "FATAL: no partition accepted the job" >&2
    exit 1
fi

echo
echo "run root     : $RUN_ROOT"
echo "watch with   : squeue -u $USER --name=$JOB_NAME"
