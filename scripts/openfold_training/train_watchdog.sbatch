#!/bin/bash -l
#SBATCH --job-name=of_train_watchdog
#SBATCH --partition=cpu-l40s-low
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:10:00
#SBATCH --output=/mnt/fastscratch/users/sgsrehm1/openfold_cache/logs/watchdog_%j.out
#SBATCH --error=/mnt/fastscratch/users/sgsrehm1/openfold_cache/logs/watchdog_%j.out

# Heartbeat for the OpenFold retraining lineage.
#
# Why this exists. On 2026-08-29 the run hit its 3-day wall clock and requeued
# exactly as designed, checkpointing cleanly at step 7296. The requeued job then
# died in 89 seconds, and nothing noticed for five days -- five days of an H100
# the run had already queued for. The proximate cause was a fused-kernel build
# failure, but the expensive part was not the crash: it was the silence after
# it. Lightning's SIGUSR1 requeue only covers the wall clock. Any *other* way of
# stopping -- a hard failure, a node going down, an operator cancel, a bad
# stored script -- leaves an empty queue and no alarm.
#
# So this job asks one question on a timer: is a training job of this run either
# running or waiting to run? If not, and the run is neither finished nor stopped
# on purpose, it resubmits from the current on-disk scripts and carries on.
#
# Resubmitting from what is on disk matters. Slurm requeues the batch script as
# it was *stored at submission*, so a bug fixed afterwards never reaches a
# running lineage -- which is precisely what happened in August. The watchdog
# breaks that: within one heartbeat of a stale-script death, the run comes back
# on the current code.
#
# To stop the chain deliberately:  touch $RUN_ROOT/.stop_training
# It also stops itself once the run reaches its final epoch.

set -uo pipefail   # deliberately not -e: a failed probe must never break the chain

RUN_ROOT="${1:?run root required}"
INTERVAL="${2:-1800}"       # seconds between heartbeats
TOTAL_EPOCHS="${3:-400}"    # must match the --max_epochs the run was submitted with

REPO=/mnt/fastscratch/users/sgsrehm1/COMP702_pdbclean_pipeline
SELF="$REPO/scripts/openfold_training/train_watchdog.sbatch"
RESUBMIT="$RUN_ROOT/.resubmit_cmd"
LOG="$RUN_ROOT/watchdog.log"

say() { printf '%s  %s\n' "$(date -Is)" "$*" | tee -a "$LOG"; }

say "heartbeat job ${SLURM_JOB_ID:-?} on $(hostname)"

# -- 1. deliberate stop wins over everything else -----------------------------
if [[ -f "$RUN_ROOT/.stop_training" ]]; then
    say "stop: $RUN_ROOT/.stop_training exists; chain ends here"
    exit 0
fi

# -- 2. has the run finished? -------------------------------------------------
# epoch_summary.csv is appended once per completed epoch, epochs 0-indexed, so
# the run is done when the last row reads TOTAL_EPOCHS-1. A missing or
# unparseable file means "not finished", never "finished": guessing wrong in
# that direction would end the chain on a run still needing it.
LAST_EPOCH=-1
LAST_STEP=-1
SUMMARY="$RUN_ROOT/epoch_summary.csv"

if [[ -s "$SUMMARY" ]]; then
    read -r LAST_EPOCH LAST_STEP < <(
        tail -1 "$SUMMARY" | awk -F, 'NF>=2 && $1 ~ /^[0-9]+$/ {print $1, $2}'
    ) || true
    LAST_EPOCH="${LAST_EPOCH:--1}"
    LAST_STEP="${LAST_STEP:--1}"
fi

say "progress: epoch $LAST_EPOCH of $TOTAL_EPOCHS (final index $((TOTAL_EPOCHS - 1))), global_step $LAST_STEP"

if (( LAST_EPOCH >= TOTAL_EPOCHS - 1 )); then
    say "done: run reached epoch $LAST_EPOCH of $TOTAL_EPOCHS; chain ends here"
    exit 0
fi

# -- 2b. is there room to checkpoint? -----------------------------------------
# Fastscratch is quota'd on file *count*, not just bytes, and it ran out on
# 2026-09-03 at 614k files against a 500k soft limit -- with the grace period
# expired, every new file is refused. A run in that state cannot write its next
# checkpoint, and no amount of resubmitting will help: the rescue would fail the
# same way. Resubmitting is still worth doing (training itself continues until
# the next checkpoint), but the warning has to be loud, because this is the one
# failure the watchdog cannot repair on its own.
# Counted from the end of the row -- lfs wraps the filesystem onto its own line
# when the path is long, which shifts every field if counted from the front.
# The last four are always: files-used, files-quota, files-limit, grace.
INODES="$(lfs quota -u "$USER" /mnt/fastscratch 2>/dev/null | tr -d '*' \
          | awk 'NF>=8 && $(NF-3) ~ /^[0-9]+$/ {print $(NF-3), $(NF-2); exit}')"
if [[ -n "$INODES" ]]; then
    read -r USED SOFT <<<"$INODES"
    if [[ "$USED" =~ ^[0-9]+$ && "$SOFT" =~ ^[0-9]+$ ]] && (( SOFT > 0 )); then
        # Headroom, not percentage. Once the caches are cleared the floor is the
        # MSAs and the conda envs -- around 96% of the soft quota and not
        # reducible -- so a percentage alarm would fire on every heartbeat
        # forever. What actually matters is whether the few files a day this run
        # creates still fit, so alarm on the absolute count that remains.
        LEFT=$(( SOFT - USED ))
        say "fastscratch inodes: $USED / $SOFT ($LEFT free)"
        (( LEFT < 5000 )) && say "WARNING: under 5000 files of quota left; checkpoint writes are at risk"
    fi
fi

# -- 3. schedule the next heartbeat before doing any work ---------------------
# Ordering is deliberate. If the resubmit logic below were to die -- a Slurm
# hiccup, a full filesystem -- a chain scheduled afterwards would die with it,
# and the watchdog would fail in exactly the silent way it exists to prevent.
# Scheduled first, the chain survives any failure of this instance's own body.
#
# Everything from here needs Slurm, and Slurm is not always there: this session
# saw "slurm_load_partitions: Unable to contact slurm controller" from an
# ordinary squeue. An unreachable controller makes squeue print nothing and exit
# non-zero, which is indistinguishable from "no such job" if only the output is
# read -- and reading it that way is the worst mistake this script could make.
# It would conclude the run had died and submit a second trainer onto a live run
# root, the double-writer the claim file exists to prevent.
#
# So probe first, and treat a blip as something to wait out rather than act on.
# A chain that snaps on a transient outage is a watchdog that has stopped
# watching, so the probe retries before giving up.
controller_up() { squeue -u "$USER" -h -o '%i' >/dev/null 2>&1; }

for attempt in 1 2 3 4 5; do
    controller_up && break
    say "Slurm controller unreachable (attempt $attempt/5); waiting 30s"
    sleep 30
done

if ! controller_up; then
    say "ERROR: Slurm controller still unreachable. No successor scheduled and no"
    say "       rescue attempted -- sbatch would fail too. Resubmit by hand:"
    say "       sbatch $SELF $RUN_ROOT $INTERVAL $TOTAL_EPOCHS"
    exit 1
fi

# One chain only. Two would race to resubmit and could put two training jobs on
# the same run root, which the claim file resolves only by killing one.
OTHER_WD="$(squeue -u "$USER" --name=of_train_watchdog -h -o '%i' 2>/dev/null \
            | grep -v "^${SLURM_JOB_ID:-none}$" | head -1)"

if [[ -n "$OTHER_WD" ]]; then
    say "another heartbeat ($OTHER_WD) is already queued; not scheduling a second chain"
else
    # Retried, because this one call is the whole chain. Every later heartbeat
    # exists only if this succeeds, so a single transient refusal here would end
    # the watching silently -- the exact failure mode this script was written
    # for, reproduced in the script itself.
    NEXT=""
    for attempt in 1 2 3; do
        NEXT="$(sbatch --parsable --begin="now+${INTERVAL}" \
                    "$SELF" "$RUN_ROOT" "$INTERVAL" "$TOTAL_EPOCHS" 2>&1)"
        [[ "$NEXT" =~ ^[0-9]+$ ]] && break
        say "scheduling attempt $attempt/3 failed: $NEXT"
        sleep 30
    done

    if [[ "$NEXT" =~ ^[0-9]+$ ]]; then
        say "next heartbeat : job $NEXT in ${INTERVAL}s"
    else
        say "WARNING: could not schedule next heartbeat: $NEXT"
    fi
fi

# -- 4. is the run already covered? -------------------------------------------
# Any of_train_full job counts, pending as much as running. A pending one is the
# normal state while the GPU tier ratchet waits for faster hardware, and a job
# that has just been superseded is expected to disappear -- neither is a reason
# to submit more.
ALIVE="$(squeue -u "$USER" --name=of_train_full -h -o '%i %T' 2>/dev/null)"

if [[ -n "$ALIVE" ]]; then
    say "training covered by: $(echo "$ALIVE" | tr '\n' ';')"
    exit 0
fi

# An empty answer is the only one that triggers action, so it gets checked
# twice. The controller could have gone down in the seconds since the probe, and
# a false empty costs a duplicate trainer on a live run root -- far worse than
# waiting one more heartbeat to be sure.
sleep 10
if ! controller_up; then
    say "controller went away mid-check; treating this heartbeat as inconclusive"
    exit 0
fi
ALIVE="$(squeue -u "$USER" --name=of_train_full -h -o '%i %T' 2>/dev/null)"
if [[ -n "$ALIVE" ]]; then
    say "training covered on recheck: $(echo "$ALIVE" | tr '\n' ';')"
    exit 0
fi

say "NO training job is running or queued (confirmed twice) -- the run is stalled"

# -- 5. stall guard -----------------------------------------------------------
# A run that dies instantly every time would otherwise be resubmitted forever,
# burning queue priority and filling the logs. Progress is the test: if three
# consecutive rescues have not advanced global_step by even one, the fault is
# not transient and a human needs to look at it.
STALLS=0
[[ -f "$RUN_ROOT/.watchdog_stalls" ]] && STALLS="$(cat "$RUN_ROOT/.watchdog_stalls")"
PREV_STEP=-1
[[ -f "$RUN_ROOT/.watchdog_last_step" ]] && PREV_STEP="$(cat "$RUN_ROOT/.watchdog_last_step")"

if [[ "$LAST_STEP" == "$PREV_STEP" ]]; then
    STALLS=$((STALLS + 1))
else
    STALLS=0
fi

echo "$STALLS"    > "$RUN_ROOT/.watchdog_stalls"
echo "$LAST_STEP" > "$RUN_ROOT/.watchdog_last_step"

if (( STALLS >= 3 )); then
    say "GIVING UP: $STALLS rescues in a row made no progress past step $LAST_STEP."
    say "           Not resubmitting. Check the newest log in"
    say "           /mnt/fastscratch/users/sgsrehm1/openfold_cache/logs/ then"
    say "           rm $RUN_ROOT/.watchdog_stalls to re-arm."
    exit 0
fi

# -- 6. rescue ----------------------------------------------------------------
if [[ ! -x "$RESUBMIT" ]]; then
    say "ERROR: $RESUBMIT is missing or not executable; cannot rescue the run"
    exit 1
fi

# The dead job's claim would otherwise make the rescuer wait on a corpse.
if [[ -f "$RUN_ROOT/.active_job" ]]; then
    HOLDER="$(cat "$RUN_ROOT/.active_job")"
    if ! squeue -j "$HOLDER" -h -o '%T' 2>/dev/null | grep -q .; then
        say "clearing stale claim: job $HOLDER is gone"
        rm -f "$RUN_ROOT/.active_job" "$RUN_ROOT/.active_gpu_tier"
    fi
fi

say "resubmitting from $RESUBMIT (rescue $((STALLS + 1)) at step $LAST_STEP)"
"$RESUBMIT" 2>&1 | tee -a "$LOG"
say "rescue submitted"
