#!/bin/bash -l
#
# Move the OpenFold run's write paths off /mnt/fastscratch and onto /mnt/scratch,
# leaving symlinks behind.
#
# Why symlinks rather than new paths. Slurm requeues a batch script as it was
# stored at submission, so editing gpu_train_full.sbatch cannot redirect a
# lineage already in flight -- job 10413722 was submitted in September and has
# requeued repeatedly since. A symlink at the old location needs no script to
# change its mind: every hard-coded path in the stored script keeps resolving,
# and lands on the new filesystem.
#
# Why this waits for a gap. The run root is written by a live trainer. Renaming
# a directory under a process that holds descriptors into it strands whatever it
# still has open, and a checkpoint caught mid-write would be torn. Requeue leaves
# a gap of roughly half an hour between one allocation ending and the next
# starting; the rename and symlink take under a second, so the whole swap fits
# inside it many times over.
#
# Nothing is deleted. The originals are renamed to *.relocated and left in place
# for a human to remove once the resumed run is confirmed healthy.

# Portability: every path below comes from here, derived from $USER.
# Override any of them in the environment; nothing names an individual.
source "${PDBCLEAN_ENV_FILE:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/../../config/pdbclean/pipeline_env.sh}"

set -uo pipefail

SRC=${PDBCLEAN_FASTSCRATCH_ROOT}
DST=${PDBCLEAN_SCRATCH_ROOT}
DIRS=(openfold_runs openfold_cache)
LOG="$DST/relocate.log"

say() { printf '%s  %s\n' "$(date -Is)" "$*" | tee -a "$LOG"; }

training_running() {
    squeue -u "$USER" --name=of_train_full -h -o '%T' 2>/dev/null | grep -q RUNNING
}

controller_up() { squeue -u "$USER" -h -o '%i' >/dev/null 2>&1; }

say "waiting for the requeue gap (training must not be RUNNING)"

# An unreachable controller reports nothing, which looks exactly like "no job is
# running" -- the one reading that would make this swap unsafe. Treat silence as
# "keep waiting", never as permission to proceed.
while true; do
    if ! controller_up; then
        say "slurm controller unreachable; holding"
    elif ! training_running; then
        say "no RUNNING trainer seen; confirming"
        sleep 20
        if controller_up && ! training_running; then
            break
        fi
        say "trainer reappeared or controller went away; back to waiting"
    fi
    sleep 60
done

say "gap detected -- starting swap"

for d in "${DIRS[@]}"; do
    [[ -L "$SRC/$d" ]] && { say "$d is already a symlink; skipping"; continue; }
    say "delta sync: $d"
    rsync -a --delete "$SRC/$d/" "$DST/$d/" 2>&1 | tail -2 | tee -a "$LOG"
done

# Last check before the irreversible part. The delta sync above can take a
# minute on a fresh 1.5 GB checkpoint, which is long enough for the next
# allocation to start; swapping underneath it would be the one way this script
# could do harm.
if training_running; then
    say "ABORT: a trainer started during the sync. Nothing was renamed."
    say "       The pre-seeded copy under $DST is still valid; rerun this later."
    exit 1
fi

for d in "${DIRS[@]}"; do
    [[ -L "$SRC/$d" ]] && continue
    if mv "$SRC/$d" "$SRC/$d.relocated"; then
        if ln -s "$DST/$d" "$SRC/$d"; then
            say "swapped: $SRC/$d -> $DST/$d  (original kept as $d.relocated)"
        else
            mv "$SRC/$d.relocated" "$SRC/$d"
            say "ERROR: symlink failed for $d; original restored"
            exit 1
        fi
    else
        say "ERROR: could not rename $d; left untouched"
        exit 1
    fi
done

say "verifying"
RR="$SRC/openfold_runs/pdbclean_dedup_v1_scratch"
for check in \
    "$SRC/openfold_runs" "$SRC/openfold_cache" \
    "$RR/epoch_summary.csv" "$RR/lr_schedule.json" "$SRC/openfold_cache/logs"; do
    if [[ -e "$check" ]]; then
        say "  ok: $check -> $(readlink -f "$check")"
    else
        say "  MISSING: $check"
    fi
done

say "writability test through the symlink"
T="$RR/.relocate_write_test"
if echo ok > "$T" 2>/dev/null; then
    say "  write OK, landed on $(df -h --output=target "$T" 2>/dev/null | tail -1)"
    rm -f "$T"
else
    say "  WRITE FAILED -- investigate before the next requeue"
fi

say "done. Originals preserved as $SRC/openfold_runs.relocated and openfold_cache.relocated"
