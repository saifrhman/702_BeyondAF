#!/bin/bash -l
#SBATCH --job-name=pdbclean_stage
#SBATCH --nodes=1
#SBATCH --ntasks=1
#
# The generic batch wrapper for one stage of one frozen run.
#
# Partition, time, memory and CPU count are supplied by the submitter from the
# stage registry (overridable per stage under `execution.slurm` in the run
# configuration), so they are deliberately absent here.
#
#   sbatch task_scripts/run_stage.sbatch \
#       <RUN_DIR> <STAGE_ID> <COMMIT> <RESOLVED_SHA256> <SCIENTIFIC_SHA256> \
#       [--allow-dirty-worktree]
#
# This script contains no scientific value and builds no scientific argument.
# It hands the run directory and the stage id to `pdbclean.stage_runner`, which
# verifies the frozen configuration on this node and only then builds and
# executes the stage's own command.
#
# The three digests are what the *submission* believed about the run. The
# runner re-reads them from the run directory and aborts if they differ, so a
# job can never quietly execute a run other than the one that was submitted.

set -euo pipefail

if [[ $# -lt 5 ]]; then
    echo "Usage: $0 RUN_DIR STAGE_ID COMMIT RESOLVED_SHA256 SCIENTIFIC_SHA256 [--allow-dirty-worktree]" >&2
    exit 2
fi

RUN_DIR="$1"
STAGE_ID="$2"
EXPECT_COMMIT="$3"
EXPECT_RESOLVED_SHA="$4"
EXPECT_SCIENTIFIC_SHA="$5"
shift 5

ALLOW_DIRTY=()

for argument in "$@"; do
    case "$argument" in
        --allow-dirty-worktree)
            ALLOW_DIRTY+=("--allow-dirty-worktree")
            ;;
        *)
            echo "Unexpected argument: $argument" >&2
            exit 2
            ;;
    esac
done

# Slurm copies the batch script into its spool directory, so the repository
# cannot be found from BASH_SOURCE inside a running job.
if [[ -z "${PDBCLEAN_REPO_ROOT:-}" ]]; then
    PDBCLEAN_REPO_ROOT="${SLURM_SUBMIT_DIR:-$PWD}"
fi

PDBCLEAN_REPO_ROOT="$(cd "$PDBCLEAN_REPO_ROOT" && pwd)"
export PDBCLEAN_REPO_ROOT

if [[ ! -f "$PDBCLEAN_REPO_ROOT/config/pdbclean/pipeline_env.sh" ]]; then
    echo "ERROR: not a PDBClean repository root: $PDBCLEAN_REPO_ROOT" >&2
    exit 2
fi

source "$PDBCLEAN_REPO_ROOT/config/pdbclean/pipeline_env.sh"

if [[ ! -d "$RUN_DIR" ]]; then
    echo "ERROR: run directory does not exist: $RUN_DIR" >&2
    exit 2
fi

# The pinned scientific environment. `module load` is unavailable on some
# hosts (and unnecessary when the interpreter is already absolute), so its
# absence is not fatal.
if command -v module >/dev/null 2>&1; then
    module load miniforge3/25.3.0-python3.12.10 || true
fi

if [[ -x "${PDBCLEAN_CONDA_ENV:-}/bin/python" ]]; then
    PDBCLEAN_PYTHON="$PDBCLEAN_CONDA_ENV/bin/python"
fi

cd "$PDBCLEAN_REPO_ROOT"

export PYTHONUNBUFFERED=1

echo "============================================================"
echo "PDBClean stage job"
echo "Job:        ${SLURM_JOB_ID:-manual}"
echo "Host:       $(hostname)"
echo "Repository: $PDBCLEAN_REPO_ROOT"
echo "Run:        $RUN_DIR"
echo "Stage:      $STAGE_ID"
echo "Python:     $PDBCLEAN_PYTHON"
echo "Started:    $(date --iso-8601=seconds)"
echo "============================================================"

exec "$PDBCLEAN_PYTHON" -m pdbclean.stage_runner \
    --run-dir "$RUN_DIR" \
    --stage "$STAGE_ID" \
    --repo-root "$PDBCLEAN_REPO_ROOT" \
    --expect-commit "$EXPECT_COMMIT" \
    --expect-resolved-config-sha256 "$EXPECT_RESOLVED_SHA" \
    --expect-scientific-config-sha256 "$EXPECT_SCIENTIFIC_SHA" \
    ${ALLOW_DIRTY[@]+"${ALLOW_DIRTY[@]}"}
