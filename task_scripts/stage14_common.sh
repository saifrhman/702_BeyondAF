# Shared preamble for the Stage-14 task scripts.
#
# Source this from a Slurm batch script:
#
#     source "$PDBCLEAN_REPO_ROOT/task_scripts/stage14_common.sh"
#
# It resolves the repository, sources the pipeline environment, selects the
# interpreter, and exposes `pdbclean_stage_command` so a wrapper never has to
# repeat a scientific argument. Every scientific value comes from the resolved
# run configuration; nothing scientific is written here.

# ---------------------------------------------------------------------------
# Repository root
# ---------------------------------------------------------------------------
# Slurm copies the batch script into its spool directory, so BASH_SOURCE cannot
# locate the repository from inside a running job. Resolve it from the
# environment or from the submission directory instead.
if [[ -z "${PDBCLEAN_REPO_ROOT:-}" ]]; then
    PDBCLEAN_REPO_ROOT="${SLURM_SUBMIT_DIR:-$PWD}"
fi

PDBCLEAN_REPO_ROOT="$(cd "$PDBCLEAN_REPO_ROOT" && pwd)"
export PDBCLEAN_REPO_ROOT

if [[ ! -f "$PDBCLEAN_REPO_ROOT/config/pdbclean/pipeline_env.sh" ]]; then
    echo "ERROR: not a PDBClean repository root: $PDBCLEAN_REPO_ROOT" >&2
    echo "Submit from the repository root, or export PDBCLEAN_REPO_ROOT." >&2
    exit 2
fi

source "$PDBCLEAN_REPO_ROOT/config/pdbclean/pipeline_env.sh"

cd "$PDBCLEAN_REPO_ROOT"

# ---------------------------------------------------------------------------
# Interpreter
# ---------------------------------------------------------------------------
PDBCLEAN_PY="${PDBCLEAN_PYTHON}"

if [[ ! -x "$PDBCLEAN_PY" ]]; then
    PDBCLEAN_PY="$(command -v "$PDBCLEAN_PYTHON" || true)"
fi

if [[ -z "$PDBCLEAN_PY" ]]; then
    echo "ERROR: no usable interpreter: ${PDBCLEAN_PYTHON}" >&2
    echo "Set PDBCLEAN_PYTHON to a Python with pdbclean importable." >&2
    exit 2
fi

export PDBCLEAN_PY

# ---------------------------------------------------------------------------
# Run configuration
# ---------------------------------------------------------------------------
# PDBCLEAN_PROFILE names the configuration profile the job runs under. It
# defaults to the frozen COMP702 20260101 profile so an unconfigured
# submission behaves exactly as the previous hard-coded wrappers did.
export PDBCLEAN_PROFILE="${PDBCLEAN_PROFILE:-config/pdbclean/profiles/comp702_frozen_20260101.yaml}"

if [[ ! -f "$PDBCLEAN_PROFILE" ]]; then
    echo "ERROR: configuration profile not found: $PDBCLEAN_PROFILE" >&2
    exit 2
fi

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# Print the banner every Stage-14 job shares.
pdbclean_banner() {
    echo "============================================================"
    echo "$1"
    echo "Job:        ${SLURM_JOB_ID:-manual}"
    echo "Host:       $(hostname)"
    echo "Repository: $PDBCLEAN_REPO_ROOT"
    echo "Git HEAD:   $(git rev-parse HEAD)"
    echo "Git dirty:  $(git status --porcelain | wc -l) path(s)"
    echo "Python:     $PDBCLEAN_PY"
    echo "Profile:    $PDBCLEAN_PROFILE"
    echo "Started:    $(date --iso-8601=seconds)"
    echo "============================================================"

    "$PDBCLEAN_PY" -m pdbclean.cli config --config "$PDBCLEAN_PROFILE"
}

# Emit the argv for one stage, one part per line, derived from the resolved
# run configuration. This is the same builder the CLI and the UI use.
pdbclean_stage_command() {
    "$PDBCLEAN_PY" -m pdbclean.cli stage-command \
        --stage "$1" \
        --config "$PDBCLEAN_PROFILE"
}

# Print the value the resolved configuration gives to one stage flag, e.g.
#     pdbclean_stage_argument redundancy_graph --output-dir
pdbclean_stage_argument() {
    local stage_id="$1"
    local flag="$2"

    local -a command
    mapfile -t command < <(pdbclean_stage_command "$stage_id")

    local index
    for index in "${!command[@]}"; do
        if [[ "${command[$index]}" == "$flag" ]]; then
            echo "${command[$((index + 1))]}"
            return 0
        fi
    done

    echo "ERROR: stage $stage_id has no $flag argument" >&2
    return 1
}

# Run one stage exactly as `pdbclean run` would run it.
pdbclean_run_stage() {
    local stage_id="$1"

    echo
    echo "===== STAGE: $stage_id ====="

    local -a command
    mapfile -t command < <(pdbclean_stage_command "$stage_id")

    if [[ ${#command[@]} -eq 0 ]]; then
        echo "ERROR: could not build a command for stage $stage_id" >&2
        return 2
    fi

    printf '%q ' "${command[@]}"
    echo

    "${command[@]}"
}
