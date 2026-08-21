# Environment derivation for the COMP702 PDBClean pipeline.
#
# Source this from any pipeline shell/Slurm wrapper:
#
#     source "$(dirname "$0")/../config/pdbclean/pipeline_env.sh"
#
# Every value is overridable from the environment.  Nothing here is a
# scientific parameter: scientific values live in the resolved run
# configuration (see `pdbclean.defaults` and `docs/CONFIGURATION.md`).
#
# Defaults are chosen so that an unconfigured shell behaves exactly as the
# previous hard-coded wrappers did on Barkla.

# ---------------------------------------------------------------------------
# Repository root
# ---------------------------------------------------------------------------
# Derived from this file's own location, so a clone anywhere works without
# editing. PDBCLEAN_REPO_ROOT still wins if it is already exported.
if [[ -z "${PDBCLEAN_REPO_ROOT:-}" ]]; then
    _pdbclean_env_file="${BASH_SOURCE[0]}"
    _pdbclean_config_dir="$(cd "$(dirname "$_pdbclean_env_file")" && pwd)"
    PDBCLEAN_REPO_ROOT="$(cd "$_pdbclean_config_dir/../.." && pwd)"
    unset _pdbclean_env_file _pdbclean_config_dir
fi

export PDBCLEAN_REPO_ROOT

# ---------------------------------------------------------------------------
# Python interpreter and environment
# ---------------------------------------------------------------------------
# Conda environment activated by the Slurm array wrappers.
export PDBCLEAN_CONDA_ENV="${PDBCLEAN_CONDA_ENV:-$HOME/fastscratch/envs/bri_env_1.2.2}"

# PDBCLEAN_PYTHON is the interpreter used by pipeline task scripts.
#
# A bare `python` resolves to /usr/bin/python on a Barkla compute node, which
# has none of the scientific dependencies, so prefer the pinned environment's
# interpreter when it exists.  This is an interpreter-location decision, not a
# scientific one: the environment is pinned in `reproducibility/`, and a caller
# who exports PDBCLEAN_PYTHON still wins.
if [[ -z "${PDBCLEAN_PYTHON:-}" ]]; then
    if [[ -x "$PDBCLEAN_CONDA_ENV/bin/python" ]]; then
        PDBCLEAN_PYTHON="$PDBCLEAN_CONDA_ENV/bin/python"
    else
        PDBCLEAN_PYTHON="python"
    fi
fi

export PDBCLEAN_PYTHON

# ---------------------------------------------------------------------------
# Output and log roots
# ---------------------------------------------------------------------------
export PDBCLEAN_OUTPUT_ROOT="${PDBCLEAN_OUTPUT_ROOT:-$PDBCLEAN_REPO_ROOT/outputs/pdbclean}"
export PDBCLEAN_RELEASE_ROOT="${PDBCLEAN_RELEASE_ROOT:-$PDBCLEAN_REPO_ROOT/outputs/releases}"
export PDBCLEAN_RUN_ROOT="${PDBCLEAN_RUN_ROOT:-$PDBCLEAN_REPO_ROOT/outputs/runs}"
export PDBCLEAN_LOG_ROOT="${PDBCLEAN_LOG_ROOT:-$HOME/fastscratch/pdbclean_logs}"

# ---------------------------------------------------------------------------
# Configuration defaults
# ---------------------------------------------------------------------------
export PDBCLEAN_PROTOCOL_CONFIG="${PDBCLEAN_PROTOCOL_CONFIG:-$PDBCLEAN_REPO_ROOT/config/pdbclean/protocol_3_2_comp702_v1.yaml}"
export PDBCLEAN_REPRESENTATIVE_POLICY="${PDBCLEAN_REPRESENTATIVE_POLICY:-$PDBCLEAN_REPO_ROOT/config/pdbclean/stage14_representative_policy_v1.yaml}"

# ---------------------------------------------------------------------------
# Slurm defaults
# ---------------------------------------------------------------------------
# Barkla limits the number of submitted jobs, so arrays stay small and a
# physical worker strides over logical tasks.  These are execution parameters
# only; they never change a scientific result.
export PDBCLEAN_SLURM_PARTITION="${PDBCLEAN_SLURM_PARTITION:-nodes}"
export PDBCLEAN_ARRAY_WORKERS="${PDBCLEAN_ARRAY_WORKERS:-64}"
export PDBCLEAN_ARRAY_CONCURRENCY="${PDBCLEAN_ARRAY_CONCURRENCY:-4}"

export PYTHONPATH="${PDBCLEAN_REPO_ROOT}/src${PYTHONPATH:+:$PYTHONPATH}"
