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
# Barkla filesystem roots
# ---------------------------------------------------------------------------
# Barkla gives every user three areas under their own name. Deriving them from
# $USER is what lets a second person run this pipeline without editing a single
# script: nothing below names an individual.
#
# $HOME/scratch and $HOME/fastscratch are the conventional symlinks. Where they
# are absent the /mnt paths are used directly, so both layouts work.
#
# These come first because the sections below are defined in terms of them.
export PDBCLEAN_HOME_ROOT="${PDBCLEAN_HOME_ROOT:-$HOME}"

if [[ -z "${PDBCLEAN_SCRATCH_ROOT:-}" ]]; then
    if [[ -d "$HOME/scratch" ]]; then
        PDBCLEAN_SCRATCH_ROOT="$HOME/scratch"
    else
        PDBCLEAN_SCRATCH_ROOT="/mnt/scratch/users/${USER}"
    fi
fi
export PDBCLEAN_SCRATCH_ROOT

if [[ -z "${PDBCLEAN_FASTSCRATCH_ROOT:-}" ]]; then
    if [[ -d "$HOME/fastscratch" ]]; then
        PDBCLEAN_FASTSCRATCH_ROOT="$HOME/fastscratch"
    else
        PDBCLEAN_FASTSCRATCH_ROOT="/mnt/fastscratch/users/${USER}"
    fi
fi
export PDBCLEAN_FASTSCRATCH_ROOT

# ---------------------------------------------------------------------------
# Python interpreter and environment
# ---------------------------------------------------------------------------
# Conda environment activated by the Slurm array wrappers.
export PDBCLEAN_CONDA_ENV="${PDBCLEAN_CONDA_ENV:-$PDBCLEAN_FASTSCRATCH_ROOT/envs/bri_env_1.2.2}"

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
export PDBCLEAN_LOG_ROOT="${PDBCLEAN_LOG_ROOT:-$PDBCLEAN_FASTSCRATCH_ROOT/pdbclean_logs}"

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

# ---------------------------------------------------------------------------
# OpenFold retraining layer
# ---------------------------------------------------------------------------
# These are large, mutable working areas and belong on scratch, not in the
# repository. Each is overridable on its own, so a user who keeps one of them
# somewhere unusual does not have to move the rest.
export OPENFOLD_SRC="${OPENFOLD_SRC:-$PDBCLEAN_FASTSCRATCH_ROOT/openfold_src}"
export OPENFOLD_TRAIN_ENV="${OPENFOLD_TRAIN_ENV:-$PDBCLEAN_FASTSCRATCH_ROOT/envs/openfold_train}"
export OPENFOLD_CACHE_ROOT="${OPENFOLD_CACHE_ROOT:-$PDBCLEAN_FASTSCRATCH_ROOT/openfold_cache}"
export OPENFOLD_LOG_ROOT="${OPENFOLD_LOG_ROOT:-$OPENFOLD_CACHE_ROOT/logs}"
export OPENFOLD_RUNS_ROOT="${OPENFOLD_RUNS_ROOT:-$PDBCLEAN_FASTSCRATCH_ROOT/openfold_runs}"
export OPENFOLD_MMCIF_DIR="${OPENFOLD_MMCIF_DIR:-$PDBCLEAN_SCRATCH_ROOT/COMP702_openfold_20260101/pdbclean-dedup-v1/mmcif}"
export OPENFOLD_MSA_STORE="${OPENFOLD_MSA_STORE:-$PDBCLEAN_FASTSCRATCH_ROOT/COMP702_openfold_msa/20260101/pdbclean-dedup-v1/msas}"
export MMSEQS_BIN="${MMSEQS_BIN:-$OPENFOLD_TRAIN_ENV/bin/mmseqs}"

# Scratch space for compilers and kernel caches. Left on the default /tmp these
# fill a compute node's local disk and the job dies mid-epoch.
export TMPDIR="${TMPDIR:-$OPENFOLD_CACHE_ROOT/tmp}"
export TRITON_CACHE_DIR="${TRITON_CACHE_DIR:-$OPENFOLD_CACHE_ROOT/triton}"
export TORCH_EXTENSIONS_DIR="${TORCH_EXTENSIONS_DIR:-$OPENFOLD_CACHE_ROOT/torch_extensions}"

# ---------------------------------------------------------------------------
# COMP390 legacy tree (historical analyses only; not on the PDBClean path)
# ---------------------------------------------------------------------------
export COMP702_ROOT="${COMP702_ROOT:-$PDBCLEAN_HOME_ROOT/COMP702_BeyondAF}"
export COMP390_ROOT="${COMP390_ROOT:-$COMP702_ROOT/code/COMP390_code}"
