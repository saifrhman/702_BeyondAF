#!/bin/bash -l
#
# Preflight for a fresh checkout: resolve every path the pipeline will use,
# say whether it exists, and for anything missing say what to do about it.
#
# Run it first, on a login node:
#
#     bash scripts/pdbclean_doctor.sh
#
# It changes nothing. Exit status is 0 when the PDBClean pipeline can run, 1
# when something it needs is missing. The OpenFold retraining layer is
# reported separately, because the PDBClean pipeline does not depend on it.

set -uo pipefail

source "${PDBCLEAN_ENV_FILE:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/../config/pdbclean/pipeline_env.sh}"

pass=0; warn=0; fail=0

ok()   { printf "  \033[32mOK\033[0m    %-30s %s\n" "$1" "$2"; pass=$((pass+1)); }
note() { printf "  \033[33mNOTE\033[0m  %-30s %s\n" "$1" "$2"; warn=$((warn+1)); }
bad()  { printf "  \033[31mMISS\033[0m  %-30s %s\n" "$1" "$2"; fail=$((fail+1)); }

check_dir() {   # name path required(yes|no) hint
    local n="$1" p="$2" req="$3" hint="$4"
    if [[ -d "$p" ]]; then
        if [[ -w "$p" ]]; then ok "$n" "$p"; else note "$n" "$p  (not writable)"; fi
    elif [[ "$req" == yes ]]; then bad "$n" "$p"$'\n'"        -> $hint"
    else note "$n" "$p  (absent)"$'\n'"        -> $hint"; fi
}

echo "PDBClean preflight"
echo "  user         ${USER}"
echo "  repository   ${PDBCLEAN_REPO_ROOT}"
echo "  env file     ${PDBCLEAN_ENV_FILE:-${PDBCLEAN_REPO_ROOT}/config/pdbclean/pipeline_env.sh}"
echo
echo "Filesystem roots (derived from \$USER; override any in the environment)"
check_dir PDBCLEAN_SCRATCH_ROOT     "$PDBCLEAN_SCRATCH_ROOT"     yes \
  "Barkla gives every user /mnt/scratch/users/\$USER. If yours is elsewhere: export PDBCLEAN_SCRATCH_ROOT=..."
check_dir PDBCLEAN_FASTSCRATCH_ROOT "$PDBCLEAN_FASTSCRATCH_ROOT" yes \
  "Barkla gives every user /mnt/fastscratch/users/\$USER. If yours is elsewhere: export PDBCLEAN_FASTSCRATCH_ROOT=..."

echo
echo "PDBClean pipeline"
if [[ -x "$PDBCLEAN_PYTHON" ]] || command -v "$PDBCLEAN_PYTHON" >/dev/null 2>&1; then
    v="$("$PDBCLEAN_PYTHON" -c 'import sys;print(".".join(map(str,sys.version_info[:3])))' 2>/dev/null)"
    if "$PDBCLEAN_PYTHON" -c 'import pyarrow, numpy, yaml' >/dev/null 2>&1; then
        ok "python"  "$PDBCLEAN_PYTHON  (${v:-?}, pyarrow+numpy+yaml present)"
    else
        bad "python" "$PDBCLEAN_PYTHON  (${v:-?}, missing pyarrow/numpy/yaml)"$'\n'"        -> conda env create -p \$PDBCLEAN_FASTSCRATCH_ROOT/envs/bri_env_1.2.2 -f reproducibility/bri_environment.yml"$'\n'"           (or point PDBCLEAN_PYTHON at an interpreter that already has them)"
    fi
else
    bad "python" "$PDBCLEAN_PYTHON not executable"$'\n'"        -> export PDBCLEAN_PYTHON=/path/to/python, or export PDBCLEAN_CONDA_ENV=/path/to/env"
fi

if "$PDBCLEAN_PYTHON" -c "import sys;sys.path.insert(0,'$PDBCLEAN_REPO_ROOT/src');import pdbclean" >/dev/null 2>&1; then
    ok "pdbclean package" "importable from \$PDBCLEAN_REPO_ROOT/src"
else
    bad "pdbclean package" "not importable"$'\n'"        -> pip install -e .   (from the repository root)"
fi

# Output roots are created by the run itself; absent is normal on a fresh
# checkout, so report them without treating absence as a failure.
check_dir PDBCLEAN_OUTPUT_ROOT  "$PDBCLEAN_OUTPUT_ROOT"  no "created on first run; mkdir -p to pre-check writability"
check_dir PDBCLEAN_RELEASE_ROOT "$PDBCLEAN_RELEASE_ROOT" no "created on first run"
check_dir PDBCLEAN_RUN_ROOT     "$PDBCLEAN_RUN_ROOT"     no "created on first run"

if command -v sbatch >/dev/null 2>&1; then
    ok "slurm" "$(command -v sbatch)  partition ${PDBCLEAN_SLURM_PARTITION}"
else
    note "slurm" "sbatch not on PATH -- local execution only (pdbclean run --executor local)"
fi

echo
echo "OpenFold retraining layer (optional; the PDBClean pipeline does not need it)"
check_dir OPENFOLD_SRC        "$OPENFOLD_SRC"        no "git clone the OpenFold source here, or export OPENFOLD_SRC=..."
check_dir OPENFOLD_TRAIN_ENV  "$OPENFOLD_TRAIN_ENV"  no "create the training conda env here, or export OPENFOLD_TRAIN_ENV=..."
check_dir OPENFOLD_MMCIF_DIR  "$OPENFOLD_MMCIF_DIR"  no "materialise the snapshot mmCIFs here, or export OPENFOLD_MMCIF_DIR=..."
check_dir OPENFOLD_MSA_STORE  "$OPENFOLD_MSA_STORE"  no "generate MSAs here, or export OPENFOLD_MSA_STORE=..."
if [[ -x "$MMSEQS_BIN" ]]; then ok "mmseqs" "$MMSEQS_BIN  ($("$MMSEQS_BIN" version 2>/dev/null))"
else note "mmseqs" "$MMSEQS_BIN absent"$'\n'"        -> needed only by Stage 14d and MSA generation; export MMSEQS_BIN=..."; fi

echo
echo "  $pass ok, $warn note, $fail missing"
if (( fail )); then
    echo
    echo "  The PDBClean pipeline cannot run until the MISS lines are resolved."
    echo "  Every path above is overridable; nothing in this repository names a user."
    exit 1
fi
echo "  Ready. Next:  pdbclean plan --config config/pdbclean/profiles/comp702_frozen_20260101.yaml"
exit 0
