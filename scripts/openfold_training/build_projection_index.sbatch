#!/bin/bash -l
#SBATCH --job-name=of_proj_index
#SBATCH --partition=nodes
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00

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

REPO=${PDBCLEAN_REPO_ROOT}
BASE="$REPO/outputs/openfold_training/20260101/pdbclean-dedup-v1"
OUT="$BASE/training_view"

export TMPDIR=${OPENFOLD_CACHE_ROOT}/tmp
mkdir -p "$TMPDIR" "$OUT"

PYTHON=${PDBCLEAN_PYTHON}

echo "host        : $(hostname)"
echo "started     : $(date -Is)"
echo "pdbclean git: $(git -C "$REPO" rev-parse HEAD)"

"$PYTHON" "$REPO/scripts/openfold_training/build_projection_index.py" \
    --training-chains "$BASE/openfold_training_chains.parquet" \
    --output "$OUT/projection_index.json" \
    --manifest "$OUT/projection_index_manifest.json" \
    --expected-chains 499770

echo "finished    : $(date -Is)"
