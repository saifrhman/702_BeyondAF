#!/bin/bash
#SBATCH --job-name=pdbclean_release
#SBATCH --partition=nodes
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --output=logs/stage14_release_%j.out
#SBATCH --error=logs/stage14_release_%j.err

set -euo pipefail

cd "$HOME/fastscratch/COMP702_pdbclean_pipeline"

PYTHON_ENV="${COMP702_PYTHON:-$HOME/fastscratch/envs/bri_env_1.2.2/bin/python}"

if [[ ! -x "$PYTHON_ENV" ]]; then
    echo "ERROR: Python not executable: $PYTHON_ENV"
    exit 1
fi

mapfile -t ROOTS < <(
    find outputs/pdbclean \
      -mindepth 2 \
      -maxdepth 2 \
      -type d \
      -name 'protocol3.2-comp702-v1' \
      -print
)

if [[ ${#ROOTS[@]} -ne 1 ]]; then
    echo "ERROR: expected exactly one protocol root; found ${#ROOTS[@]}"
    printf '%s\n' "${ROOTS[@]}"
    exit 1
fi

ROOT="${ROOTS[0]}"

SNAPSHOT="$(basename "$(dirname "$ROOT")")"
PROTOCOL="$(basename "$ROOT")"

RELEASE="outputs/releases/PDBClean-${SNAPSHOT}-${PROTOCOL}-dedup-v1"

echo "===== FINAL RELEASE JOB ====="
echo "Host: $(hostname)"
echo "Started: $(date --iso-8601=seconds)"
echo "Git HEAD: $(git rev-parse HEAD)"
echo "Protocol root: $ROOT"
echo "Release: $RELEASE"
echo

"$PYTHON_ENV" scripts/build_stage14_final_release.py \
  --protocol-root "$ROOT" \
  --policy-config config/pdbclean/stage14_representative_policy_v1.yaml \
  --output-dir "$RELEASE"

echo
echo "Finished: $(date --iso-8601=seconds)"
echo "FINAL RELEASE JOB: PASS"
