#!/bin/bash
#SBATCH --job-name=stage14_reps
#SBATCH --partition=nodes
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --output=logs/stage14_reps_%j.out
#SBATCH --error=logs/stage14_reps_%j.err

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

GRAPH="$ROOT/stage14_geometric_graph"
FULL_BRI="$ROOT/full_bri_nn"
ACCEPTED="$ROOT/quality/merged/accepted.parquet"
META="$ROOT/downstream_metadata/finalized/entry_metadata.parquet"

OUT="$ROOT/stage14_representative_selection_v1"

if [[ -e "$OUT" ]]; then
    echo "ERROR: Stage-14 representative output already exists: $OUT"
    exit 1
fi

echo "===== STAGE-14 REPRESENTATIVE JOB ====="
echo "Host: $(hostname)"
echo "Started: $(date --iso-8601=seconds)"
echo "Git HEAD: $(git rev-parse HEAD)"
echo

"$PYTHON_ENV" scripts/select_stage14_representatives.py \
  --graph-dir "$GRAPH" \
  --edges "$FULL_BRI/finalized/candidate_near_duplicates.parquet" \
  --accepted "$ACCEPTED" \
  --metadata "$META" \
  --config config/pdbclean/stage14_representative_policy_v1.yaml \
  --output-dir "$OUT"

echo
echo "Finished: $(date --iso-8601=seconds)"
echo "STAGE-14 REPRESENTATIVE JOB: PASS"
