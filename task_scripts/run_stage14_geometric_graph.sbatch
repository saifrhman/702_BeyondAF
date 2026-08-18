#!/bin/bash
#SBATCH --job-name=stage14_graph
#SBATCH --partition=nodes
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=logs/stage14_graph_%j.out
#SBATCH --error=logs/stage14_graph_%j.err

set -euo pipefail

cd "$HOME/fastscratch/COMP702_pdbclean_pipeline"

PYTHON_ENV="${COMP702_PYTHON:-$HOME/fastscratch/envs/bri_env_1.2.2/bin/python}"

if [[ ! -x "$PYTHON_ENV" ]]; then
    echo "ERROR: Python not executable: $PYTHON_ENV"
    exit 1
fi

ROOT="$(
    find outputs/pdbclean \
      -path '*/protocol3.2-comp702-v1/full_bri_nn' \
      -type d \
      -print -quit
)"

if [[ -z "$ROOT" ]]; then
    echo "ERROR: frozen full_bri_nn publication not found"
    exit 1
fi

OUT="${ROOT%/full_bri_nn}/stage14_geometric_graph"

if [[ -e "$OUT/_SUCCESS" ]]; then
    echo "ERROR: finalized Stage-14 graph already exists: $OUT"
    exit 1
fi

echo "===== STAGE-14 GRAPH JOB ====="
echo "Host: $(hostname)"
echo "Started: $(date --iso-8601=seconds)"
echo "Git HEAD: $(git rev-parse HEAD)"
echo "Source: $ROOT"
echo "Output: $OUT"
echo

"$PYTHON_ENV" scripts/build_stage14_geometric_graph.py \
  --edges "$ROOT/finalized/candidate_near_duplicates.parquet" \
  --m1-edges "$ROOT/finalized/m1_near_duplicates.parquet" \
  --output-dir "$OUT" \
  --expected-edges 1068256 \
  --expected-m1-edges 4495 \
  --expected-mge2-nodes 577760

echo
echo "Finished: $(date --iso-8601=seconds)"
echo "STAGE-14 GRAPH JOB: PASS"
