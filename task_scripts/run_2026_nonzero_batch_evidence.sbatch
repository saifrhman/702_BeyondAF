#!/bin/bash
#SBATCH --job-name=nonzero_bri_2026
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --output=logs/nonzero_bri_2026_%j.out
#SBATCH --error=logs/nonzero_bri_2026_%j.err

set -euo pipefail

cd "$HOME/fastscratch/COMP702_pdbclean_pipeline"

PYTHON_ENV="${COMP702_PYTHON:-$HOME/fastscratch/envs/bri_env_1.2.2/bin/python}"

if [[ ! -x "$PYTHON_ENV" ]]; then
    echo "ERROR: environment Python not executable: $PYTHON_ENV"
    exit 1
fi

CONFIG="config/pdbclean/protocol_3_2_comp702_v1.yaml"
REVIEW="reports/acta_2026_detailed_review.csv"

STAGE11="$(
    find outputs/pdbclean \
        -path '*/protocol3.2-comp702-v1/acta_downstream_investigation_v2/finalized/detailed_inspection_candidates.parquet' \
        -print -quit
)"

if [[ -z "$STAGE11" ]]; then
    echo "ERROR: Stage-11 detailed inspection candidates not found"
    exit 1
fi

SNAPSHOT_ROOT="${STAGE11%%/protocol3.2-comp702-v1/*}"
SOURCE_MANIFEST="$SNAPSHOT_ROOT/bronze/source_manifest.parquet"

if [[ ! -f "$SOURCE_MANIFEST" ]]; then
    echo "ERROR: source manifest not found: $SOURCE_MANIFEST"
    exit 1
fi

OUT="reports/acta_2026_nonzero_batch_evidence"

echo "===== 2026 NONZERO BATCH JOB ====="
echo "Host: $(hostname)"
echo "Started: $(date --iso-8601=seconds)"
echo "Git HEAD: $(git rev-parse HEAD)"
echo "Script SHA256: $(sha256sum scripts/analyze_2026_nonzero_near_duplicates.py | awk '{print $1}')"
echo "Stage11: $STAGE11"
echo "Source manifest: $SOURCE_MANIFEST"
echo "Output: $OUT"
echo

"$PYTHON_ENV" scripts/analyze_2026_nonzero_near_duplicates.py \
    --config "$CONFIG" \
    --review-csv "$REVIEW" \
    --stage11-parquet "$STAGE11" \
    --source-manifest "$SOURCE_MANIFEST" \
    --output-dir "$OUT" \
    --expected-pairs 1753 \
    --expected-chain-pairs 3179

echo
echo "Finished: $(date --iso-8601=seconds)"
echo "2026 NONZERO BATCH JOB: PASS"
