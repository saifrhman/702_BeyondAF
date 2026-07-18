#!/bin/bash
#SBATCH --job-name=ep32_bri
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00

set -euo pipefail

source "$HOME/COMP702_BeyondAF/config/comp702_paths.sh"

PDB_ID="${PDB_ID:?PDB_ID is required}"
RUN_ROOT="${RUN_ROOT:?RUN_ROOT is required}"
REFERENCE_PDB="${REFERENCE_PDB:?REFERENCE_PDB is required}"
BRI_OUT_ROOT="${BRI_OUT_ROOT:?BRI_OUT_ROOT is required}"

PDB_ID="$(
    printf '%s' "$PDB_ID" |
    tr '[:upper:]' '[:lower:]'
)"

module load miniforge3/25.3.0-python3.12.10
eval "$(conda shell.bash hook)"
conda activate "$HOME/fastscratch/envs/bri_env_1.2.2"

BRI_SCRIPT="$COMP702_ROOT/scripts/compute_bri_from_pdb_1_2_2.py"

PREDICTION_OUT="$BRI_OUT_ROOT/${PDB_ID}_ep32_relaxed_100"
REFERENCE_OUT="$BRI_OUT_ROOT/experimental_references"

mkdir -p "$PREDICTION_OUT" "$REFERENCE_OUT"

test -f "$BRI_SCRIPT"
test -d "$RUN_ROOT"
test -s "$REFERENCE_PDB"

echo "Target: $PDB_ID"
echo "Run root: $RUN_ROOT"
echo "Prediction output: $PREDICTION_OUT"
echo "Reference: $REFERENCE_PDB"

python "$BRI_SCRIPT" \
    --pdb "$REFERENCE_PDB" \
    --output "$REFERENCE_OUT/${PDB_ID}_A_experimental.csv"

for REP in $(seq 1 100); do
    INPUT_PDB="$RUN_ROOT/repeat_${REP}/predictions/${PDB_ID}_A_model_1_relaxed.pdb"
    OUTPUT_CSV="$PREDICTION_OUT/repeat_$(printf '%03d' "$REP").csv"

    test -s "$INPUT_PDB" || {
        echo "ERROR: missing prediction $INPUT_PDB"
        exit 1
    }

    echo "[$(date)] ${PDB_ID}: repeat ${REP}/100"

    python "$BRI_SCRIPT" \
        --pdb "$INPUT_PDB" \
        --output "$OUTPUT_CSV"
done

PREDICTION_COUNT=$(
    find "$PREDICTION_OUT" \
        -maxdepth 1 \
        -type f \
        -name 'repeat_*.csv' |
    wc -l
)

REFERENCE_COUNT=$(
    find "$REFERENCE_OUT" \
        -maxdepth 1 \
        -type f \
        -name "${PDB_ID}_A_experimental.csv" |
    wc -l
)

echo "Prediction CSVs: $PREDICTION_COUNT"
echo "Reference CSVs:  $REFERENCE_COUNT"

if [ "$PREDICTION_COUNT" -ne 100 ]; then
    echo "ERROR: expected 100 prediction CSVs"
    exit 1
fi

if [ "$REFERENCE_COUNT" -ne 1 ]; then
    echo "ERROR: expected one reference CSV"
    exit 1
fi

echo "[$(date)] Completed successfully"
