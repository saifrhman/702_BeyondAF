#!/bin/bash
#SBATCH --job-name=inv_analysis
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00

set -euo pipefail

TARGET_LABEL="${TARGET_LABEL:?TARGET_LABEL is required}"
PREDICTION_DIR="${PREDICTION_DIR:?PREDICTION_DIR is required}"
REFERENCE_CSV="${REFERENCE_CSV:?REFERENCE_CSV is required}"
MAPPING_CSV="${MAPPING_CSV:?MAPPING_CSV is required}"
OUTPUT_DIR="${OUTPUT_DIR:?OUTPUT_DIR is required}"
ANALYSIS_SCRIPT="${ANALYSIS_SCRIPT:?ANALYSIS_SCRIPT is required}"

PYTHON="$HOME/fastscratch/envs/bri_env_1.2.2/bin/python"

echo "===== Runtime environment ====="
echo "Host: $(hostname)"
echo "Python executable: $PYTHON"
"$PYTHON" --version

echo
echo "===== Dependency check ====="
"$PYTHON" - <<'PY'
import sys
import matplotlib
import numpy
import pandas

print("Python:", sys.executable)
print("Matplotlib:", matplotlib.__version__)
print("NumPy:", numpy.__version__)
print("Pandas:", pandas.__version__)
print("Dependencies: OK")
PY

echo
echo "===== Starting invariant analysis ====="

"$PYTHON" "$ANALYSIS_SCRIPT" \
  --prediction-dir "$PREDICTION_DIR" \
  --reference "$REFERENCE_CSV" \
  --mapping "$MAPPING_CSV" \
  --target-label "$TARGET_LABEL" \
  --output-dir "$OUTPUT_DIR"
