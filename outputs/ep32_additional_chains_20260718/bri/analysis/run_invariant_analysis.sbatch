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

module load miniforge3/25.3.0-python3.12.10
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$HOME/fastscratch/envs/bri_env_1.2.2"

python "$ANALYSIS_SCRIPT" \
  --prediction-dir "$PREDICTION_DIR" \
  --reference "$REFERENCE_CSV" \
  --mapping "$MAPPING_CSV" \
  --target-label "$TARGET_LABEL" \
  --output-dir "$OUTPUT_DIR"
