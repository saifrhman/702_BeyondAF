#!/bin/bash -l
#SBATCH -J of_2olo_100pred_af
#SBATCH -p gpu-h100
#SBATCH -N 1
#SBATCH --gres=gpu:h100:1
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=3-00:00:00
#SBATCH --output=/users/%u/COMP702_BeyondAF/logs/slurm/of_2olo_100pred_af_%j.out
#SBATCH --error=/users/%u/COMP702_BeyondAF/logs/slurm/of_2olo_100pred_af_%j.err

set -euo pipefail

source "$HOME/COMP702_BeyondAF/config/comp702_paths.sh"

echo "===== Job started on $(date) ====="
echo "Hostname: $(hostname)"
echo "CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-unset}"

# Activate the shared OpenFold environment.
module load miniforge3/25.3.0-python3.12.10
eval "$(conda shell.bash hook)"
conda activate "$OPENFOLD_ENV"

export PYTHONPATH="$OPENFOLD_CODE_DIR:${PYTHONPATH:-}"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"
mkdir -p "$TRITON_CACHE_DIR"

OPENFOLD_SCRIPT="$OPENFOLD_CODE_DIR/run_pretrained_openfold.py"
FASTA_FILE="$OPENFOLD_INFERENCE_FASTA_DIR/2olo.fasta"
MMCIF_DIR="$OPENFOLD_INFERENCE_CIF_DIR"
ALIGN_ROOT="$OPENFOLD_INFERENCE_ALIGN_DIR"
JAX_PARAM="$OPENFOLD_AF_MODEL1_PARAM"
BASE_OUT="${OPENFOLD_2OLO_100PRED_AF_OUT:-$OPENFOLD_INFERENCE_OUT_ROOT/2olo_100pred_af}"
FASTA_DIR="${BASE_OUT}/fasta_input"
REPEAT_COUNT="${OPENFOLD_REPEAT_COUNT:-100}"

mkdir -p "$BASE_OUT" "$FASTA_DIR"

echo "[CHECK] OPENFOLD_SCRIPT=$OPENFOLD_SCRIPT"
test -f "$OPENFOLD_SCRIPT" || { echo "Missing: $OPENFOLD_SCRIPT"; exit 1; }

echo "[CHECK] FASTA_FILE=$FASTA_FILE"
test -f "$FASTA_FILE" || { echo "Missing: $FASTA_FILE"; exit 1; }

echo "[CHECK] MMCIF_DIR=$MMCIF_DIR"
test -d "$MMCIF_DIR" || { echo "Missing: $MMCIF_DIR"; exit 1; }

echo "[CHECK] ALIGN_ROOT=$ALIGN_ROOT"
test -d "$ALIGN_ROOT" || { echo "Missing: $ALIGN_ROOT"; exit 1; }

echo "[CHECK] JAX_PARAM=$JAX_PARAM"
test -f "$JAX_PARAM" || { echo "Missing: $JAX_PARAM"; exit 1; }

cp "$FASTA_FILE" "$FASTA_DIR/2olo.fasta"

echo "FASTA file : $FASTA_FILE"
echo "MMCIF dir  : $MMCIF_DIR"
echo "Align root : $ALIGN_ROOT"
echo "JAX params : $JAX_PARAM"
echo "Output dir : $BASE_OUT"
echo "Repeats    : $REPEAT_COUNT"

for i in $(seq 1 "$REPEAT_COUNT"); do
    REPEAT_DIR="$BASE_OUT/repeat_${i}"
    mkdir -p "$REPEAT_DIR"

    echo "===== Running repeat ${i}/${REPEAT_COUNT} ====="
    echo "Output -> $REPEAT_DIR"

    python3 "$OPENFOLD_SCRIPT"       "$FASTA_DIR"       "$MMCIF_DIR"       --output_dir "$REPEAT_DIR"       --use_precomputed_alignments "$ALIGN_ROOT"       --config_preset model_1       --jax_param_path "$JAX_PARAM"       --model_device cuda:0

    echo "===== Finished repeat ${i}/${REPEAT_COUNT} at $(date) ====="
done

echo "===== All predictions completed at $(date) ====="
