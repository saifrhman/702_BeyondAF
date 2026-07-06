#!/bin/bash -l
#SBATCH -J of_infer3
#SBATCH -p gpu-h100
#SBATCH -N 1
#SBATCH --gres=gpu:h100:1
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=3-00:00:00
#SBATCH --output=/users/%u/COMP702_BeyondAF/logs/slurm/of_infer3_%j.log
#SBATCH --error=/users/%u/COMP702_BeyondAF/logs/slurm/of_infer3_%j.err

set -euo pipefail

source "$HOME/COMP702_BeyondAF/config/comp702_paths.sh"

echo "======================================"
echo "Job ID: ${SLURM_JOB_ID:-manual}"
echo "Node: $(hostname)"
echo "Start time: $(date)"
echo "CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-unset}"
echo "======================================"

# Activate the shared OpenFold environment.
module load miniforge3/25.3.0-python3.12.10
eval "$(conda shell.bash hook)"
conda activate "$OPENFOLD_ENV"

export PYTHONPATH="$OPENFOLD_CODE_DIR:${PYTHONPATH:-}"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"
mkdir -p "$TRITON_CACHE_DIR"

OPENFOLD_SCRIPT="$OPENFOLD_CODE_DIR/run_pretrained_openfold.py"

FASTA_DIR="$OPENFOLD_INFERENCE_FASTA_DIR"
CIF_DIR="$OPENFOLD_INFERENCE_CIF_DIR"
ALIGN_DIR="$OPENFOLD_INFERENCE_ALIGN_DIR"

AF_PARAM="$OPENFOLD_AF_MODEL1_PARAM"
OF_CKPT="$OPENFOLD_OFFICIAL_FINETUNING_CKPT"
MY_CKPT="$OPENFOLD_EP32_CKPT"

AF_OUT="${OPENFOLD_AF_OUT:-$OPENFOLD_INFERENCE_OUT_ROOT/alphafold_params_output}"
OF_OUT="${OPENFOLD_OF_OUT:-$OPENFOLD_INFERENCE_OUT_ROOT/openfold_official_output}"
MY_OUT="${OPENFOLD_MY_OUT:-$OPENFOLD_INFERENCE_OUT_ROOT/trained_model_ep32_output}"

mkdir -p "$AF_OUT" "$OF_OUT" "$MY_OUT"

echo "[CHECK] OPENFOLD_SCRIPT=$OPENFOLD_SCRIPT"
test -f "$OPENFOLD_SCRIPT" || { echo "Missing: $OPENFOLD_SCRIPT"; exit 1; }

echo "[CHECK] FASTA_DIR=$FASTA_DIR"
test -d "$FASTA_DIR" || { echo "Missing: $FASTA_DIR"; exit 1; }

echo "[CHECK] CIF_DIR=$CIF_DIR"
test -d "$CIF_DIR" || { echo "Missing: $CIF_DIR"; exit 1; }

echo "[CHECK] ALIGN_DIR=$ALIGN_DIR"
test -d "$ALIGN_DIR" || { echo "Missing: $ALIGN_DIR"; exit 1; }

echo "[CHECK] AF_PARAM=$AF_PARAM"
test -f "$AF_PARAM" || { echo "Missing: $AF_PARAM"; exit 1; }

echo "[CHECK] OF_CKPT=$OF_CKPT"
test -f "$OF_CKPT" || { echo "Missing: $OF_CKPT"; exit 1; }

echo "[CHECK] MY_CKPT=$MY_CKPT"
test -f "$MY_CKPT" || { echo "Missing: $MY_CKPT"; exit 1; }

echo "===== Running AlphaFold parameter model ====="
python3 "$OPENFOLD_SCRIPT" \
  "$FASTA_DIR" \
  "$CIF_DIR" \
  --output_dir "$AF_OUT" \
  --use_precomputed_alignments "$ALIGN_DIR" \
  --config_preset model_1 \
  --jax_param_path "$AF_PARAM" \
  --model_device cuda:0

echo "===== Running official OpenFold model ====="
python3 "$OPENFOLD_SCRIPT" \
  "$FASTA_DIR" \
  "$CIF_DIR" \
  --output_dir "$OF_OUT" \
  --use_precomputed_alignments "$ALIGN_DIR" \
  --config_preset model_1 \
  --openfold_checkpoint_path "$OF_CKPT" \
  --model_device cuda:0

echo "===== Running archived trained model checkpoint ====="
python3 "$OPENFOLD_SCRIPT" \
  "$FASTA_DIR" \
  "$CIF_DIR" \
  --output_dir "$MY_OUT" \
  --use_precomputed_alignments "$ALIGN_DIR" \
  --config_preset model_1 \
  --openfold_checkpoint_path "$MY_CKPT" \
  --model_device cuda:0

echo "======================================"
echo "All inference finished at: $(date)"
echo "======================================"
