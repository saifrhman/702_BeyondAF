#!/bin/bash -l
#SBATCH -J smoke_2olo_ep32
#SBATCH -p gpu-h100
#SBATCH -N 1
#SBATCH --gres=gpu:h100:1
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=1-00:00:00
#SBATCH --output=/users/%u/COMP702_BeyondAF/logs/slurm/smoke_2olo_ep32_%j.out
#SBATCH --error=/users/%u/COMP702_BeyondAF/logs/slurm/smoke_2olo_ep32_%j.err

set -euo pipefail

source "$HOME/COMP702_BeyondAF/config/comp702_paths.sh"

echo "===== Job started on $(date) ====="
echo "Hostname: $(hostname)"
echo "CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-unset}"

module load miniforge3/25.3.0-python3.12.10
eval "$(conda shell.bash hook)"
conda activate "$OPENFOLD_ENV"

export PYTHONPATH="$OPENFOLD_CODE_DIR:${PYTHONPATH:-}"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"
mkdir -p "$TRITON_CACHE_DIR"

OPENFOLD_SCRIPT="$OPENFOLD_CODE_DIR/run_pretrained_openfold.py"

FASTA_DIR="$OPENFOLD_INFERENCE_FASTA_DIR"
MMCIF_DIR="$OPENFOLD_INFERENCE_CIF_DIR"
ALIGN_DIR="$OPENFOLD_INFERENCE_ALIGN_DIR"
CKPT="$OPENFOLD_EP32_CKPT"

OUT_DIR="${COMP702_ROOT}/outputs/smoke_2olo_trained_ep32_${SLURM_JOB_ID:-manual}"

mkdir -p "$OUT_DIR"

echo "[CHECK] OPENFOLD_SCRIPT=$OPENFOLD_SCRIPT"
test -f "$OPENFOLD_SCRIPT" || { echo "Missing: $OPENFOLD_SCRIPT"; exit 1; }

echo "[CHECK] FASTA_DIR=$FASTA_DIR"
test -d "$FASTA_DIR" || { echo "Missing: $FASTA_DIR"; exit 1; }

echo "[CHECK] MMCIF_DIR=$MMCIF_DIR"
test -d "$MMCIF_DIR" || { echo "Missing: $MMCIF_DIR"; exit 1; }

echo "[CHECK] ALIGN_DIR=$ALIGN_DIR"
test -d "$ALIGN_DIR" || { echo "Missing: $ALIGN_DIR"; exit 1; }

echo "[CHECK] CKPT=$CKPT"
test -f "$CKPT" || { echo "Missing: $CKPT"; exit 1; }

echo "Output dir: $OUT_DIR"

python3 "$OPENFOLD_SCRIPT" \
  "$FASTA_DIR" \
  "$MMCIF_DIR" \
  --output_dir "$OUT_DIR" \
  --use_precomputed_alignments "$ALIGN_DIR" \
  --config_preset model_1 \
  --openfold_checkpoint_path "$CKPT" \
  --model_device cuda:0

echo "===== Job finished on $(date) ====="
