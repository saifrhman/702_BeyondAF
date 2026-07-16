#!/bin/bash -l
#SBATCH -J of_infer_ep32_x5_3models
#SBATCH -p gpu-h100
#SBATCH -N 1
#SBATCH --gres=gpu:h100:1
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=1-00:00:00
#SBATCH --output=/users/%u/COMP702_BeyondAF/logs/slurm/of_infer_ep32_x5_3models_%j.log
#SBATCH --error=/users/%u/COMP702_BeyondAF/logs/slurm/of_infer_ep32_x5_3models_%j.err

set -euo pipefail

source "$HOME/COMP702_BeyondAF/config/comp702_paths.sh"

module load miniforge3/25.3.0-python3.12.10
eval "$(conda shell.bash hook)"
module load miniforge3/25.3.0-python3.12.10
source "$(conda info --base)/etc/profile.d/conda.sh"
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
MY_CKPT="${OPENFOLD_3MODEL_CKPT:-$OPENFOLD_EP32_CKPT}"

RUN_ROOT="${OPENFOLD_3MODEL_X5_ROOT:-$OPENFOLD_INFERENCE_RUNS_ROOT/epoch_32_3models_x5}"

AF_ROOT="$RUN_ROOT/alphafold_params"
OF_ROOT="$RUN_ROOT/openfold_official"
MY_ROOT="$RUN_ROOT/trained_model_epoch_32"

REPEAT_COUNT="${OPENFOLD_SWEEP_REPEAT_COUNT:-5}"

mkdir -p "$RUN_ROOT" "$AF_ROOT" "$OF_ROOT" "$MY_ROOT"

test -f "$OPENFOLD_SCRIPT" || { echo "Missing: $OPENFOLD_SCRIPT"; exit 1; }
test -d "$FASTA_DIR" || { echo "Missing: $FASTA_DIR"; exit 1; }
test -d "$CIF_DIR" || { echo "Missing: $CIF_DIR"; exit 1; }
test -d "$ALIGN_DIR" || { echo "Missing: $ALIGN_DIR"; exit 1; }
test -f "$AF_PARAM" || { echo "Missing: $AF_PARAM"; exit 1; }
test -f "$OF_CKPT" || { echo "Missing: $OF_CKPT"; exit 1; }
test -f "$MY_CKPT" || { echo "Missing: $MY_CKPT"; exit 1; }

for REP in $(seq 1 "$REPEAT_COUNT"); do
  SEED=$((1000 + REP))

  AF_OUT="$AF_ROOT/repeat_${REP}"
  OF_OUT="$OF_ROOT/repeat_${REP}"
  MY_OUT="$MY_ROOT/repeat_${REP}"

  mkdir -p "$AF_OUT" "$OF_OUT" "$MY_OUT"

  echo "[$(date)] Repeat $REP/$REPEAT_COUNT, seed $SEED"

  python3 "$OPENFOLD_SCRIPT" \
    "$FASTA_DIR" \
    "$CIF_DIR" \
    --output_dir "$AF_OUT" \
    --use_precomputed_alignments "$ALIGN_DIR" \
    --config_preset model_1 \
    --jax_param_path "$AF_PARAM" \
    --model_device cuda:0 \
    --data_random_seed "$SEED" \
    --skip_relaxation

  python3 "$OPENFOLD_SCRIPT" \
    "$FASTA_DIR" \
    "$CIF_DIR" \
    --output_dir "$OF_OUT" \
    --use_precomputed_alignments "$ALIGN_DIR" \
    --config_preset model_1 \
    --openfold_checkpoint_path "$OF_CKPT" \
    --model_device cuda:0 \
    --data_random_seed "$SEED" \
    --skip_relaxation

  python3 "$OPENFOLD_SCRIPT" \
    "$FASTA_DIR" \
    "$CIF_DIR" \
    --output_dir "$MY_OUT" \
    --use_precomputed_alignments "$ALIGN_DIR" \
    --config_preset model_1 \
    --openfold_checkpoint_path "$MY_CKPT" \
    --model_device cuda:0 \
    --data_random_seed "$SEED" \
    --skip_relaxation
done

echo "[$(date)] Three-model x5 inference finished."
