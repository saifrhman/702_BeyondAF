#!/bin/bash -l
#SBATCH -J of_infer_ep21_26_x5
#SBATCH -p gpu-h100
#SBATCH -N 1
#SBATCH --gres=gpu:h100:1
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=3-00:00:00
#SBATCH --output=/users/%u/COMP702_BeyondAF/logs/slurm/of_infer_ep21_26_x5_%j.log
#SBATCH --error=/users/%u/COMP702_BeyondAF/logs/slurm/of_infer_ep21_26_x5_%j.err

set -euo pipefail

source "$HOME/COMP702_BeyondAF/config/comp702_paths.sh"

echo "======================================"
echo "Job ID: ${SLURM_JOB_ID:-manual}"
echo "Node: $(hostname)"
echo "Start time: $(date)"
echo "CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-unset}"
echo "======================================"

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

OUT_BASE="${OPENFOLD_SWEEP_OUT_BASE:-$OPENFOLD_INFERENCE_RUNS_ROOT/sweep_ep21_26}"
REPEAT_COUNT="${OPENFOLD_SWEEP_REPEAT_COUNT:-5}"
DEFAULT_EPOCHS="21 22 23 24 25 26"

mkdir -p "$OUT_BASE"

echo "[CHECK] OPENFOLD_SCRIPT=$OPENFOLD_SCRIPT"
test -f "$OPENFOLD_SCRIPT" || { echo "Missing: $OPENFOLD_SCRIPT"; exit 1; }

echo "[CHECK] FASTA_DIR=$FASTA_DIR"
test -d "$FASTA_DIR" || { echo "Missing: $FASTA_DIR"; exit 1; }

echo "[CHECK] CIF_DIR=$CIF_DIR"
test -d "$CIF_DIR" || { echo "Missing: $CIF_DIR"; exit 1; }

echo "[CHECK] ALIGN_DIR=$ALIGN_DIR"
test -d "$ALIGN_DIR" || { echo "Missing: $ALIGN_DIR"; exit 1; }

checkpoint_for_epoch () {
  local epoch="$1"
  local step=$(( (epoch + 1) * 10000 ))
  local block=""

  if (( epoch <= 7 )); then
    block="ep10"
  elif (( epoch <= 14 )); then
    block="ep15"
  elif (( epoch <= 20 )); then
    block="ep21"
  elif (( epoch <= 26 )); then
    block="ep27"
  else
    block="ep33"
  fi

  if [[ "$epoch" == "32" && -f "${OPENFOLD_EP32_CKPT:-}" ]]; then
    echo "$OPENFOLD_EP32_CKPT"
    return
  fi

  local primary="$OPENFOLD_TRAINED_CKPT_ROOT/run_subset_h100_${block}/csv_logs/version_0/checkpoints/${epoch}-${step}.ckpt"
  local archived="$OPENFOLD_ARCHIVE_CKPT_ROOT/run_subset_gpu32_h100_${block}/csv_logs/version_0/checkpoints/${epoch}-${step}.ckpt"

  if [[ -f "$primary" ]]; then
    echo "$primary"
  elif [[ -f "$archived" ]]; then
    echo "$archived"
  else
    echo "$primary"
  fi
}

for EPOCH in $DEFAULT_EPOCHS; do
  MY_CKPT="$(checkpoint_for_epoch "$EPOCH")"

  if [[ ! -f "$MY_CKPT" ]]; then
    echo "[$(date)] WARNING: checkpoint not found, skipping epoch $EPOCH: $MY_CKPT"
    continue
  fi

  CKPT_FILE="$(basename "$MY_CKPT")"
  CKPT_TAG="${CKPT_FILE%%-*}"

  RUN_ROOT="$OUT_BASE/epoch_${CKPT_TAG}"
  MY_ROOT="$RUN_ROOT/trained_model_epoch_${CKPT_TAG}"
  mkdir -p "$MY_ROOT"

  echo "[$(date)] Running checkpoint: $MY_CKPT"
  echo "[$(date)] Epoch tag: $CKPT_TAG"
  echo "[$(date)] Output root: $RUN_ROOT"
  echo "[$(date)] Repeats: $REPEAT_COUNT"

  for REP in $(seq 1 "$REPEAT_COUNT"); do
    SEED=$((1000 + REP))
    REP_OUT="$MY_ROOT/repeat_${REP}"
    mkdir -p "$REP_OUT"

    echo "[$(date)] Checkpoint ${CKPT_TAG}, repeat ${REP}, seed ${SEED}"

    python3 "$OPENFOLD_SCRIPT" \
      "$FASTA_DIR" \
      "$CIF_DIR" \
      --output_dir "$REP_OUT" \
      --use_precomputed_alignments "$ALIGN_DIR" \
      --config_preset model_1 \
      --openfold_checkpoint_path "$MY_CKPT" \
      --model_device cuda:0 \
      --data_random_seed "$SEED" \
      --skip_relaxation
  done

  echo "[$(date)] Finished checkpoint epoch_${CKPT_TAG}"
done

echo "[$(date)] All selected checkpoints finished."
