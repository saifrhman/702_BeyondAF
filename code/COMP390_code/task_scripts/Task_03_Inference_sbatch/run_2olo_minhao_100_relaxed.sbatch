#!/bin/bash
#SBATCH --job-name=2olo_100_relaxed
#SBATCH --partition=gpu-h100
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=2-00:00:00
#SBATCH --output=/users/%u/COMP702_BeyondAF/logs/slurm/2olo_100_relaxed_%j.log
#SBATCH --error=/users/%u/COMP702_BeyondAF/logs/slurm/2olo_100_relaxed_%j.err

set -euo pipefail

source "$HOME/COMP702_BeyondAF/config/comp702_paths.sh"

module load miniforge3/25.3.0-python3.12.10
eval "$(conda shell.bash hook)"
conda activate "$OPENFOLD_ENV"

export PYTHONPATH="$OPENFOLD_CODE_DIR:${PYTHONPATH:-}"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"

OPENFOLD_SCRIPT="$OPENFOLD_CODE_DIR/run_pretrained_openfold.py"
FASTA_DIR="$OPENFOLD_INFERENCE_FASTA_DIR"
CIF_DIR="$OPENFOLD_INFERENCE_CIF_DIR"
ALIGN_DIR="$OPENFOLD_INFERENCE_ALIGN_DIR"
CHECKPOINT="$OPENFOLD_EP32_CKPT"

REPEAT_COUNT="${OPENFOLD_2OLO_REPEAT_COUNT:-100}"
RUN_ROOT="${OPENFOLD_2OLO_RELAXED_ROOT:?OPENFOLD_2OLO_RELAXED_ROOT is required}"

mkdir -p "$RUN_ROOT" "$TRITON_CACHE_DIR"

test -f "$OPENFOLD_SCRIPT"
test -f "$CHECKPOINT"
test -d "$FASTA_DIR"
test -d "$CIF_DIR"
test -d "$ALIGN_DIR"

echo "Run root: $RUN_ROOT"
echo "Repeats: $REPEAT_COUNT"
echo "Checkpoint: $CHECKPOINT"

for REP in $(seq 1 "$REPEAT_COUNT"); do
    SEED=$((1000 + REP))
    OUT_DIR="$RUN_ROOT/repeat_${REP}"

    mkdir -p "$OUT_DIR"

    echo "[$(date)] Repeat $REP/$REPEAT_COUNT, seed $SEED"

    python3 "$OPENFOLD_SCRIPT" \
        "$FASTA_DIR" \
        "$CIF_DIR" \
        --output_dir "$OUT_DIR" \
        --use_precomputed_alignments "$ALIGN_DIR" \
        --config_preset model_1 \
        --openfold_checkpoint_path "$CHECKPOINT" \
        --model_device cuda:0 \
        --data_random_seed "$SEED"
done

echo "[$(date)] Finished all $REPEAT_COUNT repeats"

echo "PDB count:"
find "$RUN_ROOT" -type f -name '*.pdb' | wc -l
