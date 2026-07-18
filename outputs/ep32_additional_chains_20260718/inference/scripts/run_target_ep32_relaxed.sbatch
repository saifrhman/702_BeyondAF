#!/bin/bash
#SBATCH --job-name=target_ep32
#SBATCH --partition=gpu-h100
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=3-00:00:00

set -euo pipefail

source "$HOME/COMP702_BeyondAF/config/comp702_paths.sh"

PDB_ID="${PDB_ID:?PDB_ID is required}"
REPEAT_COUNT="${REPEAT_COUNT:-1}"
RUN_LABEL="${RUN_LABEL:-smoke}"
SEED_BASE="${SEED_BASE:-1000}"

TEMPLATE_MMCIF_DIR="${TEMPLATE_MMCIF_DIR:?TEMPLATE_MMCIF_DIR is required}"
EP32_RUNS_DIR="${EP32_RUNS_DIR:?EP32_RUNS_DIR is required}"

PDB_ID="$(
    printf '%s' "$PDB_ID" |
    tr '[:upper:]' '[:lower:]'
)"

CHAIN_ID="${PDB_ID}_A"

module load miniforge3/25.3.0-python3.12.10
eval "$(conda shell.bash hook)"
conda activate "$OPENFOLD_ENV"

export PYTHONPATH="$OPENFOLD_CODE_DIR:${PYTHONPATH:-}"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"

OPENFOLD_SCRIPT="$OPENFOLD_CODE_DIR/run_pretrained_openfold.py"
CHECKPOINT="$OPENFOLD_EP32_CKPT"

SOURCE_FASTA="$OPENFOLD_INFERENCE_FASTA_DIR/${PDB_ID}.fasta"
ALIGNMENT_ROOT="$OPENFOLD_INFERENCE_ALIGN_DIR"
TARGET_ALIGNMENT="$ALIGNMENT_ROOT/${CHAIN_ID}"

JOB_IDENTIFIER="${SLURM_JOB_ID:-manual}"
RUN_ROOT="$EP32_RUNS_DIR/${PDB_ID}_ep32_${RUN_LABEL}_${JOB_IDENTIFIER}"
FASTA_DIR="$RUN_ROOT/fasta_input"

mkdir -p \
    "$RUN_ROOT" \
    "$FASTA_DIR" \
    "$TRITON_CACHE_DIR"

echo "===== Input checks ====="

echo "OpenFold script: $OPENFOLD_SCRIPT"
test -f "$OPENFOLD_SCRIPT" || {
    echo "ERROR: OpenFold script missing"
    exit 1
}

echo "Checkpoint: $CHECKPOINT"
test -f "$CHECKPOINT" || {
    echo "ERROR: epoch-32 checkpoint missing"
    exit 1
}

echo "FASTA: $SOURCE_FASTA"
test -s "$SOURCE_FASTA" || {
    echo "ERROR: FASTA missing or empty"
    exit 1
}

echo "Alignment directory: $TARGET_ALIGNMENT"
test -d "$TARGET_ALIGNMENT" || {
    echo "ERROR: alignment directory missing"
    exit 1
}

ALIGNMENT_COUNT=$(
    find -L "$TARGET_ALIGNMENT" -type f |
    wc -l
)

if [ "$ALIGNMENT_COUNT" -eq 0 ]; then
    echo "ERROR: no alignment files found"
    exit 1
fi

echo "Template directory: $TEMPLATE_MMCIF_DIR"
test -d "$TEMPLATE_MMCIF_DIR" || {
    echo "ERROR: template mmCIF directory missing"
    exit 1
}

TEMPLATE_COUNT=$(
    find "$TEMPLATE_MMCIF_DIR" \
        -maxdepth 1 \
        -type f \
        -name '*.cif' |
    wc -l
)

if [ "$TEMPLATE_COUNT" -eq 0 ]; then
    echo "ERROR: no template mmCIF files found"
    exit 1
fi

# Use a target-only FASTA directory so OpenFold does not run other proteins.
cp "$SOURCE_FASTA" "$FASTA_DIR/${PDB_ID}.fasta"

cat > "$RUN_ROOT/run_manifest.txt" <<MANIFEST
Target: $CHAIN_ID
Model: retrained OpenFold epoch 32
Checkpoint: $CHECKPOINT
Run label: $RUN_LABEL
Run root: $RUN_ROOT
Source FASTA: $SOURCE_FASTA
Target alignment directory: $TARGET_ALIGNMENT
Alignment file count: $ALIGNMENT_COUNT
Template mmCIF directory: $TEMPLATE_MMCIF_DIR
Template mmCIF count: $TEMPLATE_COUNT
Repeat count: $REPEAT_COUNT
Seed base: $SEED_BASE
Seed formula: seed base + repeat number
Slurm job ID: $JOB_IDENTIFIER
Started: $(date)
MANIFEST

echo "=============================================="
echo "Target:               $CHAIN_ID"
echo "Run root:             $RUN_ROOT"
echo "Alignment files:      $ALIGNMENT_COUNT"
echo "Template structures:  $TEMPLATE_COUNT"
echo "Repeats:              $REPEAT_COUNT"
echo "Seeds:                $((SEED_BASE + 1)) to $((SEED_BASE + REPEAT_COUNT))"
echo "=============================================="

for REP in $(seq 1 "$REPEAT_COUNT"); do
    SEED=$((SEED_BASE + REP))
    OUT_DIR="$RUN_ROOT/repeat_${REP}"

    mkdir -p "$OUT_DIR"

    echo "[$(date)] Repeat $REP/$REPEAT_COUNT, seed $SEED"

    python3 "$OPENFOLD_SCRIPT" \
        "$FASTA_DIR" \
        "$TEMPLATE_MMCIF_DIR" \
        --output_dir "$OUT_DIR" \
        --use_precomputed_alignments "$ALIGNMENT_ROOT" \
        --config_preset model_1 \
        --openfold_checkpoint_path "$CHECKPOINT" \
        --model_device cuda:0 \
        --data_random_seed "$SEED"
done

UNRELAXED_COUNT=$(
    find "$RUN_ROOT" \
        -type f \
        -name "${CHAIN_ID}_model_1_unrelaxed.pdb" |
    wc -l
)

RELAXED_COUNT=$(
    find "$RUN_ROOT" \
        -type f \
        -name "${CHAIN_ID}_model_1_relaxed.pdb" |
    wc -l
)

TOTAL_COUNT=$((UNRELAXED_COUNT + RELAXED_COUNT))

cat >> "$RUN_ROOT/run_manifest.txt" <<MANIFEST
Finished: $(date)
Unrelaxed PDB count: $UNRELAXED_COUNT
Relaxed PDB count: $RELAXED_COUNT
Total PDB count: $TOTAL_COUNT
MANIFEST

echo "Unrelaxed PDBs: $UNRELAXED_COUNT"
echo "Relaxed PDBs:   $RELAXED_COUNT"
echo "Total PDBs:     $TOTAL_COUNT"

if [ "$UNRELAXED_COUNT" -ne "$REPEAT_COUNT" ]; then
    echo "ERROR: incorrect number of unrelaxed PDB files"
    exit 1
fi

if [ "$RELAXED_COUNT" -ne "$REPEAT_COUNT" ]; then
    echo "ERROR: incorrect number of relaxed PDB files"
    exit 1
fi

echo "[$(date)] Completed successfully"
