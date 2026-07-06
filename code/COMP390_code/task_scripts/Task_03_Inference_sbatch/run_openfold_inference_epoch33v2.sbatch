#!/bin/bash -l
#SBATCH -J of_infer3_ep33_v2
#SBATCH -p gpu-h100
#SBATCH -N 1
#SBATCH --gres=gpu:h100:1
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=1-00:00:00
#SBATCH --output=/users/%u/COMP702_BeyondAF/logs/slurm/of_infer3_ep33_v2_%j.log
#SBATCH --error=/users/%u/COMP702_BeyondAF/logs/slurm/of_infer3_ep33_v2_%j.err

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

RUN_ROOT="${OPENFOLD_EPOCH33_RERUN4_ROOT:-$OPENFOLD_INFERENCE_OUT_ROOT/epoch33_rerun4}"

AF_OUT="$RUN_ROOT/alphafold_params_output"
OF_OUT="$RUN_ROOT/openfold_official_output"
MY_OUT="$RUN_ROOT/trained_model_ep32_output"

mkdir -p "$RUN_ROOT" "$AF_OUT" "$OF_OUT" "$MY_OUT"

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

echo "[$(date)] Running AlphaFold parameter model, unrelaxed only..."
python3 "$OPENFOLD_SCRIPT" \
  "$FASTA_DIR" \
  "$CIF_DIR" \
  --output_dir "$AF_OUT" \
  --use_precomputed_alignments "$ALIGN_DIR" \
  --config_preset model_1 \
  --jax_param_path "$AF_PARAM" \
  --model_device cuda:0 \
  --skip_relaxation

echo "[$(date)] Running official OpenFold model, unrelaxed only..."
python3 "$OPENFOLD_SCRIPT" \
  "$FASTA_DIR" \
  "$CIF_DIR" \
  --output_dir "$OF_OUT" \
  --use_precomputed_alignments "$ALIGN_DIR" \
  --config_preset model_1 \
  --openfold_checkpoint_path "$OF_CKPT" \
  --model_device cuda:0 \
  --skip_relaxation

echo "[$(date)] Running archived trained checkpoint, unrelaxed only..."
python3 "$OPENFOLD_SCRIPT" \
  "$FASTA_DIR" \
  "$CIF_DIR" \
  --output_dir "$MY_OUT" \
  --use_precomputed_alignments "$ALIGN_DIR" \
  --config_preset model_1 \
  --openfold_checkpoint_path "$MY_CKPT" \
  --model_device cuda:0 \
  --skip_relaxation

echo "[$(date)] All unrelaxed inference finished."

echo "[$(date)] Organising outputs and running 5 relaxation repeats..."

for MODEL_DIR in "$AF_OUT" "$OF_OUT" "$MY_OUT"; do
    MODEL_NAME="$(basename "$MODEL_DIR")"
    PRED_DIR="${MODEL_DIR}/predictions"

    echo "[$(date)] Processing model: $MODEL_NAME"

    if [ ! -d "$PRED_DIR" ]; then
        echo "[$(date)] No predictions directory found for $MODEL_NAME, skipping relaxation."
        continue
    fi

    find "$PRED_DIR" -maxdepth 1 -type f -name "*unrelaxed*.pdb" | while read -r pdb; do
        fname="$(basename "$pdb")"
        base="$(basename "$pdb" .pdb)"

        TARGET_DIR="${PRED_DIR}/${base}"
        mkdir -p "$TARGET_DIR"

        if [ ! -f "${TARGET_DIR}/${fname}" ]; then
            mv "$pdb" "$TARGET_DIR/"
        fi

        pdb="${TARGET_DIR}/${fname}"

        echo "[$(date)] Processing $MODEL_NAME: $base"
        echo "[$(date)] Unrelaxed file: $pdb"

        for i in 1 2 3 4 5; do
            out_pdb="${TARGET_DIR}/${base}_relax${i}.pdb"
            out_dbg="${TARGET_DIR}/${base}_relax${i}_debug.txt"

            echo "[$(date)]   Relax run $i -> $out_pdb"

            PDB_IN="$pdb" PDB_OUT="$out_pdb" DBG_OUT="$out_dbg" python3 - <<'PYRELAX' || echo "[$(date)] WARNING: relax failed for ${MODEL_NAME} ${base} run ${i}"
import os
from openfold.np import protein
from openfold.np.relax.relax import AmberRelaxation

pdb_in = os.environ["PDB_IN"]
pdb_out = os.environ["PDB_OUT"]
dbg_out = os.environ["DBG_OUT"]

with open(pdb_in, "r") as f:
    pdb_str = f.read()

prot = protein.from_pdb_string(pdb_str)

relaxer = AmberRelaxation(
    max_iterations=0,
    tolerance=2.39,
    stiffness=10.0,
    exclude_residues=[],
    max_outer_iterations=20,
    use_gpu=True,
)

relaxed_pdb, debug_data, violations = relaxer.process(prot=prot)

with open(pdb_out, "w") as f:
    f.write(relaxed_pdb)

with open(dbg_out, "w") as f:
    f.write(str(debug_data) + "\n")
    f.write("violations_sum=" + str(violations.sum()) + "\n")
PYRELAX
        done
    done
done

echo "[$(date)] All 5x relaxation finished."
