#!/bin/bash -l
#SBATCH --job-name=of_eval
#SBATCH --partition=gpu-h100
#SBATCH --gres=gpu:h100:1
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=3-00:00:00
#SBATCH --array=0-1
#SBATCH --output=/users/%u/COMP702_BeyondAF/logs/slurm/of_eval_%A_%a.out
#SBATCH --error=/users/%u/COMP702_BeyondAF/logs/slurm/of_eval_%A_%a.err

# Unified COMP702 production pipeline for the non-2OLO evaluation targets.
# Array task 0 = 4BP8-A
# Array task 1 = 7DSQ-A
#
# Each task performs:
#   1) 100 official OpenFold predictions with unique deterministic seeds
#   2) preservation of unrelaxed PDBs
#   3) OpenFold Amber relaxation and preservation of relaxed PDBs
#   4) experimental model-1 chain-A reference extraction
#   5) BRI/LAI invariant calculation with pdb_assess
#   6) native pdb_assess plotting, preserving the appearance of the
#      plots already produced earlier in this project

# If someone accidentally sources this file, return instead of killing
# the interactive SSH shell.
if [[ "${BASH_SOURCE[0]}" != "$0" ]]; then
    echo "ERROR: do not source this file; submit it with sbatch"
    return 2
fi

# Never execute the heavy workload directly on a login node.
if [[ -z "${SLURM_JOB_ID:-}" ]]; then
    echo "ERROR: this is a Slurm batch script"
    echo "Submit with: sbatch $0"
    exit 2
fi

set -euo pipefail

source "$HOME/COMP702_BeyondAF/config/comp702_paths.sh"

TARGETS=(4bp8 7dsq)
CHAINS=(A A)
TARGET="${TARGETS[$SLURM_ARRAY_TASK_ID]}"
CHAIN="${CHAINS[$SLURM_ARRAY_TASK_ID]}"
TARGET_UPPER="${TARGET^^}"
REPEAT_COUNT="${OPENFOLD_REPEAT_COUNT:-100}"

RUN_ROOT="${COMP702_FASTSCRATCH_ROOT}/openfold_eval/${TARGET}_openfold_official_100"
RAW_ROOT="$RUN_ROOT/raw_openfold"
FASTA_RUN_DIR="$RUN_ROOT/fasta_input"
PDB_ROOT="$RUN_ROOT/pdbs"
UNRELAXED_ROOT="$PDB_ROOT/unrelaxed"
RELAXED_ROOT="$PDB_ROOT/relaxed"
ASSESS_INPUT="$RUN_ROOT/pdb_assess_input"
INVARIANT_ROOT="$RUN_ROOT/invariants"
PLOT_ROOT="$RUN_ROOT/plots"
PROV_ROOT="$RUN_ROOT/provenance"
BRI_ENV="${BRI_ENV:-$HOME/fastscratch/envs/bri_env_1.2.2}"

mkdir -p \
    "$RAW_ROOT" \
    "$FASTA_RUN_DIR" \
    "$UNRELAXED_ROOT" \
    "$RELAXED_ROOT" \
    "$ASSESS_INPUT" \
    "$INVARIANT_ROOT" \
    "$PLOT_ROOT" \
    "$PROV_ROOT" \
    "$TRITON_CACHE_DIR" \
    "$COMP702_LOG_DIR"

echo "============================================================"
echo "COMP702 OFFICIAL OPENFOLD -> PDB -> RELAX -> BRI/LAI"
echo "Target:       ${TARGET_UPPER}-${CHAIN}"
echo "Array task:   ${SLURM_ARRAY_TASK_ID}"
echo "Repeats:      ${REPEAT_COUNT}"
echo "Node:         $(hostname)"
echo "Started:      $(date --iso-8601=seconds)"
echo "Run root:     ${RUN_ROOT}"
echo "============================================================"

module load miniforge3/25.3.0-python3.12.10
eval "$(conda shell.bash hook)"
conda activate "$OPENFOLD_ENV"

export PYTHONPATH="$OPENFOLD_CODE_DIR:${PYTHONPATH:-}"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"

OPENFOLD_SCRIPT="$OPENFOLD_CODE_DIR/run_pretrained_openfold.py"
FASTA_FILE="$OPENFOLD_INFERENCE_FASTA_DIR/${TARGET}.fasta"
MMCIF_DIR="$OPENFOLD_INFERENCE_CIF_DIR"
MMCIF_FILE="$MMCIF_DIR/${TARGET}.cif"
ALIGN_ROOT="$OPENFOLD_INFERENCE_ALIGN_DIR"
ALIGN_TARGET="$ALIGN_ROOT/${TARGET}_${CHAIN}"
OF_CHECKPOINT="$OPENFOLD_OFFICIAL_FINETUNING_CKPT"

# Some source trees may use uppercase PDB filenames.
if [[ ! -f "$MMCIF_FILE" && -f "$MMCIF_DIR/${TARGET_UPPER}.cif" ]]; then
    MMCIF_FILE="$MMCIF_DIR/${TARGET_UPPER}.cif"
fi

for f in "$OPENFOLD_SCRIPT" "$FASTA_FILE" "$MMCIF_FILE" "$OF_CHECKPOINT"; do
    [[ -f "$f" ]] || { echo "FATAL: missing required file: $f"; exit 1; }
done
[[ -d "$ALIGN_TARGET" ]] || { echo "FATAL: missing alignment directory: $ALIGN_TARGET"; exit 1; }

# Isolate the target FASTA so OpenFold cannot accidentally infer other targets.
rm -f "$FASTA_RUN_DIR"/*.fasta
cp "$FASTA_FILE" "$FASTA_RUN_DIR/${TARGET}.fasta"
[[ "$(find "$FASTA_RUN_DIR" -maxdepth 1 -type f -name '*.fasta' | wc -l)" -eq 1 ]] || {
    echo "FATAL: target FASTA directory is not isolated"; exit 1;
}

echo "OpenFold script : $OPENFOLD_SCRIPT"
echo "FASTA           : $FASTA_FILE"
echo "mmCIF           : $MMCIF_FILE"
echo "Alignments      : $ALIGN_TARGET"
echo "OF checkpoint   : $OF_CHECKPOINT"
echo "Python          : $(which python)"
nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv,noheader

sha256sum "$FASTA_FILE" "$MMCIF_FILE" "$OF_CHECKPOINT" > "$PROV_ROOT/input_sha256.txt"
git -C "$COMP702_ROOT" rev-parse HEAD > "$PROV_ROOT/repository_git_head.txt"

# ------------------------------------------------------------------
# 1-3. OpenFold prediction, PDB output and relaxation
# ------------------------------------------------------------------
for REP in $(seq 1 "$REPEAT_COUNT"); do
    SEED=$((100000 * SLURM_ARRAY_TASK_ID + 1000 + REP))
    REP_DIR="$RAW_ROOT/repeat_${REP}"
    PRED_DIR="$REP_DIR/predictions"
    UNRELAXED="$PRED_DIR/${TARGET}_${CHAIN}_model_1_unrelaxed.pdb"
    RELAXED="$PRED_DIR/${TARGET}_${CHAIN}_model_1_relaxed.pdb"

    if [[ -s "$UNRELAXED" && -s "$RELAXED" ]]; then
        echo "[$(date --iso-8601=seconds)] ${TARGET} repeat ${REP}/${REPEAT_COUNT}: complete, skipping"
    else
        mkdir -p "$REP_DIR"
        echo "[$(date --iso-8601=seconds)] ${TARGET} repeat ${REP}/${REPEAT_COUNT}, seed=${SEED}"

        python "$OPENFOLD_SCRIPT" \
            "$FASTA_RUN_DIR" \
            "$MMCIF_DIR" \
            --output_dir "$REP_DIR" \
            --use_precomputed_alignments "$ALIGN_ROOT" \
            --config_preset model_1 \
            --openfold_checkpoint_path "$OF_CHECKPOINT" \
            --model_device cuda:0 \
            --data_random_seed "$SEED"
    fi

    [[ -s "$UNRELAXED" ]] || { echo "FATAL: missing unrelaxed PDB for repeat $REP"; exit 1; }
    [[ -s "$RELAXED" ]] || { echo "FATAL: missing relaxed PDB for repeat $REP"; exit 1; }

    cp -f "$UNRELAXED" "$UNRELAXED_ROOT/${TARGET}_${CHAIN}_model_1_unrelaxed_repeat${REP}.pdb"
    cp -f "$RELAXED" "$RELAXED_ROOT/${TARGET}_${CHAIN}_model_1_relaxed_repeat${REP}.pdb"
done

UNRELAXED_COUNT="$(find "$UNRELAXED_ROOT" -maxdepth 1 -type f -name '*.pdb' | wc -l)"
RELAXED_COUNT="$(find "$RELAXED_ROOT" -maxdepth 1 -type f -name '*.pdb' | wc -l)"
[[ "$UNRELAXED_COUNT" -eq "$REPEAT_COUNT" ]] || { echo "FATAL: expected $REPEAT_COUNT unrelaxed PDBs, found $UNRELAXED_COUNT"; exit 1; }
[[ "$RELAXED_COUNT" -eq "$REPEAT_COUNT" ]] || { echo "FATAL: expected $REPEAT_COUNT relaxed PDBs, found $RELAXED_COUNT"; exit 1; }

# ------------------------------------------------------------------
# Experimental reference: model 1, chain A
# ------------------------------------------------------------------
REFERENCE_PDB="$ASSESS_INPUT/${TARGET_UPPER}-1-${CHAIN}.pdb"
export MMCIF_FILE REFERENCE_PDB TARGET_UPPER CHAIN

python <<'PY'
import os
from Bio.PDB import MMCIFParser, PDBIO, Select

src = os.environ["MMCIF_FILE"]
dst = os.environ["REFERENCE_PDB"]
chain_id = os.environ["CHAIN"]

structure = MMCIFParser(QUIET=True).get_structure(os.environ["TARGET_UPPER"], src)
models = list(structure.get_models())
if not models:
    raise RuntimeError("mmCIF contains no models")
model = models[0]
if chain_id not in model:
    raise RuntimeError(f"model 1 contains no chain {chain_id}")

class FirstModelChain(Select):
    def accept_model(self, m):
        return int(m.id == model.id)
    def accept_chain(self, c):
        return int(c.id == chain_id)

io = PDBIO()
io.set_structure(structure)
io.save(dst, FirstModelChain())
if os.path.getsize(dst) == 0:
    raise RuntimeError("reference PDB is empty")
print(dst)
PY

# Copy only relaxed predictions into the invariant-analysis input; the
# unrelaxed predictions remain preserved separately for audit/use later.
find "$ASSESS_INPUT" -maxdepth 1 -type f -name "${TARGET}_${CHAIN}_model_1_relaxed_repeat*.pdb" -delete
cp "$RELAXED_ROOT"/*.pdb "$ASSESS_INPUT/"

# ------------------------------------------------------------------
# 4-5. BRI / LAI and plots
# ------------------------------------------------------------------
conda deactivate || true
conda activate "$BRI_ENV"

command -v pdb_assess >/dev/null || { echo "FATAL: pdb_assess not found in $BRI_ENV"; exit 1; }

rm -rf "$INVARIANT_ROOT" "$PLOT_ROOT"
mkdir -p "$INVARIANT_ROOT" "$PLOT_ROOT"

# Same invariant and plotting commands used for the earlier project plots.
# Reusing pdb_assess itself preserves their visual appearance/style instead
# of recreating the figures with a new matplotlib implementation.
pdb_assess inv "$ASSESS_INPUT" "$INVARIANT_ROOT" -n "${SLURM_CPUS_PER_TASK:-16}"

COMPUTE_DIR="$INVARIANT_ROOT/compute_invariants"
[[ -d "$COMPUTE_DIR" ]] || { echo "FATAL: compute_invariants directory missing"; exit 1; }

pdb_assess plot "$COMPUTE_DIR" "$PLOT_ROOT" -s "${TARGET_UPPER}-1-${CHAIN}"

find "$PLOT_ROOT" -type f | sort > "$PROV_ROOT/plot_files.txt"
find "$PLOT_ROOT" -type f -iname '*bri*.png' | grep -q . || { echo "FATAL: BRI plot not produced"; exit 1; }
find "$PLOT_ROOT" -type f -iname '*lai*.png' | grep -q . || { echo "FATAL: LAI plot not produced"; exit 1; }

sha256sum "$UNRELAXED_ROOT"/*.pdb > "$PROV_ROOT/unrelaxed_pdb_sha256.txt"
sha256sum "$RELAXED_ROOT"/*.pdb > "$PROV_ROOT/relaxed_pdb_sha256.txt"
find "$PLOT_ROOT" -type f -exec sha256sum {} \; > "$PROV_ROOT/plot_sha256.txt"

cat > "$RUN_ROOT/RUN_SUMMARY.txt" <<EOF
Target: ${TARGET_UPPER}-${CHAIN}
Model: official OpenFold finetuning_2.pt
Repeats: ${REPEAT_COUNT}
Unrelaxed PDBs: ${UNRELAXED_COUNT}
Relaxed PDBs: ${RELAXED_COUNT}
BRI/LAI input: relaxed predictions plus experimental ${TARGET_UPPER} model 1 chain ${CHAIN}
Invariant engine: pdb_assess
Plot engine: native pdb_assess plot
Plot style: intentionally identical mechanism to the project's earlier BRI/LAI plots
Run root: ${RUN_ROOT}
Completed: $(date --iso-8601=seconds)
EOF

echo "============================================================"
echo "PIPELINE PASS: ${TARGET_UPPER}-${CHAIN}"
echo "Unrelaxed PDBs: $UNRELAXED_COUNT"
echo "Relaxed PDBs:   $RELAXED_COUNT"
echo "BRI plot(s):"
find "$PLOT_ROOT" -type f -iname '*bri*.png' -print
echo "LAI plot(s):"
find "$PLOT_ROOT" -type f -iname '*lai*.png' -print
echo "Run root: $RUN_ROOT"
echo "Finished: $(date --iso-8601=seconds)"
echo "============================================================"
