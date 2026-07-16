#!/bin/bash
#SBATCH -J of_4bp8_100pred_ep32
#SBATCH -p gpu-h100
#SBATCH -N 1
#SBATCH -w gpu32
#SBATCH --gres=gpu:h100:1
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=3-00:00:00
#SBATCH -o /users/sgmwu14/scratch/logs/of_4bp8_100pred_ep32_%j.out
#SBATCH -e /users/sgmwu14/scratch/logs/of_4bp8_100pred_ep32_%j.err

set -euo pipefail

echo "===== Job started on $(date) ====="
echo "Hostname: $(hostname)"
echo "CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-unset}"

# -----------------------------
# Paths
# -----------------------------
OPENFOLD_SCRIPT="/users/sgmwu14/openfold/run_pretrained_openfold.py"
FASTA_FILE="/users/sgmwu14/scratch/openfold_inference/fasta_dir/4bp8.fasta"
MMCIF_DIR="/users/sgmwu14/localscratch/pdb_data_full/mmcif_files"
ALIGN_DIR="/users/sgmwu14/localscratch/alignment_openfold_full_20260210/alignment_data/alignments"
CKPT="/users/sgmwu14/scratch/openfold_runs/run_subset_gpu32_h100_ep33/csv_logs/version_0/checkpoints/32-330000.ckpt"

# -----------------------------
# Output directories
# -----------------------------
BASE_OUT="/users/sgmwu14/scratch/openfold_inference_4bp8_ep32_100pred"
FASTA_DIR="${BASE_OUT}/fasta_input"
LOG_DIR="/users/sgmwu14/scratch/logs"

mkdir -p "${BASE_OUT}"
mkdir -p "${FASTA_DIR}"
mkdir -p "${LOG_DIR}"

# only keep 4bp8 fasta in a dedicated folder
cp "${FASTA_FILE}" "${FASTA_DIR}/4bp8.fasta"

echo "Output base: ${BASE_OUT}"
echo "FASTA dir:    ${FASTA_DIR}"
echo "Checkpoint:   ${CKPT}"

# -----------------------------
# Optional: activate env
# Uncomment and edit if needed
# -----------------------------
# source /users/sgmwu14/.bashrc
# conda activate openfold

# -----------------------------
# Run 100 predictions
# Each repeat gets its own output folder
# -----------------------------
for i in $(seq 1 100); do
    REPEAT_DIR="${BASE_OUT}/repeat_${i}"
    mkdir -p "${REPEAT_DIR}"

    echo "===== Running repeat ${i}/100 ====="
    echo "Output -> ${REPEAT_DIR}"

    python3 "${OPENFOLD_SCRIPT}" \
      "${FASTA_DIR}" \
      "${MMCIF_DIR}" \
      --output_dir "${REPEAT_DIR}" \
      --use_precomputed_alignments "${ALIGN_DIR}" \
      --config_preset model_1 \
      --openfold_checkpoint_path "${CKPT}" \
      --model_device cuda:0

    echo "===== Finished repeat ${i}/100 at $(date) ====="
done

echo "===== All 100 predictions completed at $(date) ====="