#!/bin/bash -l
#SBATCH -J of_7dsq_100pred_af
#SBATCH -p gpu-h100
#SBATCH -N 1
#SBATCH --gres=gpu:h100:1
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=3-00:00:00
#SBATCH -o /users/%u/scratch/logs/of_7dsq_100pred_af_%j.out
#SBATCH -e /users/%u/scratch/logs/of_7dsq_100pred_af_%j.err

set -euo pipefail

echo "===== Job started on $(date) ====="
echo "Hostname: $(hostname)"
echo "CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-unset}"

# -----------------------------
# Activate environment
# -----------------------------
set +u
source /mnt/data2/users/sgmwu14/tools/miniforge3/etc/profile.d/conda.sh
conda activate openfold_env
set -u

# -----------------------------
# Environment variables
# -----------------------------
export PYTHONPATH=/users/sgmwu14/openfold:$PYTHONPATH
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export TRITON_CACHE_DIR=/users/sgmwu14/fastscratch/triton_cache
mkdir -p "$TRITON_CACHE_DIR"

# -----------------------------
# Input paths
# -----------------------------
OPENFOLD_SCRIPT="/users/sgmwu14/openfold/run_pretrained_openfold.py"
FASTA_FILE="/users/sgmwu14/scratch/openfold_inference/fasta_dir/7dsq.fasta"
MMCIF_DIR="/users/sgmwu14/scratch/taskcif"
ALIGN_ROOT="/users/sgmwu14/scratch/openfold_inference/alignments_cf"
JAX_PARAM="/users/sgmwu14/openfold/openfold/resources/params/params_model_1.npz"

# -----------------------------
# Output paths
# -----------------------------
BASE_OUT="/users/sgmwu14/scratch/openfold_inference_7dsq_100pred_af"
FASTA_DIR="${BASE_OUT}/fasta_input"
LOG_DIR="/users/sgmwu14/scratch/logs"

mkdir -p "${BASE_OUT}" "${FASTA_DIR}" "${LOG_DIR}"

# Keep only 7dsq fasta in dedicated folder
cp "${FASTA_FILE}" "${FASTA_DIR}/7dsq.fasta"

echo "FASTA file : ${FASTA_FILE}"
echo "MMCIF dir  : ${MMCIF_DIR}"
echo "Align root : ${ALIGN_ROOT}"
echo "JAX params : ${JAX_PARAM}"
echo "Output dir : ${BASE_OUT}"

# -----------------------------
# Run 100 predictions
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
      --use_precomputed_alignments "${ALIGN_ROOT}" \
      --config_preset model_1 \
      --jax_param_path "${JAX_PARAM}" \
      --model_device cuda:0

    echo "===== Finished repeat ${i}/100 at $(date) ====="
done

echo "===== All 100 predictions completed at $(date) ====="