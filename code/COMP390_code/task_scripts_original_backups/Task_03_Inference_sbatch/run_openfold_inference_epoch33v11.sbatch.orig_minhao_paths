#!/bin/bash -l
#SBATCH -J of_infer_ep15_20_x5
#SBATCH -p gpu-h100
#SBATCH -N 1
#SBATCH -w gpu32
#SBATCH --gres=gpu:h100:1
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=2-00:00:00
#SBATCH --output=/users/%u/fastscratch/of_infer_ep15_20_x5_%j.log
#SBATCH --error=/users/%u/fastscratch/of_infer_ep15_20_x5_%j.err

set -euo pipefail

########################################
# Activate environment
########################################
set +u
source /mnt/data2/users/sgmwu14/tools/miniforge3/etc/profile.d/conda.sh
conda activate openfold_env
set -u

########################################
# Environment variables
########################################
export PYTHONPATH=/users/sgmwu14/openfold:$PYTHONPATH
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export TRITON_CACHE_DIR=/users/sgmwu14/fastscratch/triton_cache
mkdir -p "$TRITON_CACHE_DIR"

########################################
# Shared input paths
########################################
FASTA_DIR=/users/sgmwu14/scratch/openfold_inference/fasta_dir
CIF_DIR=/users/sgmwu14/scratch/taskcif
ALIGN_DIR=/users/sgmwu14/scratch/openfold_inference/alignments_cf
OUT_BASE=/users/sgmwu14/scratch/openfold_inference_runs

mkdir -p "$OUT_BASE"

########################################
# Checkpoint list
########################################
CKPTS=(
  /users/sgmwu14/scratch/openfold_runs/run_subset_gpu32_h100_ep21/csv_logs/version_0/checkpoints/15-160000.ckpt
  /users/sgmwu14/scratch/openfold_runs/run_subset_gpu32_h100_ep21/csv_logs/version_0/checkpoints/16-170000.ckpt
  /users/sgmwu14/scratch/openfold_runs/run_subset_gpu32_h100_ep21/csv_logs/version_0/checkpoints/17-180000.ckpt
  /users/sgmwu14/scratch/openfold_runs/run_subset_gpu32_h100_ep21/csv_logs/version_0/checkpoints/18-190000.ckpt
  /users/sgmwu14/scratch/openfold_runs/run_subset_gpu32_h100_ep21/csv_logs/version_0/checkpoints/19-200000.ckpt
  /users/sgmwu14/scratch/openfold_runs/run_subset_gpu32_h100_ep21/csv_logs/version_0/checkpoints/20-210000.ckpt
  /users/sgmwu14/scratch/openfold_runs/run_subset_gpu32_h100_ep27/csv_logs/version_0/checkpoints/21-220000.ckpt
)

########################################
# Main loop
########################################
for MY_CKPT in "${CKPTS[@]}"; do
    if [[ ! -f "$MY_CKPT" ]]; then
        echo "[$(date)] ERROR: checkpoint not found: $MY_CKPT"
        exit 1
    fi

    CKPT_FILE=$(basename "$MY_CKPT")   # e.g. 15-160000.ckpt
    CKPT_TAG=${CKPT_FILE%%-*}          # e.g. 15

    RUN_ROOT=${OUT_BASE}/epoch_${CKPT_TAG}
    MY_ROOT=${RUN_ROOT}/your_model_epoch_${CKPT_TAG}
    mkdir -p "$MY_ROOT"

    echo "[$(date)] ========================================"
    echo "[$(date)] Running checkpoint: $MY_CKPT"
    echo "[$(date)] Epoch tag  : $CKPT_TAG"
    echo "[$(date)] Output dir : $MY_ROOT"
    echo "[$(date)] ========================================"

    ########################################
    # Repeat inference 5 times
    ########################################
    for REP in 1 2 3 4 5; do
        SEED=$((100 + REP))
        MY_OUT=${MY_ROOT}/repeat_${REP}
        mkdir -p "$MY_OUT"

        echo "[$(date)] ---- Checkpoint ${CKPT_TAG}, repeat ${REP}, seed ${SEED} ----"

        python3 /users/sgmwu14/openfold/run_pretrained_openfold.py \
          "$FASTA_DIR" \
          "$CIF_DIR" \
          --output_dir "$MY_OUT" \
          --use_precomputed_alignments "$ALIGN_DIR" \
          --config_preset model_1 \
          --openfold_checkpoint_path "$MY_CKPT" \
          --model_device cuda:0 \
          --data_random_seed "$SEED" 
    done

    echo "[$(date)] ===== FINISHED checkpoint epoch_${CKPT_TAG} ====="
done

echo "[$(date)] ===== ALL CHECKPOINTS FINISHED ====="