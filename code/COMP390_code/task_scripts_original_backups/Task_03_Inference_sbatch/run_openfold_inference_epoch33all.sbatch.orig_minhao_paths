#!/bin/bash -l
#SBATCH -J of_infer_ep0_32_x5
#SBATCH -p gpu-h100
#SBATCH -N 1
#SBATCH -w gpu32
#SBATCH --gres=gpu:h100:1
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=3-00:00:00
#SBATCH --output=/users/%u/fastscratch/of_infer_ep0_32_x5_%j.log
#SBATCH --error=/users/%u/fastscratch/of_infer_ep0_32_x5_%j.err

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
OUT_BASE=/users/sgmwu14/scratch/openfold_inference_runs2

mkdir -p "$OUT_BASE"

########################################
# Checkpoint list
########################################
CKPTS=(
  /users/sgmwu14/scratch/openfold_runs/run_subset_gpu32_h100_ep10/csv_logs/version_0/checkpoints/0-10000.ckpt
  /users/sgmwu14/scratch/openfold_runs/run_subset_gpu32_h100_ep10/csv_logs/version_0/checkpoints/1-20000.ckpt
  /users/sgmwu14/scratch/openfold_runs/run_subset_gpu32_h100_ep10/csv_logs/version_0/checkpoints/2-30000.ckpt
  /users/sgmwu14/scratch/openfold_runs/run_subset_gpu32_h100_ep10/csv_logs/version_0/checkpoints/3-40000.ckpt
  /users/sgmwu14/scratch/openfold_runs/run_subset_gpu32_h100_ep10/csv_logs/version_0/checkpoints/4-50000.ckpt
  /users/sgmwu14/scratch/openfold_runs/run_subset_gpu32_h100_ep10/csv_logs/version_0/checkpoints/5-60000.ckpt
  /users/sgmwu14/scratch/openfold_runs/run_subset_gpu32_h100_ep10/csv_logs/version_0/checkpoints/6-70000.ckpt
  /users/sgmwu14/scratch/openfold_runs/run_subset_gpu32_h100_ep10/csv_logs/version_0/checkpoints/7-80000.ckpt

  /users/sgmwu14/scratch/openfold_runs/run_subset_gpu32_h100_ep15/csv_logs/version_0/checkpoints/8-90000.ckpt
  /users/sgmwu14/scratch/openfold_runs/run_subset_gpu32_h100_ep15/csv_logs/version_0/checkpoints/9-100000.ckpt
  /users/sgmwu14/scratch/openfold_runs/run_subset_gpu32_h100_ep15/csv_logs/version_0/checkpoints/10-110000.ckpt
  /users/sgmwu14/scratch/openfold_runs/run_subset_gpu32_h100_ep15/csv_logs/version_0/checkpoints/11-120000.ckpt
  /users/sgmwu14/scratch/openfold_runs/run_subset_gpu32_h100_ep15/csv_logs/version_0/checkpoints/12-130000.ckpt
  /users/sgmwu14/scratch/openfold_runs/run_subset_gpu32_h100_ep15/csv_logs/version_0/checkpoints/13-140000.ckpt
  /users/sgmwu14/scratch/openfold_runs/run_subset_gpu32_h100_ep15/csv_logs/version_0/checkpoints/14-150000.ckpt

  /users/sgmwu14/scratch/openfold_runs/run_subset_gpu32_h100_ep21/csv_logs/version_0/checkpoints/15-160000.ckpt
  /users/sgmwu14/scratch/openfold_runs/run_subset_gpu32_h100_ep21/csv_logs/version_0/checkpoints/16-170000.ckpt
  /users/sgmwu14/scratch/openfold_runs/run_subset_gpu32_h100_ep21/csv_logs/version_0/checkpoints/17-180000.ckpt
  /users/sgmwu14/scratch/openfold_runs/run_subset_gpu32_h100_ep21/csv_logs/version_0/checkpoints/18-190000.ckpt
  /users/sgmwu14/scratch/openfold_runs/run_subset_gpu32_h100_ep21/csv_logs/version_0/checkpoints/19-200000.ckpt
  /users/sgmwu14/scratch/openfold_runs/run_subset_gpu32_h100_ep21/csv_logs/version_0/checkpoints/20-210000.ckpt

  /users/sgmwu14/scratch/openfold_runs/run_subset_gpu32_h100_ep27/csv_logs/version_0/checkpoints/21-220000.ckpt
  /users/sgmwu14/scratch/openfold_runs/run_subset_gpu32_h100_ep27/csv_logs/version_0/checkpoints/22-230000.ckpt
  /users/sgmwu14/scratch/openfold_runs/run_subset_gpu32_h100_ep27/csv_logs/version_0/checkpoints/23-240000.ckpt
  /users/sgmwu14/scratch/openfold_runs/run_subset_gpu32_h100_ep27/csv_logs/version_0/checkpoints/24-250000.ckpt
  /users/sgmwu14/scratch/openfold_runs/run_subset_gpu32_h100_ep27/csv_logs/version_0/checkpoints/25-260000.ckpt
  /users/sgmwu14/scratch/openfold_runs/run_subset_gpu32_h100_ep27/csv_logs/version_0/checkpoints/26-270000.ckpt

  /users/sgmwu14/scratch/openfold_runs/run_subset_gpu32_h100_ep33/csv_logs/version_0/checkpoints/27-280000.ckpt
  /users/sgmwu14/scratch/openfold_runs/run_subset_gpu32_h100_ep33/csv_logs/version_0/checkpoints/28-290000.ckpt
  /users/sgmwu14/scratch/openfold_runs/run_subset_gpu32_h100_ep33/csv_logs/version_0/checkpoints/29-300000.ckpt
  /users/sgmwu14/scratch/openfold_runs/run_subset_gpu32_h100_ep33/csv_logs/version_0/checkpoints/30-310000.ckpt
  /users/sgmwu14/scratch/openfold_runs/run_subset_gpu32_h100_ep33/csv_logs/version_0/checkpoints/31-320000.ckpt
  /users/sgmwu14/scratch/openfold_runs/run_subset_gpu32_h100_ep33/csv_logs/version_0/checkpoints/32-330000.ckpt
)

########################################
# Main loop
########################################
for MY_CKPT in "${CKPTS[@]}"; do
    if [[ ! -f "$MY_CKPT" ]]; then
        echo "[$(date)] ERROR: checkpoint not found: $MY_CKPT"
        exit 1
    fi

    CKPT_FILE=$(basename "$MY_CKPT")   # e.g. 27-280000.ckpt
    CKPT_TAG=${CKPT_FILE%%-*}          # e.g. 27

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