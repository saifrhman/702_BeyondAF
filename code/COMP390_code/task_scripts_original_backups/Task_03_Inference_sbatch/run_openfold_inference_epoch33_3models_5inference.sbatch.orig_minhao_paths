#!/bin/bash -l
#SBATCH -J of_infer_ep33_x5
#SBATCH -p gpu-h100
#SBATCH -N 1
#SBATCH -w gpu32
#SBATCH --gres=gpu:h100:1
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=1-00:00:00
#SBATCH --output=/users/%u/fastscratch/of_infer_ep33_x5_%j.log
#SBATCH --error=/users/%u/fastscratch/of_infer_ep33_x5_%j.err

set -euo pipefail

########################################
# Activate environment
########################################
set +u
source /mnt/data2/users/sgmwu14/tools/miniforge3/etc/profile.d/conda.sh
conda activate openfold_env
set -u

export PYTHONPATH=/users/sgmwu14/openfold:$PYTHONPATH
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}

export TRITON_CACHE_DIR=/users/sgmwu14/fastscratch/triton_cache
mkdir -p "$TRITON_CACHE_DIR"

########################################
# Paths
########################################
FASTA_DIR=/users/sgmwu14/scratch/openfold_inference/fasta_dir
CIF_DIR=/users/sgmwu14/scratch/taskcif
ALIGN_DIR=/users/sgmwu14/scratch/openfold_inference/alignments_cf

# 你当前要跑的 epoch
EPOCH_NUM=33

# 你的 checkpoint
MY_CKPT=/users/sgmwu14/scratch/openfold_runs/run_subset_gpu32_h100_ep33/csv_logs/version_0/checkpoints/27-280000.ckpt

# 总输出目录：按 epoch 命名
RUN_ROOT=/users/sgmwu14/scratch/openfold_inference_runs/epoch_${EPOCH_NUM}

# 三个模型目录：按 epoch / 模型命名
AF_ROOT=${RUN_ROOT}/alphafold
OF_ROOT=${RUN_ROOT}/openfold_official
MY_ROOT=${RUN_ROOT}/your_model_ep${EPOCH_NUM}

mkdir -p "$AF_ROOT" "$OF_ROOT" "$MY_ROOT"

########################################
# Repeat inference 5 times
########################################
for REP in 1 2 3 4 5; do
    SEED=$((100 + REP))

    echo "[$(date)] ========================================"
    echo "[$(date)] Starting repeat ${REP} with seed ${SEED}"
    echo "[$(date)] ========================================"

    AF_OUT=${AF_ROOT}/repeat_${REP}
    OF_OUT=${OF_ROOT}/repeat_${REP}
    MY_OUT=${MY_ROOT}/repeat_${REP}

    mkdir -p "$AF_OUT" "$OF_OUT" "$MY_OUT"

    ########################################
    # 1) AlphaFold parameters
    ########################################
    echo "[$(date)] Running AlphaFold params, repeat ${REP}"

    python3 /users/sgmwu14/openfold/run_pretrained_openfold.py \
      "$FASTA_DIR" \
      "$CIF_DIR" \
      --output_dir "$AF_OUT" \
      --use_precomputed_alignments "$ALIGN_DIR" \
      --config_preset model_1 \
      --jax_param_path /users/sgmwu14/openfold/openfold/resources/params/params_model_1.npz \
      --model_device cuda:0 \
      --data_random_seed "$SEED"

    ########################################
    # 2) OpenFold official model
    ########################################
    echo "[$(date)] Running OpenFold official, repeat ${REP}"

    python3 /users/sgmwu14/openfold/run_pretrained_openfold.py \
      "$FASTA_DIR" \
      "$CIF_DIR" \
      --output_dir "$OF_OUT" \
      --use_precomputed_alignments "$ALIGN_DIR" \
      --config_preset model_1 \
      --openfold_checkpoint_path /users/sgmwu14/openfold/openfold/resources/openfold_params/finetuning_2.pt \
      --model_device cuda:0 \
      --data_random_seed "$SEED"

    ########################################
    # 3) Your trained model
    ########################################
    echo "[$(date)] Running your model ep${EPOCH_NUM}, repeat ${REP}"

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

echo "[$(date)] ===== ALL 5 REPEATED INFERENCE RUNS FINISHED ====="