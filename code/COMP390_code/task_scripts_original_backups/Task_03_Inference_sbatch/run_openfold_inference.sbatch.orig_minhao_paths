#!/bin/bash -l
#SBATCH -J of_infer3
#SBATCH -p gpu-h100
#SBATCH -N 1
#SBATCH -w gpu32
#SBATCH --gres=gpu:h100:1
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=3-00:00:00
#SBATCH --output=/users/%u/fastscratch/of_infer3_%j.log
#SBATCH --error=/users/%u/fastscratch/of_infer3_%j.err

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

########################################
# Paths
########################################

FASTA_DIR=/users/sgmwu14/scratch/openfold_inference/fasta_dir
CIF_DIR=/users/sgmwu14/scratch/taskcif
ALIGN_DIR=/users/sgmwu14/scratch/openfold_inference/alignments_cf

AF_OUT=/users/sgmwu14/scratch/openfold_inference/alphafold_output
OF_OUT=/users/sgmwu14/scratch/openfold_inference/openfold_output
MY_OUT=/users/sgmwu14/scratch/openfold_inference/your_model_output_ep21

mkdir -p $AF_OUT
mkdir -p $OF_OUT
mkdir -p $MY_OUT

########################################
# 1️⃣ AlphaFold parameters
########################################

echo "===== Running AlphaFold params ====="

python3 /users/sgmwu14/openfold/run_pretrained_openfold.py \
  $FASTA_DIR \
  $CIF_DIR \
  --output_dir $AF_OUT \
  --use_precomputed_alignments $ALIGN_DIR \
  --config_preset model_1 \
  --jax_param_path /users/sgmwu14/openfold/openfold/resources/params/params_model_1.npz \
  --model_device cuda:0

########################################
# 2️⃣ OpenFold official model
########################################

echo "===== Running OpenFold official ====="

python3 /users/sgmwu14/openfold/run_pretrained_openfold.py \
  $FASTA_DIR \
  $CIF_DIR \
  --output_dir $OF_OUT \
  --use_precomputed_alignments $ALIGN_DIR \
  --config_preset model_1 \
  --openfold_checkpoint_path /users/sgmwu14/openfold/openfold/resources/openfold_params/finetuning_2.pt \
  --model_device cuda:0

########################################
# 3️⃣ Your trained model (ep21)
########################################

echo "===== Running Your Model ep21 ====="

python3 /users/sgmwu14/openfold/run_pretrained_openfold.py \
  $FASTA_DIR \
  $CIF_DIR \
  --output_dir $MY_OUT \
  --use_precomputed_alignments $ALIGN_DIR \
  --config_preset model_1 \
  --openfold_checkpoint_path /users/sgmwu14/scratch/openfold_runs/run_subset_gpu32_h100_ep21/csv_logs/version_0/checkpoints/20-210000.ckpt \
  --model_device cuda:0

echo "===== ALL INFERENCE FINISHED ====="