#!/bin/bash
#SBATCH --job-name=2olo_of_proj
#SBATCH --partition=gpu-l40s
#SBATCH --gres=gpu:1
#SBATCH --time=03:00:00
#SBATCH --mem=64G
#SBATCH --output=/users/sgsrehm1/COMP702_BeyondAF/logs/%x_%j.out
#SBATCH --error=/users/sgsrehm1/COMP702_BeyondAF/logs/%x_%j.err

set -euo pipefail

module purge
module load miniforge3/25.3.0-python3.12.10

eval "$(conda shell.bash hook)"
conda activate /users/sgsrehm1/COMP702_BeyondAF/envs/openfold_env

cd /users/sgsrehm1/COMP702_BeyondAF/code/COMP390_code

export PYTHONPATH=$(pwd)/openfold:${PYTHONPATH:-}
export TRITON_CACHE_DIR=/users/sgsrehm1/fastscratch/COMP702_BeyondAF/triton_cache
mkdir -p "$TRITON_CACHE_DIR"

mkdir -p demo_run_outputs_2olo/openfold_official
mkdir -p demo_run_outputs_2olo/your_model_ep33

echo "Running on node: $(hostname)"
echo "Job ID: ${SLURM_JOB_ID}"
echo "CUDA visible devices: ${CUDA_VISIBLE_DEVICES:-not_set}"
nvidia-smi

echo "===== OpenFold official checkpoint run ====="

python openfold/run_pretrained_openfold.py \
  demo_single_2olo/fasta_dir \
  demo_single_2olo/taskcif \
  --output_dir demo_run_outputs_2olo/openfold_official \
  --use_precomputed_alignments demo_single_2olo/alignments_cf \
  --config_preset model_1 \
  --openfold_checkpoint_path openfold/resources/openfold_params/finetuning_2.pt \
  --model_device cuda:0 \
  --data_random_seed 101

echo "===== Project-trained epoch-33 checkpoint run ====="

python openfold/run_pretrained_openfold.py \
  demo_single_2olo/fasta_dir \
  demo_single_2olo/taskcif \
  --output_dir demo_run_outputs_2olo/your_model_ep33 \
  --use_precomputed_alignments demo_single_2olo/alignments_cf \
  --config_preset model_1 \
  --openfold_checkpoint_path checkpoints/run_subset_gpu32_h100_ep33/csv_logs/version_0/checkpoints/32-330000.ckpt \
  --model_device cuda:0 \
  --data_random_seed 101

echo "All 2OLO demo runs completed."
