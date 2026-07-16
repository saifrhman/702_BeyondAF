#!/bin/bash
#SBATCH --job-name=2olo_af
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

export PYTHONPATH=$(pwd)/openfold:$PYTHONPATH
export TRITON_CACHE_DIR=/users/sgsrehm1/fastscratch/COMP702_BeyondAF/triton_cache
mkdir -p "$TRITON_CACHE_DIR"

mkdir -p demo_run_outputs_2olo/alphafold

echo "Running on node: $(hostname)"
echo "Job ID: ${SLURM_JOB_ID}"
echo "CUDA visible devices: ${CUDA_VISIBLE_DEVICES:-not_set}"
nvidia-smi

python openfold/run_pretrained_openfold.py \
  demo_single_2olo/fasta_dir \
  demo_single_2olo/taskcif \
  --output_dir demo_run_outputs_2olo/alphafold \
  --use_precomputed_alignments demo_single_2olo/alignments_cf \
  --config_preset model_1 \
  --jax_param_path openfold/resources/params/params_model_1.npz \
  --model_device cuda:0 \
  --data_random_seed 101
