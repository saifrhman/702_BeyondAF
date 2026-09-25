#!/bin/bash
#SBATCH --job-name=2olo_of_proj
#SBATCH --partition=gpu-l40s
#SBATCH --gres=gpu:1
#SBATCH --time=03:00:00
#SBATCH --mem=64G
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# Portability: every path below comes from here, derived from $USER.
# Override any of them in the environment; nothing names an individual.
# Slurm copies this script into its spool directory before running it, so
# ${BASH_SOURCE[0]} does not point into the repository under sbatch. Try the
# explicit override, the repo root the submitter exported, the directory the
# job was submitted from, and only then this file's own location.
for _pdbclean_env in "${PDBCLEAN_ENV_FILE:-}" \
                     "${PDBCLEAN_REPO_ROOT:-}/config/pdbclean/pipeline_env.sh" \
                     "${SLURM_SUBMIT_DIR:-}/config/pdbclean/pipeline_env.sh" \
                     "$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd)/../config/pdbclean/pipeline_env.sh"; do
    [[ -f "$_pdbclean_env" ]] && { source "$_pdbclean_env"; break; }
done
unset _pdbclean_env

if [[ -z "${PDBCLEAN_REPO_ROOT:-}" ]]; then
    echo "cannot locate config/pdbclean/pipeline_env.sh; set PDBCLEAN_ENV_FILE" >&2
    exit 1
fi

set -euo pipefail

module purge
module load miniforge3/25.3.0-python3.12.10

eval "$(conda shell.bash hook)"
conda activate ${COMP702_ROOT}/envs/openfold_env

cd ${COMP702_ROOT}/code/COMP390_code

export PYTHONPATH=$(pwd)/openfold:${PYTHONPATH:-}
export TRITON_CACHE_DIR=${TRITON_CACHE_DIR}
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
