#!/bin/bash
#SBATCH --job-name=prep_pdb707k
#SBATCH --partition=nodes
#SBATCH --time=00:20:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# Portability: every path below comes from here, derived from $USER.
# Override any of them in the environment; nothing names an individual.
source "${PDBCLEAN_ENV_FILE:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/../config/pdbclean/pipeline_env.sh}"

set -euo pipefail

module purge
module load miniforge3/25.3.0-python3.12.10
eval "$(conda shell.bash hook)"
conda activate ${COMP702_ROOT}/envs/comp702-data

python ${COMP702_ROOT}/scripts/01_prepare_pdb707k.py
