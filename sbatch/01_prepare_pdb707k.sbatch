#!/bin/bash
#SBATCH --job-name=prep_pdb707k
#SBATCH --partition=nodes
#SBATCH --time=00:20:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --output=/users/sgsrehm1/COMP702_BeyondAF/logs/%x_%j.out
#SBATCH --error=/users/sgsrehm1/COMP702_BeyondAF/logs/%x_%j.err

set -euo pipefail

module purge
module load miniforge3/25.3.0-python3.12.10
eval "$(conda shell.bash hook)"
conda activate /users/sgsrehm1/COMP702_BeyondAF/envs/comp702-data

python /users/sgsrehm1/COMP702_BeyondAF/scripts/01_prepare_pdb707k.py
