#!/bin/bash -l
#SBATCH --job-name=of_ckpt_verify
#SBATCH --partition=gpu-l40s-low
#SBATCH --gres=gpu:l40s:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=01:00:00

set -euo pipefail

CKPT="${1:?checkpoint required}"

REPO=/mnt/fastscratch/users/sgsrehm1/COMP702_pdbclean_pipeline
BASE="$REPO/outputs/openfold_training/20260101/pdbclean-dedup-v1"
VIEW="$BASE/training_view"
ENV=/mnt/fastscratch/users/sgsrehm1/envs/openfold_train

. "$ENV/etc/conda/activate.d/10-cache-dirs.sh"
export PYTHONPATH="$REPO/src${PYTHONPATH:+:$PYTHONPATH}"

echo "host   : $(hostname)"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader

cd /mnt/fastscratch/users/sgsrehm1/openfold_src

"$ENV/bin/python" "$REPO/scripts/openfold_training/verify_checkpoint_and_memory.py" \
    --checkpoint "$CKPT" \
    --mmcif-dir /users/sgsrehm1/scratch/COMP702_openfold_20260101/pdbclean-dedup-v1/mmcif \
    --msa-store /mnt/fastscratch/users/sgsrehm1/COMP702_openfold_msa/20260101/pdbclean-dedup-v1/msas \
    --alignment-index "$BASE/alignment_index/pdbclean_dedup_v1.index" \
    --projection-index "$VIEW/projection_index.json" \
    --train-filter "$VIEW/train_filter_view.txt" \
    --steps 2 \
    --report "$VIEW/checkpoint_reload_and_memory.json"
