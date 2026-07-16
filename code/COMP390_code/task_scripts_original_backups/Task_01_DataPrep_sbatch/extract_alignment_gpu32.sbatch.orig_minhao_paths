#!/bin/bash -l
#SBATCH -J extract_align
#SBATCH -p gpu-h100
#SBATCH -N 1
#SBATCH -w gpu32
#SBATCH --gres=gpu:h100:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=80G
#SBATCH --time=1-00:00:00
#SBATCH --output=/users/%u/fastscratch/extract_align_%j.log
#SBATCH --error=/users/%u/fastscratch/extract_align_%j.err

set -euo pipefail

echo "Starting alignment extraction on $(hostname)"

SRC=/users/sgmwu14/scratch/alignment_data_node037_20260210_083640.tar.zst
DEST_DIR=/users/sgmwu14/localscratch/alignment_openfold_full_20260210

mkdir -p "$DEST_DIR"

echo "Copying archive..."
cp "$SRC" "$DEST_DIR/"

cd "$DEST_DIR"

echo "Extracting archive..."
tar --use-compress-program=unzstd -xvf alignment_data_node037_20260210_083640.tar.zst

echo "Extraction finished"
ls -lh