#!/bin/bash -l
#SBATCH -J pack_of_subset
#SBATCH -p long
#SBATCH -N 1
#SBATCH -w node037
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=1-00:00:00
#SBATCH --output=/users/%u/scratch/pack_of_subset_%j.log
#SBATCH --error=/users/%u/scratch/pack_of_subset_%j.err

set -euo pipefail

echo "======================================"
echo "Node: $(hostname)"
echo "Start: $(date)"
echo "======================================"

SRC=/users/$USER/localscratch
DST=/users/$USER/scratch/openfold_transfer
TAR=$DST/openfold_subset_node037.tar.gz

mkdir -p "$DST"

echo "Packing: alignment_subset + lists + pdb_data"
echo "From: $SRC"
echo "To:   $TAR"

# 使用 pigz（如果有）以加速
if command -v pigz >/dev/null 2>&1; then
  echo "Using pigz compression"
  tar -C "$SRC" -cf - alignment_subset lists pdb_data | pigz -p 16 > "$TAR"
else
  echo "Using gzip compression"
  tar -C "$SRC" -czf "$TAR" alignment_subset lists pdb_data
fi

echo "--------------------------------------"
echo "Package created:"
ls -lh "$TAR"
echo "Finished: $(date)"
echo "======================================"