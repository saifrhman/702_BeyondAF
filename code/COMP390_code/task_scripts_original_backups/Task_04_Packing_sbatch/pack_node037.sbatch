#!/bin/bash -l
#SBATCH -J pack_of_subset
#SBATCH -p long
#SBATCH -N 1
#SBATCH -w node037
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=1-00:00:00
#SBATCH --output=/users/%u/fastscratch/pack_of_subset_%j.log
#SBATCH --error=/users/%u/fastscratch/pack_of_subset_%j.err

set -euo pipefail

SRC=/users/$USER/localscratch
DST=/users/$USER/fastscratch/openfold_transfer
TAR=$DST/openfold_subset_node037.tar.gz

mkdir -p "$DST"

echo "Packing from $SRC ..."
# pigz 可用就用 pigz（更快）
if command -v pigz >/dev/null 2>&1; then
  tar -C "$SRC" -cf - pdb_data alignment_subset | pigz -p 16 > "$TAR"
else
  tar -C "$SRC" -czf "$TAR" pdb_data alignment_subset
fi

ls -lh "$TAR"
echo "Done: $(date)"