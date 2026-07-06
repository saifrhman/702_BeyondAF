#!/bin/bash -l
#SBATCH -J pack_of_subset
#SBATCH -p long
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=1-00:00:00
#SBATCH --output=/users/%u/COMP702_BeyondAF/logs/slurm/pack_of_subset_to_scratch_%j.log
#SBATCH --error=/users/%u/COMP702_BeyondAF/logs/slurm/pack_of_subset_to_scratch_%j.err

set -euo pipefail

source "$HOME/COMP702_BeyondAF/config/comp702_paths.sh"

SRC="${OPENFOLD_PACK_SRC:-$OPENFOLD_GPU_DATA_ROOT}"
DST="${OPENFOLD_PACK_DST:-$COMP702_SCRATCH_ROOT/openfold_transfer}"
TAR="$DST/openfold_subset.tar.gz"

mkdir -p "$DST"

echo "======================================"
echo "Node: $(hostname)"
echo "Start: $(date)"
echo "Source: $SRC"
echo "Archive: $TAR"
echo "======================================"

test -d "$SRC" || { echo "Missing source directory: $SRC"; exit 1; }

INCLUDES=()
for item in alignment_subset lists pdb_data; do
  if [ -e "$SRC/$item" ]; then
    INCLUDES+=( "$item" )
  else
    echo "WARNING: missing item, skipping: $SRC/$item"
  fi
done

if [ "${#INCLUDES[@]}" -eq 0 ]; then
  echo "No packable items found under $SRC"
  exit 1
fi

if command -v pigz >/dev/null 2>&1; then
  echo "Using pigz compression"
  tar -C "$SRC" -cf - "${INCLUDES[@]}" | pigz -p "${SLURM_CPUS_PER_TASK:-16}" > "$TAR"
else
  echo "Using gzip compression"
  tar -C "$SRC" -czf "$TAR" "${INCLUDES[@]}"
fi

echo "Package created:"
ls -lh "$TAR"
echo "Finished: $(date)"
echo "======================================"
