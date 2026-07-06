#!/bin/bash -l
#SBATCH -J pack_of_testset
#SBATCH -p long
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=1-00:00:00
#SBATCH --output=/users/%u/COMP702_BeyondAF/logs/slurm/pack_of_testset_%j.log
#SBATCH --error=/users/%u/COMP702_BeyondAF/logs/slurm/pack_of_testset_%j.err

set -euo pipefail

source "$HOME/COMP702_BeyondAF/config/comp702_paths.sh"

SRC="${OPENFOLD_PACK_SRC:-$OPENFOLD_GPU_DATA_ROOT}"
DST="${OPENFOLD_PACK_DST:-$COMP702_SCRATCH_ROOT/openfold_transfer}"
TAR="$DST/openfold_testset_$(date +%F_%H%M%S).tar.gz"

mkdir -p "$DST"

TESTSET_REL="alignment_testset"
MMCIF_REL="pdb_data/mmcif/mmcif_files"
CACHES_REL="pdb_data/data_caches"
OBSOLETE_REL="pdb_data/obsolete.dat"

echo "Node: $(hostname)"
echo "Start: $(date)"
echo "Source: $SRC"
echo "Archive: $TAR"

test -d "$SRC" || { echo "Missing source directory: $SRC"; exit 1; }

INCLUDES=()

for item in "$TESTSET_REL" "$MMCIF_REL" "$CACHES_REL"; do
  if [ -e "$SRC/$item" ]; then
    INCLUDES+=( "$item" )
  else
    echo "WARNING: missing item, skipping: $SRC/$item"
  fi
done

if [ -f "$SRC/$OBSOLETE_REL" ]; then
  INCLUDES+=( "$OBSOLETE_REL" )
fi

if [ "${#INCLUDES[@]}" -eq 0 ]; then
  echo "No packable testset items found under $SRC"
  exit 1
fi

if command -v pigz >/dev/null 2>&1; then
  tar -C "$SRC" -cf - "${INCLUDES[@]}" | pigz -p "${SLURM_CPUS_PER_TASK:-16}" > "$TAR"
else
  tar -C "$SRC" -czf "$TAR" "${INCLUDES[@]}"
fi

ls -lh "$TAR"
echo "Done: $(date)"
