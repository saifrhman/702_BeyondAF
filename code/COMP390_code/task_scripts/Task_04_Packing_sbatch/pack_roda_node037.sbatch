#!/bin/bash -l
#SBATCH -J pack_roda
#SBATCH -p nodes
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=8
#SBATCH --time=1-00:00:00
#SBATCH --output=/users/%u/COMP702_BeyondAF/logs/slurm/pack_roda_%j.log
#SBATCH --error=/users/%u/COMP702_BeyondAF/logs/slurm/pack_roda_%j.err

set -euo pipefail

source "$HOME/COMP702_BeyondAF/config/comp702_paths.sh"

SRC="${RODA_PACK_SRC:-$OPENFOLD_GPU_LOCAL_BASE/alignment_data}"
DST="${RODA_PACK_DST:-$COMP702_SCRATCH_ROOT}"
ARCHIVE="alignment_data_$(date +%Y%m%d_%H%M%S).tar.zst"

mkdir -p "$DST"

echo "Node: $(hostname)"
echo "Start: $(date)"
echo "USER=$USER"
echo "Source: $SRC"
echo "Destination: $DST/$ARCHIVE"

if [[ ! -d "$SRC" ]]; then
  echo "ERROR: source directory does not exist: $SRC"
  exit 1
fi

echo "Checking source size and file count..."
du -sh "$SRC" || true
find "$SRC" -type f | wc -l || true

if [[ "$(find "$SRC" -type f | head -n 1 | wc -l)" -eq 0 ]]; then
  echo "ERROR: source directory contains no files."
  exit 1
fi

cd "$(dirname "$SRC")"

if command -v zstd >/dev/null 2>&1; then
  tar -I "zstd -T${SLURM_CPUS_PER_TASK:-8}" -cf "$ARCHIVE" "$(basename "$SRC")"
else
  echo "ERROR: zstd is not available on this node."
  exit 1
fi

ls -lh "$ARCHIVE"
mv "$ARCHIVE" "$DST/"

echo "Finished: $(date)"
echo "Saved to: $DST/$ARCHIVE"
