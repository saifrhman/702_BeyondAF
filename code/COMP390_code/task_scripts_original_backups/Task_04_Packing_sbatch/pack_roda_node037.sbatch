#!/bin/bash -l
#SBATCH -J pack_roda
#SBATCH -p nodes
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=8
#SBATCH --time=1-00:00:00
#SBATCH -w node037
#SBATCH --output=/users/%u/fastscratch/pack_roda_%j.log

set -euo pipefail

echo "Node: $(hostname)"
echo "Start: $(date)"
echo "USER=$USER"

SRC="/users/$USER/localscratch/alignment_data"
DST="/users/$USER/scratch"
ARCHIVE="alignment_data_node037_$(date +%Y%m%d_%H%M%S).tar.zst"

mkdir -p "$DST"

# --- sanity checks ---
if [[ ! -d "$SRC" ]]; then
  echo "ERROR: $SRC does not exist on this node."
  exit 1
fi

echo "Checking size & file count before packing..."
du -sh "$SRC" || true
find "$SRC" -type f | wc -l || true

# Make sure it's not empty
if [[ "$(find "$SRC" -type f | head -n 1 | wc -l)" -eq 0 ]]; then
  echo "ERROR: $SRC contains no files on this node (would pack an empty archive)."
  exit 1
fi

cd "$(dirname "$SRC")"

echo "Compressing (tar + zstd, threads=8) ..."
tar -I "zstd -T8" -cf "$ARCHIVE" "$(basename "$SRC")"

echo "Archive size:"
ls -lh "$ARCHIVE"

echo "Moving archive to scratch..."
mv "$ARCHIVE" "$DST/"

echo "Finished: $(date)"
echo "Saved to: $DST/$ARCHIVE"