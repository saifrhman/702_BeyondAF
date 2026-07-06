#!/bin/bash -l
#SBATCH -J extract_align
#SBATCH -p long
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=80G
#SBATCH --time=1-00:00:00
#SBATCH --output=/users/%u/COMP702_BeyondAF/logs/slurm/extract_align_%j.log
#SBATCH --error=/users/%u/COMP702_BeyondAF/logs/slurm/extract_align_%j.err

set -euo pipefail

source "$HOME/COMP702_BeyondAF/config/comp702_paths.sh"

echo "Starting alignment extraction on $(hostname) at $(date)"
echo "USER=$USER"

SRC="$ALIGNMENT_ARCHIVE"
DEST_DIR="$ALIGNMENT_EXTRACT_DIR"

mkdir -p "$DEST_DIR"

echo "Archive source:      $SRC"
echo "Extraction location: $DEST_DIR"

if [ ! -f "$SRC" ]; then
  echo "ERROR: alignment archive not found: $SRC"
  exit 1
fi

echo "Copying archive..."
cp "$SRC" "$DEST_DIR/"

cd "$DEST_DIR"

ARCHIVE_NAME="$(basename "$SRC")"

echo "Extracting archive: $ARCHIVE_NAME"
tar --use-compress-program=unzstd -xvf "$ARCHIVE_NAME"

echo "Extraction finished at $(date)"
ls -lh
