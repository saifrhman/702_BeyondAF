#!/bin/bash -l
#SBATCH -J pack_of_testset
#SBATCH -p long
#SBATCH -N 1
#SBATCH -w node037
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=1-00:00:00
#SBATCH --output=/users/%u/fastscratch/pack_of_testset_%j.log
#SBATCH --error=/users/%u/fastscratch/pack_of_testset_%j.err

set -euo pipefail

SRC=/users/$USER/localscratch
DST=/users/$USER/fastscratch/openfold_transfer
TAR=$DST/openfold_testset_$(date +%F_%H%M%S).tar.gz

mkdir -p "$DST"

# 你 testset 的根目录（你已经用它生成了 all-seqs.fasta）
TESTSET_REL=alignment_testset

# 训练/模板需要的 mmcif（你当前训练用的是这个路径）
MMCIF_REL=pdb_data/mmcif/mmcif_files

# 训练会用到的 cache（你是下载好的，直接打包即可）
CACHES_REL=pdb_data/data_caches

# obsolete（有就打，没有就跳过）
OBSOLETE_REL=pdb_data/obsolete.dat

echo "Packing from $SRC ..."
echo "Including:"
echo "  - $TESTSET_REL (alignments + alignment_dbs + all-seqs*.txt/fasta if present)"
echo "  - $MMCIF_REL"
echo "  - $CACHES_REL"
echo "  - $OBSOLETE_REL (if exists)"
echo

# 组装要打包的文件清单（相对 SRC）
INCLUDES=(
  "$TESTSET_REL"
  "$MMCIF_REL"
  "$CACHES_REL"
)

# obsolete 可选
if [ -f "$SRC/$OBSOLETE_REL" ]; then
  INCLUDES+=( "$OBSOLETE_REL" )
fi

# pigz 可用就用 pigz（更快）
if command -v pigz >/dev/null 2>&1; then
  tar -C "$SRC" -cf - "${INCLUDES[@]}" | pigz -p ${SLURM_CPUS_PER_TASK:-16} > "$TAR"
else
  tar -C "$SRC" -czf "$TAR" "${INCLUDES[@]}"
fi

ls -lh "$TAR"
echo "Done: $(date)"