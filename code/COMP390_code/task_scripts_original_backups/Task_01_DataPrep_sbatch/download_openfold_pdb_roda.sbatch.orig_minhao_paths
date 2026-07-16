#!/bin/bash -l
#SBATCH -J of_roda_dl
#SBATCH -p long
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=7-00:00:00
#SBATCH --output=/users/%u/fastscratch/of_roda_dl_%j.log

set -euo pipefail

echo "Job $SLURM_JOB_ID started on $(hostname) at $(date)"
echo "USER=$USER"

# 你现在 aws 在 conda 里就能用，不要 module load awscli
which aws
aws --version

SRC="s3://openfold/pdb"
DST="/mnt/scratch/users/$USER/alignment_data/alignment_dir_roda"   # 等价于 ~/scratch/...

mkdir -p "$DST"

echo "Syncing from: $SRC"
echo "Syncing to:   $DST"

# 建议加 --only-show-errors，日志不会爆炸
aws s3 sync \
  "$SRC" \
  "$DST" \
  --no-sign-request \
  --only-show-errors

echo "Job $SLURM_JOB_ID finished at $(date)"