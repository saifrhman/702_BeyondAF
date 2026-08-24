#!/bin/bash -l
#SBATCH --job-name=of26_msa_audit
#SBATCH --partition=nodes
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/users/%u/fastscratch/COMP702_pdbclean_pipeline/outputs/openfold_training/20260101/pdbclean-dedup-v1/msa/audit/msa_audit_%j.out
#SBATCH --error=/users/%u/fastscratch/COMP702_pdbclean_pipeline/outputs/openfold_training/20260101/pdbclean-dedup-v1/msa/audit/msa_audit_%j.err

set -euo pipefail

ROOT="$HOME/fastscratch/COMP702_pdbclean_pipeline"
HOME_REPO="$HOME/COMP702_BeyondAF"

TRAIN="$ROOT/outputs/openfold_training/20260101/pdbclean-dedup-v1"
MSA="$TRAIN/msa"

OPENFOLD="$HOME_REPO/code/COMP390_code/openfold"
OPENFOLD_ENV="$HOME_REPO/envs/openfold_env"

export PYTHONUNBUFFERED=1
export PYTHONDONTWRITEBYTECODE=1
export XDG_CACHE_HOME="$HOME/fastscratch/.cache"

cd "$ROOT"

echo "============================================================"
echo "OPENFOLD / MSA PIPELINE AUDIT"
echo "Job:     ${SLURM_JOB_ID}"
echo "Node:    $(hostname)"
echo "Started: $(date --iso-8601=seconds)"
echo "============================================================"

echo
echo "===== 1. FROZEN MSA QUERY POPULATION ====="

if [[ -f "$MSA/_QUERIES_SUCCESS" ]]; then
    cat "$MSA/_QUERIES_SUCCESS"
else
    echo "MISSING: $MSA/_QUERIES_SUCCESS"
fi

echo
echo "Query files:"
ls -lh \
    "$MSA/msa_queries.parquet" \
    "$MSA/msa_queries.fasta" \
    "$MSA/msa_query_manifest.json" \
    2>/dev/null || true

echo
echo "Query manifest:"
cat "$MSA/msa_query_manifest.json" 2>/dev/null || true


echo
echo "===== 2. HISTORICAL COMP702 PATH CONFIG ====="

if [[ -f "$HOME_REPO/config/comp702_paths.sh" ]]; then
    grep -nEi \
        'openfold|align|uniref|mgnify|bfd|uniclust|pdb70|database|hhblits|jackhmmer' \
        "$HOME_REPO/config/comp702_paths.sh" \
        || true
else
    echo "MISSING comp702_paths.sh"
fi


echo
echo "===== 3. HISTORICAL TRAINING / ALIGNMENT SCRIPTS ====="

find \
    "$HOME_REPO/code/COMP390_code/task_scripts" \
    "$HOME_REPO/code/COMP390_code/openfold/scripts" \
    -type f \
    \( -name '*.sh' -o -name '*.sbatch' -o -name '*.py' \) \
    -print0 2>/dev/null \
| xargs -0 grep -nEi \
    'alignment|jackhmmer|hhblits|hhsearch|uniref|mgnify|bfd|uniclust|pdb70|pdb_seqres|alignment_db' \
    2>/dev/null \
| head -n 500 \
|| true


echo
echo "===== 4. OPENFOLD DATA-PIPELINE SEARCH TOOLS ====="

grep -RniE \
    'jackhmmer|hhblits|hhsearch|uniref90|mgnify|bfd|uniclust30|pdb70|pdb_seqres' \
    "$OPENFOLD/openfold/data" \
    "$OPENFOLD/scripts" \
    2>/dev/null \
| head -n 500 \
|| true


echo
echo "===== 5. AVAILABLE BINARIES ====="

echo "System/PATH:"
for TOOL in \
    jackhmmer \
    hmmbuild \
    hmmsearch \
    hhblits \
    hhsearch \
    mmseqs
do
    printf '%-12s ' "$TOOL"
    command -v "$TOOL" || echo "NOT FOUND"
done

echo
echo "OpenFold environment:"
if [[ -x "$OPENFOLD_ENV/bin/python" ]]; then
    for TOOL in \
        jackhmmer \
        hmmbuild \
        hmmsearch \
        hhblits \
        hhsearch \
        mmseqs
    do
        if [[ -x "$OPENFOLD_ENV/bin/$TOOL" ]]; then
            echo "$TOOL -> $OPENFOLD_ENV/bin/$TOOL"
        else
            echo "$TOOL -> NOT IN OPENFOLD ENV"
        fi
    done
else
    echo "OpenFold environment missing: $OPENFOLD_ENV"
fi


echo
echo "===== 6. KNOWN OLD ALIGNMENT LOCATIONS ====="

for P in \
    "$HOME/scratch/openfold_data" \
    "$HOME/fastscratch/openfold_data" \
    "$HOME_REPO/openfold_data" \
    "$HOME_REPO/data" \
    "/mnt/data1/users/$USER/openfold_data" \
    "/mnt/data2/users/$USER/openfold_data"
do
    echo
    echo "--- $P ---"

    if [[ -e "$P" ]]; then
        du -sh "$P" 2>/dev/null || true
        find "$P" -maxdepth 3 -type d 2>/dev/null | head -n 100
    else
        echo "NOT PRESENT"
    fi
done


echo
echo "===== 7. TARGETED DATABASE DISCOVERY ====="

SEARCH_ROOTS=(
    "$HOME/scratch"
    "$HOME/fastscratch"
    "/mnt/data1/users/$USER"
    "/mnt/data2/users/$USER"
)

for BASE in "${SEARCH_ROOTS[@]}"; do
    echo
    echo "--- searching $BASE ---"

    if [[ ! -d "$BASE" ]]; then
        echo "NOT PRESENT"
        continue
    fi

    find "$BASE" \
        -maxdepth 5 \
        \( \
            -iname '*uniref*' \
            -o -iname '*mgnify*' \
            -o -iname '*bfd*' \
            -o -iname '*uniclust*' \
            -o -iname '*pdb70*' \
            -o -iname '*pdb_seqres*' \
            -o -iname '*alignment_db*' \
            -o -iname '*alignment*dbs*' \
        \) \
        -print \
        2>/dev/null \
        | head -n 300
done


echo
echo "===== 8. EXISTING ALIGNMENT INDEX / DB FILES ====="

for BASE in "${SEARCH_ROOTS[@]}"; do
    if [[ ! -d "$BASE" ]]; then
        continue
    fi

    find "$BASE" \
        -maxdepth 6 \
        -type f \
        \( \
            -name 'alignment_db.index' \
            -o -name '*alignment*index*.json' \
            -o -name '*alignment*.index' \
        \) \
        -print \
        2>/dev/null \
        | head -n 200
done


echo
echo "===== 9. OPENFOLD SCRIPT INVENTORY ====="

find "$OPENFOLD/scripts" \
    -maxdepth 2 \
    -type f \
    -printf '%f\n' \
    2>/dev/null \
    | sort \
    | grep -Ei \
        'align|msa|precompute|database|cache' \
    || true


echo
echo "============================================================"
echo "MSA PIPELINE AUDIT COMPLETE"
echo "Completed: $(date --iso-8601=seconds)"
echo "============================================================"
