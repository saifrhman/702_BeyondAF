#!/bin/bash -l
#SBATCH --job-name=of_msa_prod
#SBATCH --partition=nodes
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=400G
#SBATCH --time=24:00:00

# Stage 17 production: MSAs for the OpenFold training set.
#
# Scientifically identical to run_msa_pilot.sbatch
# ------------------------------------------------
# The mmseqs command sequence, every parameter, and the index-free fallback
# below are copied from the pilot without alteration:
#
#   createdb -> search (3 iterations) -> mvdb/lndb profile handoff
#            -> expandaln -> align -> filterresult -> result2msa -> unpackdb
#
# The pilot is the control that proves this protocol; production must not
# diverge from it, so a divergence is a bug rather than a tuning opportunity.
# tests/pdbclean/test_msa_production.py asserts the two remain in step.
#
# What differs from the pilot, and why
# ------------------------------------
# 1. It shards. One array task handles one contiguous slice of the unique
#    sequence set, so 142,056 searches fit in a bounded number of jobs.
#
# 2. It is keyed by sequence, not by chain. The query set is
#    `unique_sequences.fasta`, whose FASTA header is the sequence SHA256.
#    499,770 training chains share only 142,056 distinct sequences, and an MSA
#    is a pure function of (query sequence, reference database) -- so searching
#    a sequence once and attributing it to every chain that carries it is
#    exact, not an approximation. The pilot is what licenses this: it searched
#    a deliberately duplicated sequence pair independently and compared the
#    results.
#
#    The resulting store is content-addressed: `<sequence_sha256>.a3m`. The
#    chain -> sequence mapping already exists in `chain_to_sequence.parquet`,
#    so no per-chain copy is made. That is not only cheaper, it is necessary:
#    500k files would exceed this filesystem's inode quota, whereas 142k does
#    not.
#
# 3. It is resumable. A shard whose `_SHARD_SUCCESS` marker is present is
#    skipped, so a partially failed array is restarted by resubmitting it.
#
# Usage
# -----
#   sbatch --array=0-N%C run_msa_production.sbatch DB_ROOT QUERY_FASTA OUT_ROOT SHARD_COUNT
#
# Nothing here is snapshot-specific: every path is an argument.

set -euo pipefail

DB_ROOT="${1:?mmseqs database root required}"
QUERY_FASTA="${2:?unique-sequence FASTA required}"
OUT_ROOT="${3:?output root required}"
SHARD_COUNT="${4:?shard count required}"

DB_ROOT="$(readlink -f "$DB_ROOT")"
QUERY_FASTA="$(readlink -f "$QUERY_FASTA")"
mkdir -p "$OUT_ROOT"
OUT_ROOT="$(readlink -f "$OUT_ROOT")"

SHARD_ID="${SLURM_ARRAY_TASK_ID:?this script must be submitted as an array}"

if [[ ! "$SHARD_COUNT" =~ ^[0-9]+$ ]] || (( SHARD_COUNT < 1 )); then
    echo "ERROR: invalid shard count: $SHARD_COUNT" >&2
    exit 2
fi

if (( SHARD_ID >= SHARD_COUNT )); then
    echo "ERROR: shard $SHARD_ID outside 0..$((SHARD_COUNT - 1))" >&2
    exit 2
fi

DB_NAME="uniref30_2103_db"
DB="$DB_ROOT/db/$DB_NAME"
DB_MARKER="$DB_ROOT/db/_UNIREF30_DB_SUCCESS"

MMSEQS="${MMSEQS_BIN:-/users/sgsrehm1/COMP702_BeyondAF/envs/openfold_env/bin/mmseqs}"
THREADS="${SLURM_CPUS_PER_TASK:-16}"

SHARD_DIR="$OUT_ROOT/shards/$(printf '%05d' "$SHARD_ID")"
WORK="$SHARD_DIR/work"
MMSEQS_TMP="$WORK/tmp"
SHARD_FASTA="$SHARD_DIR/queries.fasta"
SHARD_SUCCESS="$SHARD_DIR/_SHARD_SUCCESS"

# The MSA store is shared by every shard: content-addressed, so two shards can
# never contend for the same name unless they hold the same sequence, which
# the sharding below makes impossible.
MSA_STORE="$OUT_ROOT/msas"

mkdir -p "$SHARD_DIR" "$MSA_STORE"

echo "============================================================"
echo "Stage 17 production MSA shard"
echo "Shard:      $SHARD_ID of $SHARD_COUNT"
echo "Job:        ${SLURM_JOB_ID:-manual} (array ${SLURM_ARRAY_JOB_ID:-none})"
echo "Node:       $(hostname)"
echo "Threads:    $THREADS"
echo "Database:   $DB"
echo "Query set:  $QUERY_FASTA"
echo "MSA store:  $MSA_STORE"
echo "Started:    $(date --iso-8601=seconds)"
echo "============================================================"

if [[ -f "$SHARD_SUCCESS" ]]; then
    echo "Shard already complete; nothing to do."
    exit 0
fi

if [[ ! -x "$MMSEQS" ]]; then
    echo "ERROR: mmseqs not executable at $MMSEQS" >&2
    exit 1
fi

MMSEQS_VERSION="$("$MMSEQS" version)"
echo "mmseqs version: $MMSEQS_VERSION"

if [[ ! -f "$DB_MARKER" ]]; then
    echo "ERROR: no $DB_MARKER." >&2
    echo "       Production refuses to run against an unverified database." >&2
    exit 1
fi

# Index-free operation, exactly as the pilot: prefer the .idx when one exists,
# otherwise the _seq/_aln sub-databases that tsv2exprofiledb produces.
if [[ -f "${DB}.idx" ]]; then
    DB_LOAD_MODE=2
    SUFFIX_SEQ=".idx"
    SUFFIX_ALN=".idx"
    INDEX_USED=".idx"
else
    DB_LOAD_MODE=0
    SUFFIX_SEQ="_seq"
    SUFFIX_ALN="_aln"
    INDEX_USED="_seq"
fi

echo "db_load_mode: $DB_LOAD_MODE (index: $INDEX_USED)"

# ---------------------------------------------------------------------------
# Shard the query set
# ---------------------------------------------------------------------------
# Contiguous slices by record index. Deterministic and disjoint, so a shard is
# reproducible on its own and no sequence is searched twice.

echo
echo "===== SHARD QUERY SET ====="

TOTAL_QUERIES="$(grep -c '^>' "$QUERY_FASTA")"

python3 - "$QUERY_FASTA" "$SHARD_FASTA" "$SHARD_ID" "$SHARD_COUNT" <<'PYTHON'
import sys

source, target, shard_id, shard_count = sys.argv[1:5]
shard_id, shard_count = int(shard_id), int(shard_count)

written = 0
index = -1
emit = False

with open(source) as handle, open(target, "w") as out:
    for line in handle:
        if line.startswith(">"):
            index += 1
            emit = (index % shard_count) == shard_id
            if emit:
                written += 1
        if emit:
            out.write(line)

print(f"shard sequences: {written}")
PYTHON

SHARD_QUERIES="$(grep -c '^>' "$SHARD_FASTA" || true)"

echo "total sequences: $TOTAL_QUERIES"
echo "this shard:      $SHARD_QUERIES"

if (( SHARD_QUERIES == 0 )); then
    echo "Empty shard; nothing to search."
    date --iso-8601=seconds > "$SHARD_SUCCESS"
    exit 0
fi

rm -rf "$WORK"
mkdir -p "$WORK" "$MMSEQS_TMP"

# ---------------------------------------------------------------------------
# The protocol -- identical to the pilot
# ---------------------------------------------------------------------------

echo
echo "===== QUERY DATABASE ====="

"$MMSEQS" createdb "$SHARD_FASTA" "$WORK/qdb"

QDB_ENTRIES="$(wc -l < "$WORK/qdb.index")"

echo "query database entries: $QDB_ENTRIES"

if [[ "$QDB_ENTRIES" != "$SHARD_QUERIES" ]]; then
    echo "ERROR: query database has $QDB_ENTRIES entries, expected $SHARD_QUERIES" >&2
    exit 1
fi

echo
echo "===== SEARCH (3 profile iterations) ====="
echo "started: $(date --iso-8601=seconds)"

"$MMSEQS" search \
    "$WORK/qdb" \
    "$DB" \
    "$WORK/res" \
    "$MMSEQS_TMP" \
    --num-iterations 3 \
    --db-load-mode "$DB_LOAD_MODE" \
    -a \
    -e 0.1 \
    --max-seqs 10000 \
    --threads "$THREADS"

echo "finished: $(date --iso-8601=seconds)"

echo
echo "===== PROFILE HANDOFF ====="

"$MMSEQS" mvdb "$MMSEQS_TMP/latest/profile_1" "$WORK/prof_res"
"$MMSEQS" lndb "$WORK/qdb_h" "$WORK/prof_res_h"

echo
echo "===== EXPAND / ALIGN / FILTER ====="

"$MMSEQS" expandaln \
    "$WORK/qdb" \
    "${DB}${SUFFIX_SEQ}" \
    "$WORK/res" \
    "${DB}${SUFFIX_ALN}" \
    "$WORK/res_exp" \
    --db-load-mode "$DB_LOAD_MODE" \
    --expansion-mode 0 \
    -e inf \
    --expand-filter-clusters 0 \
    --max-seq-id 0.95 \
    --threads "$THREADS"

"$MMSEQS" align \
    "$WORK/prof_res" \
    "${DB}${SUFFIX_SEQ}" \
    "$WORK/res_exp" \
    "$WORK/res_exp_realign" \
    --db-load-mode "$DB_LOAD_MODE" \
    -e 10 \
    --max-accept 1000000 \
    --alt-ali 10 \
    -a \
    --threads "$THREADS"

"$MMSEQS" filterresult \
    "$WORK/qdb" \
    "${DB}${SUFFIX_SEQ}" \
    "$WORK/res_exp_realign" \
    "$WORK/res_exp_realign_filter" \
    --db-load-mode "$DB_LOAD_MODE" \
    --qid 0 \
    --qsc 0.8 \
    --diff 0 \
    --max-seq-id 1.0 \
    --filter-min-enable 100 \
    --threads "$THREADS"

echo
echo "===== MSA CONSTRUCTION ====="

"$MMSEQS" result2msa \
    "$WORK/qdb" \
    "${DB}${SUFFIX_SEQ}" \
    "$WORK/res_exp_realign_filter" \
    "$WORK/uniref.a3m" \
    --msa-format-mode 6 \
    --db-load-mode "$DB_LOAD_MODE" \
    --filter-msa 1 \
    --filter-min-enable 1000 \
    --diff 3000 \
    --qid "0.0,0.2,0.4,0.6,0.8,1.0" \
    --qsc 0 \
    --max-seq-id 0.95 \
    --threads "$THREADS"

# ---------------------------------------------------------------------------
# Unpack, named by sequence identity
# ---------------------------------------------------------------------------
# `--unpack-name-mode 1` resolves names through the `.lookup` of the database
# being unpacked -- the RESULT database, not the query database. `result2msa`
# emits none, so mode 1 falls back to numeric DB keys and the identity is lost.
# This is the defect the pilot caught; see commit "Fix MSA output naming".

RESULT_KEYS="$(cut -f1 "$WORK/uniref.a3m.index" | sort -n | tr '\n' ' ')"
QUERY_KEYS="$(cut -f1 "$WORK/qdb.lookup" | sort -n | tr '\n' ' ')"

if [[ "$RESULT_KEYS" != "$QUERY_KEYS" ]]; then
    echo "ERROR: result database keys do not match the query database keys;" >&2
    echo "       qdb.lookup would mislabel the unpacked MSAs." >&2
    exit 1
fi

LOOKUP_TOTAL="$(wc -l < "$WORK/qdb.lookup")"
LOOKUP_UNIQUE="$(cut -f2 "$WORK/qdb.lookup" | sort -u | wc -l)"

if [[ "$LOOKUP_TOTAL" != "$LOOKUP_UNIQUE" ]]; then
    echo "ERROR: duplicate accessions in qdb.lookup; outputs would collide" >&2
    cut -f2 "$WORK/qdb.lookup" | sort | uniq -d | sed 's/^/  duplicate: /' >&2
    exit 1
fi

cp "$WORK/qdb.lookup" "$WORK/uniref.a3m.lookup"

# Unpack into the shard first, then publish. A killed job must never leave a
# half-written .a3m in the shared store where a later stage would read it.
UNPACK="$SHARD_DIR/unpacked"
rm -rf "$UNPACK"
mkdir -p "$UNPACK"

"$MMSEQS" unpackdb \
    "$WORK/uniref.a3m" \
    "$UNPACK" \
    --unpack-name-mode 1 \
    --unpack-suffix .a3m

PRODUCED="$(find "$UNPACK" -maxdepth 1 -name '*.a3m' | wc -l)"

echo "a3m produced: $PRODUCED"
echo "expected:     $SHARD_QUERIES"

if (( PRODUCED != SHARD_QUERIES )); then
    echo "ERROR: produced $PRODUCED MSAs for $SHARD_QUERIES queries" >&2
    exit 1
fi

# Every name must be a sequence SHA256 that this shard actually asked for.
python3 - "$SHARD_FASTA" "$UNPACK" <<'PYTHON'
import os
import sys

fasta, unpacked = sys.argv[1:3]

expected = {
    line[1:].split()[0]
    for line in open(fasta)
    if line.startswith(">")
}

produced = {
    name[:-4]
    for name in os.listdir(unpacked)
    if name.endswith(".a3m")
}

missing = expected - produced
orphan = produced - expected

if missing or orphan:
    for name in sorted(missing)[:10]:
        print(f"  missing MSA for query: {name}", file=sys.stderr)
    for name in sorted(orphan)[:10]:
        print(f"  orphan output not traceable to a query: {name}", file=sys.stderr)
    sys.exit(1)

print(f"identity check: {len(produced)} MSAs, all traceable to a query")
PYTHON

echo
echo "===== PUBLISH ====="

mv -f "$UNPACK"/*.a3m "$MSA_STORE"/
rmdir "$UNPACK"

# Intermediates are large and reproducible; the MSAs are what matters.
rm -rf "$WORK"

{
    printf 'stage\topenfold_msa_production_shard\n'
    printf 'shard_id\t%s\n' "$SHARD_ID"
    printf 'shard_count\t%s\n' "$SHARD_COUNT"
    printf 'shard_queries\t%s\n' "$SHARD_QUERIES"
    printf 'a3m_produced\t%s\n' "$PRODUCED"
    printf 'mmseqs_version\t%s\n' "$MMSEQS_VERSION"
    printf 'mmseqs_binary\t%s\n' "$MMSEQS"
    printf 'reference_database\t%s\n' "$DB"
    printf 'db_load_mode\t%s\n' "$DB_LOAD_MODE"
    printf 'index_used\t%s\n' "$INDEX_USED"
    printf 'query_fasta\t%s\n' "$QUERY_FASTA"
    printf 'threads\t%s\n' "$THREADS"
    printf 'slurm_job_id\t%s\n' "${SLURM_JOB_ID:-manual}"
    printf 'slurm_array_job_id\t%s\n' "${SLURM_ARRAY_JOB_ID:-none}"
    printf 'node\t%s\n' "$(hostname)"
    printf 'completed_at\t%s\n' "$(date --iso-8601=seconds)"
} > "$SHARD_DIR/shard_provenance.tsv"

date --iso-8601=seconds > "$SHARD_SUCCESS"

echo
echo "============================================================"
echo "SHARD $SHARD_ID: PASS  ($PRODUCED MSAs)"
echo "Completed: $(date --iso-8601=seconds)"
echo "============================================================"
