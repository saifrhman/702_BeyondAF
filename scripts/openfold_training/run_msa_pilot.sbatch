#!/bin/bash -l
#SBATCH --job-name=of_msa_pilot
#SBATCH --partition=nodes
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=400G
#SBATCH --time=12:00:00

# Stage 17 pilot: prove the MSA machinery on a small representative chain set
# before ~500k production searches are committed.
#
# Protocol
# --------
# `colabfold_search` is not installed in this environment, so the ColabFold
# expandable-profile protocol is issued directly against `mmseqs`.  The command
# sequence and parameters below mirror ColabFold's own uniref path so the MSAs
# this pipeline produces are the ones OpenFold expects:
#
#   createdb -> search (3 iterations) -> expandaln -> align
#            -> filterresult -> result2msa -> unpackdb
#
# Index-free operation is deliberate, not a shortcut.  ColabFold selects
# `${db}.idx` when an index exists and falls back to the `${db}_seq` /
# `${db}_aln` sub-databases when it does not.  `tsv2exprofiledb` produces those
# sub-databases, so the pilot is scientifically equivalent and merely slower.
# Whether to spend the quota on `createindex` is a throughput decision for
# production, and it is made after this pilot reports real timings.
#
# One search per chain
# --------------------
# Every pilot chain enters the query database as its own entry and leaves as
# its own `.a3m`.  Two of the chains carry a byte-identical sequence.  They are
# searched independently anyway.  The validation below refuses to pass if the
# outputs were shared, linked, or fanned out from a single search.
#
# Nothing here is snapshot-specific: every path is an argument.

set -euo pipefail

DB_ROOT="${1:?mmseqs database root required}"
PILOT_DIR="${2:?pilot directory required}"

DB_ROOT="$(readlink -f "$DB_ROOT")"
PILOT_DIR="$(readlink -f "$PILOT_DIR")"

DB_NAME="uniref30_2103_db"
DB="$DB_ROOT/db/$DB_NAME"
DB_MARKER="$DB_ROOT/db/_UNIREF30_DB_SUCCESS"

PILOT_FASTA="$PILOT_DIR/pilot_queries.fasta"
PILOT_MANIFEST="$PILOT_DIR/pilot_selection_manifest.json"

WORK="$PILOT_DIR/work"
OUT="$PILOT_DIR/msas"
MMSEQS_TMP="$WORK/tmp"

MMSEQS="${MMSEQS_BIN:-/users/sgsrehm1/COMP702_BeyondAF/envs/openfold_env/bin/mmseqs}"
THREADS="${SLURM_CPUS_PER_TASK:-16}"

echo "============================================================"
echo "STAGE 17 MSA PILOT"
echo "Job:      ${SLURM_JOB_ID:-none}"
echo "Node:     $(hostname)"
echo "CPUs:     $THREADS"
echo "Started:  $(date --iso-8601=seconds)"
echo "============================================================"

echo
echo "===== PREREQUISITES ====="

if [[ ! -x "$MMSEQS" ]]; then
    echo "ERROR: mmseqs not executable at $MMSEQS"
    exit 1
fi

MMSEQS_VERSION="$("$MMSEQS" version)"
echo "mmseqs version: $MMSEQS_VERSION"

if [[ ! -f "$DB_MARKER" ]]; then
    echo "ERROR: UniRef30 database is not marked ready."
    echo "       expected: $DB_MARKER"
    echo "       The pilot refuses to run against an unverified database."
    exit 1
fi

echo "database marker:"
sed 's/^/  /' "$DB_MARKER"

# tsv2exprofiledb produces these; an expandable-profile search cannot run
# without them, and their absence means the preparation job did not finish
# whatever its exit status claimed.
for suffix in "" "_seq" "_aln" "_h"; do
    found=0
    for probe in "${DB}${suffix}" "${DB}${suffix}.dbtype" "${DB}${suffix}.index"; do
        if [[ -e "$probe" ]]; then
            found=1
            break
        fi
    done

    if (( found == 0 )); then
        echo "ERROR: required sub-database missing: ${DB_NAME}${suffix}"
        exit 1
    fi

    echo "present: ${DB_NAME}${suffix}"
done

if [[ -f "${DB}.idx" || -f "${DB}.idx.dbtype" ]]; then
    echo "index:          present -> using .idx (db-load-mode 2)"
    DB_LOAD_MODE=2
    SUFFIX_SEQ=".idx"
    SUFFIX_ALN=".idx"
else
    echo "index:          absent  -> using _seq/_aln (db-load-mode 0)"
    DB_LOAD_MODE=0
    SUFFIX_SEQ="_seq"
    SUFFIX_ALN="_aln"
fi

if [[ ! -f "$PILOT_FASTA" ]]; then
    echo "ERROR: pilot FASTA not found: $PILOT_FASTA"
    exit 1
fi

EXPECTED_QUERIES="$(grep -c '^>' "$PILOT_FASTA")"

echo "pilot FASTA:    $PILOT_FASTA"
echo "pilot queries:  $EXPECTED_QUERIES"
echo "pilot manifest: $PILOT_MANIFEST"

echo
echo "===== WORKSPACE ====="

# A previous partial pilot must not be mistaken for a fresh result.
rm -rf "$WORK" "$OUT"
mkdir -p "$WORK" "$OUT" "$MMSEQS_TMP"

echo "work: $WORK"
echo "out:  $OUT"

echo
echo "===== QUERY DATABASE ====="

"$MMSEQS" createdb "$PILOT_FASTA" "$WORK/qdb"

QDB_ENTRIES="$(wc -l < "$WORK/qdb.index")"

echo "query database entries: $QDB_ENTRIES"

if [[ "$QDB_ENTRIES" != "$EXPECTED_QUERIES" ]]; then
    echo "ERROR: query database has $QDB_ENTRIES entries, expected $EXPECTED_QUERIES"
    echo "       Every chain must be its own independent query."
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

# The third iteration's profile becomes the query for realignment, exactly as
# ColabFold does; the header database is linked so the profile keeps identities.
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

# --unpack-name-mode 1 names each output by its accession, resolved through
# qdb.lookup -- i.e. the FASTA header, which is the openfold_chain_id.  Mode 0
# would name files by numeric DB key (0.a3m, 1.a3m, ...) and lose the chain
# identity entirely.  ColabFold uses mode 0 and renames afterwards from its own
# ordering; naming directly by accession removes that reordering step and the
# chance of getting it wrong.
#
# Mode 1 overwrites on accession collision, so uniqueness is asserted first --
# silently losing a chain's MSA to a namesake is precisely the failure this
# stage must never have.
LOOKUP_TOTAL="$(wc -l < "$WORK/qdb.lookup")"
LOOKUP_UNIQUE="$(cut -f2 "$WORK/qdb.lookup" | sort -u | wc -l)"

echo "qdb.lookup accessions: $LOOKUP_TOTAL total, $LOOKUP_UNIQUE unique"

if [[ "$LOOKUP_TOTAL" != "$LOOKUP_UNIQUE" ]]; then
    echo "ERROR: duplicate accessions in qdb.lookup; outputs would collide"
    cut -f2 "$WORK/qdb.lookup" | sort | uniq -d | sed 's/^/  duplicate: /'
    exit 1
fi

if [[ "$LOOKUP_TOTAL" != "$EXPECTED_QUERIES" ]]; then
    echo "ERROR: qdb.lookup has $LOOKUP_TOTAL rows, expected $EXPECTED_QUERIES"
    exit 1
fi

# `--unpack-name-mode 1` resolves names through the `.lookup` of the database
# being unpacked -- the RESULT database, not the query database. `result2msa`
# does not emit one, so mode 1 silently fell back to DB keys and produced
# 0.a3m .. 8.a3m with the chain identity lost. The validation below caught it.
#
# `result2msa` preserves the query keys, which the assertion below re-checks,
# so `qdb.lookup` *is* the correct lookup for the result database. Installing
# it is what makes mode 1 resolve to the openfold_chain_id.
RESULT_KEYS="$(cut -f1 "$WORK/uniref.a3m.index" | sort -n | tr '\n' ' ')"
QUERY_KEYS="$(cut -f1 "$WORK/qdb.lookup" | sort -n | tr '\n' ' ')"

if [[ "$RESULT_KEYS" != "$QUERY_KEYS" ]]; then
    echo "ERROR: result database keys do not match the query database keys;"
    echo "       qdb.lookup would mislabel the unpacked MSAs."
    echo "  result: $RESULT_KEYS"
    echo "  query:  $QUERY_KEYS"
    exit 1
fi

cp "$WORK/qdb.lookup" "$WORK/uniref.a3m.lookup"

echo "installed result-db lookup: $(wc -l < "$WORK/uniref.a3m.lookup") accessions"

"$MMSEQS" unpackdb \
    "$WORK/uniref.a3m" \
    "$OUT" \
    --unpack-name-mode 1 \
    --unpack-suffix .a3m

echo "unpacked: $(date --iso-8601=seconds)"

echo
echo "===== OUTPUTS ====="

ls -l "$OUT" | sed 's/^/  /'

PRODUCED="$(find "$OUT" -maxdepth 1 -name '*.a3m' | wc -l)"

echo
echo "a3m files produced: $PRODUCED"
echo "queries expected:   $EXPECTED_QUERIES"

echo
echo "===== VALIDATION ====="

# The pilot is only meaningful if the outputs are per-chain, correctly
# identified, and independently produced.  Each check below can fail the job.
# The validation exit code must be captured, not acted on immediately: the
# provenance record below is written whether the pilot passes or fails, because
# a failed pilot is evidence too.
set +e

python3 - "$PILOT_FASTA" "$OUT" "$PILOT_DIR/pilot_validation.json" <<'PYTHON'
import hashlib
import json
import os
import sys

fasta_path, out_dir, report_path = sys.argv[1:4]

# --- the exact queries that were submitted ---------------------------------
queries = {}
name = None
chunks = []

with open(fasta_path, encoding="utf-8") as handle:
    for line in handle:
        line = line.rstrip("\n")
        if line.startswith(">"):
            if name is not None:
                queries[name] = "".join(chunks)
            name = line[1:].strip()
            chunks = []
        else:
            chunks.append(line.strip())

if name is not None:
    queries[name] = "".join(chunks)

failures = []
records = {}

for chain_id, sequence in queries.items():
    path = os.path.join(out_dir, f"{chain_id}.a3m")

    if not os.path.isfile(path):
        failures.append(f"{chain_id}: no .a3m produced")
        continue

    stat = os.stat(path)

    with open(path, encoding="utf-8") as handle:
        content = handle.read()

    lines = content.splitlines()

    # The first record of a ColabFold a3m is the query itself.
    first_header = None
    first_seq_parts = []

    for index, line in enumerate(lines):
        if line.startswith(">"):
            if first_header is None:
                first_header = line[1:].strip()
                continue
            break
        if first_header is not None:
            first_seq_parts.append(line.strip())

    first_sequence = "".join(first_seq_parts)

    # a3m lower-case letters denote insertions relative to the query; the
    # query row itself must reproduce the submitted sequence exactly.
    ungapped = "".join(c for c in first_sequence if c.isupper() or c == "-")
    ungapped = ungapped.replace("-", "")

    if ungapped != sequence:
        failures.append(
            f"{chain_id}: a3m query row does not match the submitted sequence "
            f"(a3m {len(ungapped)} aa vs query {len(sequence)} aa)"
        )

    depth = sum(1 for line in lines if line.startswith(">"))

    records[chain_id] = {
        "a3m_path": path,
        "bytes": stat.st_size,
        "inode": stat.st_ino,
        "hard_links": stat.st_nlink,
        "is_symlink": os.path.islink(path),
        "msa_depth": depth,
        "query_length": len(sequence),
        "query_row_matches": ungapped == sequence,
        "content_sha256": hashlib.sha256(content.encode()).hexdigest(),
    }

    if os.path.islink(path):
        failures.append(f"{chain_id}: output is a symlink, not an independent MSA")

    if stat.st_nlink != 1:
        failures.append(
            f"{chain_id}: output has {stat.st_nlink} hard links; "
            "outputs must not be shared between chains"
        )

# --- the sequence-sharing pair --------------------------------------------
# Identical inputs legitimately yield identical MSA *content* under a
# deterministic search.  What must NOT happen is one search being reused for
# both chains.  Independence is therefore asserted on identity -- separate
# files, separate inodes, no linking -- never on content differing.
by_sequence = {}

for chain_id, sequence in queries.items():
    by_sequence.setdefault(sequence, []).append(chain_id)

shared_groups = {
    hashlib.sha256(seq.encode()).hexdigest(): ids
    for seq, ids in by_sequence.items()
    if len(ids) > 1
}

sharing_report = []

for seq_hash, ids in shared_groups.items():
    present = [i for i in ids if i in records]

    inodes = {records[i]["inode"] for i in present}

    independent = len(inodes) == len(present) and len(present) == len(ids)

    if not independent:
        failures.append(
            f"sequence {seq_hash[:16]}: chains {ids} do not have independent "
            "output identities"
        )

    sharing_report.append(
        {
            "sequence_sha256": seq_hash,
            "chains": ids,
            "distinct_inodes": len(inodes),
            "independent_outputs": independent,
            "content_identical": len(
                {records[i]["content_sha256"] for i in present}
            ) == 1,
        }
    )

report = {
    "stage": "openfold_msa_pilot_validation",
    "queries_submitted": len(queries),
    "a3m_produced": len(records),
    "missing": sorted(set(queries) - set(records)),
    "orphan_outputs": sorted(
        f[:-4]
        for f in os.listdir(out_dir)
        if f.endswith(".a3m") and f[:-4] not in queries
    ),
    "sequence_sharing_groups": sharing_report,
    "chains": records,
    "failures": failures,
    "status": "PASS" if not failures else "FAIL",
}

if report["orphan_outputs"]:
    report["failures"].append(
        f"orphan outputs not traceable to a query: {report['orphan_outputs']}"
    )
    report["status"] = "FAIL"

with open(report_path, "w", encoding="utf-8") as handle:
    json.dump(report, handle, indent=2, sort_keys=True)
    handle.write("\n")

print(f"queries submitted: {report['queries_submitted']}")
print(f"a3m produced:      {report['a3m_produced']}")
print(f"missing:           {report['missing']}")
print(f"orphan outputs:    {report['orphan_outputs']}")
print()

for chain_id in sorted(records):
    record = records[chain_id]
    print(
        f"  {chain_id:12s} depth={record['msa_depth']:7d} "
        f"len={record['query_length']:5d} "
        f"query_row_ok={record['query_row_matches']} "
        f"bytes={record['bytes']}"
    )

print()

for group in sharing_report:
    print(
        f"  shared sequence {group['sequence_sha256'][:16]}: "
        f"chains={group['chains']} "
        f"independent_outputs={group['independent_outputs']} "
        f"content_identical={group['content_identical']}"
    )

print()

if failures:
    print("FAILURES:")
    for failure in failures:
        print(f"  - {failure}")
    sys.exit(1)

print("All pilot validations passed.")
PYTHON

VALIDATION_STATUS=$?

set -e

echo
echo "===== TIMING / RESOURCE ====="

echo "working directory sizes:"
du -sh "$WORK" "$OUT" 2>/dev/null | sed 's/^/  /'

echo
echo "===== PROVENANCE ====="

{
    printf 'stage\topenfold_msa_pilot\n'
    printf 'mmseqs_version\t%s\n' "$MMSEQS_VERSION"
    printf 'mmseqs_binary\t%s\n' "$MMSEQS"
    printf 'reference_database\t%s\n' "$DB"
    printf 'db_load_mode\t%s\n' "$DB_LOAD_MODE"
    printf 'index_used\t%s\n' "$SUFFIX_SEQ"
    printf 'pilot_fasta\t%s\n' "$PILOT_FASTA"
    printf 'pilot_fasta_sha256\t%s\n' "$(sha256sum "$PILOT_FASTA" | awk '{print $1}')"
    printf 'queries\t%s\n' "$EXPECTED_QUERIES"
    printf 'a3m_produced\t%s\n' "$PRODUCED"
    printf 'threads\t%s\n' "$THREADS"
    printf 'slurm_job_id\t%s\n' "${SLURM_JOB_ID:-none}"
    printf 'node\t%s\n' "$(hostname)"
    printf 'completed_at\t%s\n' "$(date --iso-8601=seconds)"
} > "$PILOT_DIR/pilot_provenance.tsv"

cat "$PILOT_DIR/pilot_provenance.tsv"

if (( VALIDATION_STATUS != 0 )); then
    echo
    echo "STATUS: FAIL"
    exit 1
fi

{
    printf 'STATUS=PASS\n'
    printf 'QUERIES=%s\n' "$EXPECTED_QUERIES"
    printf 'A3M_PRODUCED=%s\n' "$PRODUCED"
    printf 'MMSEQS_VERSION=%s\n' "$MMSEQS_VERSION"
    printf 'REFERENCE_DB=%s\n' "$DB"
    printf 'SLURM_JOB_ID=%s\n' "${SLURM_JOB_ID:-none}"
    printf 'COMPLETED_AT=%s\n' "$(date --iso-8601=seconds)"
} > "$PILOT_DIR/_PILOT_SUCCESS"

echo
echo "============================================================"
echo "STATUS: PASS"
echo "Stage 17 pilot succeeded on $EXPECTED_QUERIES representative chains."
echo "No production MSAs generated.  No run-specific data modified."
echo "Completed: $(date --iso-8601=seconds)"
echo "============================================================"
