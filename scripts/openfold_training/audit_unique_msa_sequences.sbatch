#!/bin/bash -l
#SBATCH --job-name=of_msa_unique
#SBATCH --partition=nodes
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00

set -euo pipefail

TRAIN_ROOT="${1:?Usage: sbatch audit_unique_msa_sequences.sbatch <training_root> <python_executable>}"
PY="${2:?Usage: sbatch audit_unique_msa_sequences.sbatch <training_root> <python_executable>}"

TRAIN_ROOT="$(readlink -f "$TRAIN_ROOT")"
PY="$(readlink -f "$PY")"

MSA_DIR="$TRAIN_ROOT/msa"
QUERY_PARQUET="$MSA_DIR/msa_queries.parquet"
QUERY_MANIFEST="$MSA_DIR/msa_query_manifest.json"
OUT_DIR="$MSA_DIR/sequence_dedup"

mkdir -p "$OUT_DIR"

export PYTHONUNBUFFERED=1
export PYTHONDONTWRITEBYTECODE=1

echo "============================================================"
echo "OPENFOLD MSA EXACT-SEQUENCE AUDIT"
echo "Job:          ${SLURM_JOB_ID}"
echo "Node:         $(hostname)"
echo "Training root: $TRAIN_ROOT"
echo "Started:      $(date --iso-8601=seconds)"
echo "============================================================"

"$PY" - "$QUERY_PARQUET" "$QUERY_MANIFEST" "$OUT_DIR" <<'PY'
from pathlib import Path
import hashlib
import json
import sys

import pyarrow as pa
import pyarrow.parquet as pq


query_path = Path(sys.argv[1])
manifest_path = Path(sys.argv[2])
out_dir = Path(sys.argv[3])

if not query_path.is_file():
    raise RuntimeError(f"Missing query parquet: {query_path}")

if not manifest_path.is_file():
    raise RuntimeError(f"Missing query manifest: {manifest_path}")


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


manifest = json.loads(manifest_path.read_text())

expected_query_count = int(manifest["query_count"])
expected_query_sha = manifest["manifest_sha256"]

actual_query_sha = sha256_file(query_path)

if actual_query_sha != expected_query_sha:
    raise RuntimeError(
        "msa_queries.parquet SHA256 does not match its manifest"
    )

table = pq.read_table(
    query_path,
    columns=[
        "query_index",
        "openfold_chain_id",
        "pdb_id",
        "sequence_length",
        "sequence",
    ],
)

rows = table.to_pylist()

if len(rows) != expected_query_count:
    raise RuntimeError(
        f"Manifest says {expected_query_count} queries, "
        f"but parquet contains {len(rows)}"
    )

# sequence -> record
unique = {}

# SHA256 -> sequence, used to explicitly guard against a hash collision
hash_to_sequence = {}

mapping_rows = []

for row in rows:
    seq = row["sequence"]

    if len(seq) != int(row["sequence_length"]):
        raise RuntimeError(
            f"Length mismatch for {row['openfold_chain_id']}"
        )

    digest = hashlib.sha256(seq.encode("ascii")).hexdigest()

    previous = hash_to_sequence.get(digest)
    if previous is not None and previous != seq:
        raise RuntimeError(
            f"SHA256 collision detected for digest {digest}"
        )

    hash_to_sequence[digest] = seq

    rec = unique.get(seq)

    if rec is None:
        rec = {
            "sequence_sha256": digest,
            "representative_openfold_chain_id":
                row["openfold_chain_id"],
            "representative_query_index":
                int(row["query_index"]),
            "sequence_length":
                int(row["sequence_length"]),
            "multiplicity": 0,
            "sequence": seq,
        }
        unique[seq] = rec

    rec["multiplicity"] += 1

    mapping_rows.append({
        "query_index":
            int(row["query_index"]),
        "openfold_chain_id":
            row["openfold_chain_id"],
        "pdb_id":
            row["pdb_id"],
        "sequence_sha256":
            digest,
        "sequence_length":
            int(row["sequence_length"]),
    })

unique_rows = sorted(
    unique.values(),
    key=lambda r: r["representative_query_index"],
)

mapping_rows = sorted(
    mapping_rows,
    key=lambda r: r["query_index"],
)

unique_count = len(unique_rows)
query_count = len(rows)
reused_chain_count = query_count - unique_count

singleton_sequences = sum(
    r["multiplicity"] == 1
    for r in unique_rows
)

shared_sequences = sum(
    r["multiplicity"] > 1
    for r in unique_rows
)

max_multiplicity = max(
    r["multiplicity"]
    for r in unique_rows
)

total_unique_residues = sum(
    r["sequence_length"]
    for r in unique_rows
)

unique_schema = pa.schema([
    pa.field("sequence_sha256", pa.string(), nullable=False),
    pa.field(
        "representative_openfold_chain_id",
        pa.string(),
        nullable=False,
    ),
    pa.field(
        "representative_query_index",
        pa.int64(),
        nullable=False,
    ),
    pa.field("sequence_length", pa.int32(), nullable=False),
    pa.field("multiplicity", pa.int32(), nullable=False),
    pa.field("sequence", pa.string(), nullable=False),
])

mapping_schema = pa.schema([
    pa.field("query_index", pa.int64(), nullable=False),
    pa.field("openfold_chain_id", pa.string(), nullable=False),
    pa.field("pdb_id", pa.string(), nullable=False),
    pa.field("sequence_sha256", pa.string(), nullable=False),
    pa.field("sequence_length", pa.int32(), nullable=False),
])

unique_path = out_dir / "unique_sequences.parquet"
mapping_path = out_dir / "chain_to_sequence.parquet"
fasta_path = out_dir / "unique_sequences.fasta"

pq.write_table(
    pa.Table.from_pylist(
        unique_rows,
        schema=unique_schema,
    ),
    unique_path,
    compression="zstd",
)

pq.write_table(
    pa.Table.from_pylist(
        mapping_rows,
        schema=mapping_schema,
    ),
    mapping_path,
    compression="zstd",
)

with fasta_path.open("w") as f:
    for row in unique_rows:
        f.write(
            f">{row['sequence_sha256']}\n"
            f"{row['sequence']}\n"
        )

stats = {
    "query_count": query_count,
    "unique_sequence_count": unique_count,
    "chains_reusing_an_identical_sequence":
        reused_chain_count,
    "singleton_sequence_count":
        singleton_sequences,
    "shared_sequence_count":
        shared_sequences,
    "maximum_sequence_multiplicity":
        max_multiplicity,
    "total_query_residues":
        int(manifest["total_residues"]),
    "total_unique_sequence_residues":
        total_unique_residues,
    "searches_avoided_by_exact_sequence_reuse":
        reused_chain_count,
    "unique_sequence_fraction":
        unique_count / query_count,
    "query_manifest_sha256":
        expected_query_sha,
    "unique_sequences_parquet":
        str(unique_path),
    "unique_sequences_parquet_sha256":
        sha256_file(unique_path),
    "chain_mapping_parquet":
        str(mapping_path),
    "chain_mapping_parquet_sha256":
        sha256_file(mapping_path),
    "unique_fasta":
        str(fasta_path),
    "unique_fasta_sha256":
        sha256_file(fasta_path),
}

stats_path = out_dir / "sequence_dedup_manifest.json"

stats_path.write_text(
    json.dumps(
        stats,
        indent=2,
        sort_keys=True,
    ) + "\n"
)

print()
print("===== EXACT-SEQUENCE RESULTS =====")
print("training chains:", query_count)
print("unique retained sequences:", unique_count)
print(
    "chains reusing an identical sequence:",
    reused_chain_count,
)
print("singleton sequences:", singleton_sequences)
print("shared sequences:", shared_sequences)
print("maximum multiplicity:", max_multiplicity)
print(
    "searches required at most:",
    unique_count,
)
print(
    "MSA searches avoided by exact reuse:",
    reused_chain_count,
)
print(
    "unique fraction:",
    f"{100.0 * unique_count / query_count:.2f}%",
)
print()
print("unique parquet:", unique_path)
print("mapping parquet:", mapping_path)
print("unique FASTA:", fasta_path)
print("manifest:", stats_path)
print()
print("STATUS: PASS")
PY

echo
echo "Completed: $(date --iso-8601=seconds)"
