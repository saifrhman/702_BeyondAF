#!/bin/bash -l
#SBATCH --job-name=of26_msa_queries
#SBATCH --partition=nodes
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=01:00:00
#SBATCH --output=/users/%u/fastscratch/COMP702_pdbclean_pipeline/outputs/openfold_training/20260101/pdbclean-dedup-v1/msa/build_queries_%j.out
#SBATCH --error=/users/%u/fastscratch/COMP702_pdbclean_pipeline/outputs/openfold_training/20260101/pdbclean-dedup-v1/msa/build_queries_%j.err

set -euo pipefail

ROOT="$HOME/fastscratch/COMP702_pdbclean_pipeline"
TRAIN="$ROOT/outputs/openfold_training/20260101/pdbclean-dedup-v1"
MSA="$TRAIN/msa"

PY="$HOME/fastscratch/envs/bri_env_1.2.2/bin/python"

export PYTHONUNBUFFERED=1
export PYTHONDONTWRITEBYTECODE=1
export XDG_CACHE_HOME="$HOME/fastscratch/.cache"

cd "$ROOT"
mkdir -p "$MSA"

echo "========================================"
echo "OpenFold 20260101 MSA query construction"
echo "Job:     ${SLURM_JOB_ID}"
echo "Node:    $(hostname)"
echo "Started: $(date --iso-8601=seconds)"
echo "========================================"

"$PY" - "$TRAIN" "$MSA" <<'PY'
from pathlib import Path
import hashlib
import json
import sys

import pyarrow as pa
import pyarrow.parquet as pq


EXPECTED = 499770

train = Path(sys.argv[1])
msa = Path(sys.argv[2])

src = train / "openfold_training_chains.parquet"

table = pq.read_table(
    src,
    columns=[
        "pdb_id",
        "openfold_chain_id",
        "retained_residue_count",
        "retained_sequence",
    ],
)

rows = table.to_pylist()

if len(rows) != EXPECTED:
    raise RuntimeError(
        f"Expected {EXPECTED} training chains, got {len(rows)}"
    )

seen = set()
bad_lengths = []
empty = []
alphabet = set()

out_rows = []

for i, row in enumerate(rows):
    key = row["openfold_chain_id"]
    seq = row["retained_sequence"]
    expected_len = int(row["retained_residue_count"])

    if key in seen:
        raise RuntimeError(
            f"Duplicate OpenFold key: {key}"
        )

    seen.add(key)

    if not seq:
        empty.append(key)
        continue

    if len(seq) != expected_len:
        bad_lengths.append(
            (key, len(seq), expected_len)
        )

    alphabet.update(seq)

    out_rows.append({
        "query_index": i,
        "openfold_chain_id": key,
        "pdb_id": row["pdb_id"],
        "sequence_length": len(seq),
        "sequence": seq,
    })

print("===== MSA QUERY AUDIT =====")
print("training chains:", len(rows))
print("unique OpenFold keys:", len(seen))
print("empty sequences:", len(empty))
print("length mismatches:", len(bad_lengths))
print("observed sequence alphabet:", "".join(sorted(alphabet)))

if empty:
    print("first empty:", empty[:20])

if bad_lengths:
    print("first length mismatches:", bad_lengths[:20])

if empty or bad_lengths:
    raise RuntimeError(
        "MSA query population validation failed"
    )

schema = pa.schema([
    pa.field("query_index", pa.int64(), nullable=False),
    pa.field("openfold_chain_id", pa.string(), nullable=False),
    pa.field("pdb_id", pa.string(), nullable=False),
    pa.field("sequence_length", pa.int32(), nullable=False),
    pa.field("sequence", pa.string(), nullable=False),
])

query_table = pa.Table.from_pylist(
    out_rows,
    schema=schema,
)

manifest_path = msa / "msa_queries.parquet"

pq.write_table(
    query_table,
    manifest_path,
    compression="zstd",
)

fasta_path = msa / "msa_queries.fasta"

with fasta_path.open("w") as f:
    for row in out_rows:
        f.write(
            f'>{row["openfold_chain_id"]}\n'
            f'{row["sequence"]}\n'
        )


def sha256(path):
    h = hashlib.sha256()

    with path.open("rb") as f:
        for chunk in iter(
            lambda: f.read(1024 * 1024),
            b"",
        ):
            h.update(chunk)

    return h.hexdigest()


manifest_sha = sha256(manifest_path)
fasta_sha = sha256(fasta_path)

lengths = [
    r["sequence_length"]
    for r in out_rows
]

metadata = {
    "snapshot": "20260101",
    "query_count": len(out_rows),
    "unique_openfold_keys": len(seen),
    "sequence_source": "PDBClean retained_sequence",
    "minimum_length": min(lengths),
    "maximum_length": max(lengths),
    "total_residues": sum(lengths),
    "observed_sequence_alphabet":
        "".join(sorted(alphabet)),
    "manifest":
        str(manifest_path),
    "manifest_sha256":
        manifest_sha,
    "fasta":
        str(fasta_path),
    "fasta_sha256":
        fasta_sha,
}

meta_path = msa / "msa_query_manifest.json"

meta_path.write_text(
    json.dumps(
        metadata,
        indent=2,
        sort_keys=True,
    ) + "\n"
)

success_path = msa / "_QUERIES_SUCCESS"

success_path.write_text(
    "\n".join([
        "STATUS=PASS",
        f"QUERY_COUNT={len(out_rows)}",
        f"MANIFEST_SHA256={manifest_sha}",
        f"FASTA_SHA256={fasta_sha}",
        "",
    ])
)

print()
print("===== OUTPUT =====")
print("queries:", len(out_rows))
print("FASTA:", fasta_path)
print("FASTA sha256:", fasta_sha)
print("manifest:", manifest_path)
print("manifest sha256:", manifest_sha)
print("minimum length:", min(lengths))
print("maximum length:", max(lengths))
print("total residues:", sum(lengths))

print()
print("STATUS: PASS")
print(
    "All 499770 authoritative PDBClean sequences "
    "are frozen as OpenFold MSA queries."
)
PY

echo
echo "Completed: $(date --iso-8601=seconds)"
