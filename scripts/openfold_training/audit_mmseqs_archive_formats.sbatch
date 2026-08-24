#!/bin/bash -l
#SBATCH --job-name=of_db_format
#SBATCH --partition=nodes
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:15:00

set -euo pipefail

DBROOT="${1:?database root required}"
PY="${2:?python executable required}"

DBROOT="$(readlink -f "$DBROOT")"
PY="$(readlink -f "$PY")"

ARCHIVES="$DBROOT/archives"

if [[ ! -d "$ARCHIVES" ]]; then
    echo "ERROR: archive directory missing: $ARCHIVES"
    exit 1
fi

if [[ ! -x "$PY" ]]; then
    echo "ERROR: Python missing: $PY"
    exit 1
fi

echo "============================================================"
echo "MMSEQS ARCHIVE FORMAT AUDIT"
echo "Job:      ${SLURM_JOB_ID}"
echo "Node:     $(hostname)"
echo "Root:     $DBROOT"
echo "Started:  $(date --iso-8601=seconds)"
echo "============================================================"

"$PY" - "$ARCHIVES" <<'PY'
from pathlib import Path
import collections
import sys
import tarfile

root = Path(sys.argv[1])

archives = sorted(root.glob("*.tar.gz"))

if len(archives) != 2:
    raise RuntimeError(
        f"Expected exactly 2 .tar.gz archives, found {len(archives)}"
    )

for archive in archives:
    print()
    print("=" * 72)
    print("ARCHIVE:", archive.name)
    print("SIZE_BYTES:", archive.stat().st_size)
    print("=" * 72)

    names = []
    suffixes = collections.Counter()

    # Streaming mode: do not scan or extract the entire archive.
    with tarfile.open(archive, mode="r|gz") as tf:
        for i, member in enumerate(tf):
            names.append(member.name)

            name = member.name.lower()

            for marker in (
                ".tsv",
                ".ffdata",
                ".ffindex",
                "_db",
                "_dbtype",
                ".index",
                "_seq",
                "_aln",
                "_h",
                ".fasta",
                ".fa",
            ):
                if marker in name:
                    suffixes[marker] += 1

            if i >= 199:
                break

    if not names:
        raise RuntimeError(f"No tar members found in {archive}")

    print()
    print("FIRST MEMBERS:")
    for name in names[:80]:
        print(name)

    print()
    print("MARKERS IN FIRST", len(names), "MEMBERS:")
    for marker, count in sorted(suffixes.items()):
        print(f"{marker}: {count}")

    joined = "\n".join(names).lower()

    looks_hhsuite = (
        ".ffdata" in joined
        or ".ffindex" in joined
    )

    looks_tsv = ".tsv" in joined

    looks_mmseqs = any(
        x in joined
        for x in (
            "_dbtype",
            "_seq",
            "_aln",
        )
    )

    print()
    print("FORMAT SIGNALS:")
    print("HHsuite-style files:", looks_hhsuite)
    print("TSV source files:", looks_tsv)
    print("MMseqs-style files:", looks_mmseqs)

print()
print("=" * 72)
print("STATUS: PASS")
print("Archive headers read successfully; nothing extracted.")
PY

echo
echo "Completed: $(date --iso-8601=seconds)"
