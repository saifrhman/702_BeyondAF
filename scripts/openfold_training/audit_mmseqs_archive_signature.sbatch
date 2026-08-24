#!/bin/bash -l
#SBATCH --job-name=of_db_sig
#SBATCH --partition=nodes
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:10:00

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
    echo "ERROR: Python executable missing: $PY"
    exit 1
fi

export PYTHONUNBUFFERED=1
export PYTHONDONTWRITEBYTECODE=1

echo "============================================================"
echo "MMSEQS ARCHIVE SIGNATURE AUDIT"
echo "Job:      ${SLURM_JOB_ID}"
echo "Node:     $(hostname)"
echo "Root:     $DBROOT"
echo "Started:  $(date --iso-8601=seconds)"
echo "============================================================"

"$PY" - "$ARCHIVES" <<'PY'
from pathlib import Path
import sys
import tarfile

root = Path(sys.argv[1])
archives = sorted(root.glob("*.tar.gz"))

if len(archives) != 2:
    raise RuntimeError(
        f"Expected exactly 2 archives, found {len(archives)}"
    )

for archive in archives:
    print()
    print("=" * 72)
    print("ARCHIVE:", archive.name, flush=True)
    print("COMPRESSED_BYTES:", archive.stat().st_size, flush=True)

    first_regular = None
    sample = b""
    seen = []

    with tarfile.open(archive, mode="r|gz") as tf:
        for member in tf:
            seen.append(
                (
                    member.name,
                    member.size,
                    member.type,
                )
            )

            # Directories and other zero-size metadata entries are cheap
            # to step over. Stop as soon as the first real file is found.
            if member.isfile():
                first_regular = member

                f = tf.extractfile(member)
                if f is not None:
                    sample = f.read(4096)

                break

    if first_regular is None:
        raise RuntimeError(
            f"No regular file found near start of {archive.name}"
        )

    print()
    print("MEMBERS ENCOUNTERED BEFORE FIRST FILE:")
    for name, size, kind in seen:
        print(
            f"name={name!r} size={size} type={kind!r}",
            flush=True,
        )

    name = first_regular.name
    lower = name.lower()

    print()
    print("FIRST REGULAR FILE:")
    print("name:", name)
    print("uncompressed_size:", first_regular.size)
    print("sample_bytes_read:", len(sample))

    nul_count = sample.count(b"\x00")
    tab_count = sample.count(b"\t")
    newline_count = sample.count(b"\n")

    printable = sum(
        1
        for b in sample
        if b in (9, 10, 13) or 32 <= b <= 126
    )

    printable_fraction = (
        printable / len(sample)
        if sample
        else 0.0
    )

    print()
    print("SAMPLE CHARACTERISTICS:")
    print("NUL bytes:", nul_count)
    print("TAB bytes:", tab_count)
    print("newlines:", newline_count)
    print(
        "printable_fraction:",
        f"{printable_fraction:.4f}",
    )

    print()
    print("FILENAME SIGNALS:")
    print(
        "HHsuite ffdata/ffindex:",
        ".ffdata" in lower
        or ".ffindex" in lower,
    )
    print(
        "MMseqs dbtype/index-style:",
        any(
            token in lower
            for token in (
                ".dbtype",
                "_dbtype",
                ".index",
                "_seq",
                "_aln",
            )
        ),
    )
    print(
        "TSV filename:",
        lower.endswith(".tsv"),
    )

    # A text/tabular source used by old ColabFold MMseqs preparation
    # should look overwhelmingly textual and normally contain tabs/newlines.
    looks_textual = (
        printable_fraction > 0.95
        and nul_count == 0
    )

    print(
        "sample_looks_textual:",
        looks_textual,
    )

    print()
    print("FIRST SAMPLE:")
    if looks_textual:
        print(
            sample.decode(
                "utf-8",
                errors="replace",
            )[:1000]
        )
    else:
        print(sample[:128].hex())

print()
print("=" * 72)
print("STATUS: PASS")
print(
    "Only the beginning of each archive was inspected; "
    "nothing was extracted."
)
PY

echo
echo "Completed: $(date --iso-8601=seconds)"
