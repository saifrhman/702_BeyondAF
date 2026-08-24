#!/bin/bash -l
#SBATCH --job-name=of_roda_audit
#SBATCH --partition=nodes
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:30:00

set -euo pipefail

TRAIN_ROOT="${1:?Usage: sbatch audit_roda_remote.sbatch <training_root> <config> <python>}"
CONFIG="${2:?Usage: sbatch audit_roda_remote.sbatch <training_root> <config> <python>}"
PY="${3:?Usage: sbatch audit_roda_remote.sbatch <training_root> <config> <python>}"

TRAIN_ROOT="$(readlink -f "$TRAIN_ROOT")"
CONFIG="$(readlink -f "$CONFIG")"
PY="$(readlink -f "$PY")"

if [[ ! -f "$CONFIG" ]]; then
    echo "ERROR: missing config: $CONFIG"
    exit 1
fi

if [[ ! -x "$PY" ]]; then
    echo "ERROR: missing Python executable: $PY"
    exit 1
fi

source "$CONFIG"

: "${RODA_S3_SRC:?RODA_S3_SRC is not defined}"
: "${OPENFOLD_CODE_DIR:?OPENFOLD_CODE_DIR is not defined}"

OUT="$TRAIN_ROOT/msa/audit"
mkdir -p "$OUT"

export PYTHONUNBUFFERED=1
export PYTHONDONTWRITEBYTECODE=1

echo "============================================================"
echo "RODA / OPENPROTEINSET REMOTE AUDIT"
echo "Job:      ${SLURM_JOB_ID}"
echo "Node:     $(hostname)"
echo "Source:   $RODA_S3_SRC"
echo "Started:  $(date --iso-8601=seconds)"
echo "============================================================"

echo
echo "===== PUBLIC S3 LISTING ====="

"$PY" - "$RODA_S3_SRC" <<'PY'
import sys
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

uri = sys.argv[1]

parsed = urllib.parse.urlparse(uri)

if parsed.scheme != "s3":
    raise RuntimeError(
        f"Expected s3:// URI, got {uri}"
    )

bucket = parsed.netloc
prefix = parsed.path.lstrip("/").rstrip("/")

if prefix:
    prefix += "/"

params = urllib.parse.urlencode({
    "list-type": "2",
    "prefix": prefix,
    "delimiter": "/",
    "max-keys": "100",
})

urls = [
    f"https://{bucket}.s3.amazonaws.com/?{params}",
    f"https://s3.amazonaws.com/{bucket}?{params}",
]

last_error = None
data = None
used_url = None

for url in urls:
    try:
        req = urllib.request.Request(
            url,
            headers={
                "User-Agent": "COMP702-RODA-audit/1.0",
            },
        )

        with urllib.request.urlopen(
            req,
            timeout=60,
        ) as response:
            data = response.read()
            used_url = response.geturl()

        break

    except Exception as exc:
        last_error = exc

if data is None:
    raise RuntimeError(
        f"Unable to list public S3 bucket: {last_error}"
    )

print("bucket:", bucket)
print("prefix:", prefix)
print("resolved URL:", used_url)
print()

root = ET.fromstring(data)

# S3 XML responses normally use a namespace.
def text_ending(elem, suffix):
    for child in elem.iter():
        if child.tag.endswith(suffix):
            return child.text
    return None

prefixes = []

for elem in root.iter():
    if elem.tag.endswith("CommonPrefixes"):
        p = text_ending(elem, "Prefix")
        if p:
            prefixes.append(p)

objects = []

for elem in root.iter():
    if elem.tag.endswith("Contents"):
        key = text_ending(elem, "Key")
        size = text_ending(elem, "Size")

        if key:
            objects.append(
                (
                    key,
                    int(size) if size else None,
                )
            )

print("Common prefixes:")
if prefixes:
    for p in prefixes:
        print(p)
else:
    print("(none)")

print()
print("Objects directly under prefix:")
if objects:
    for key, size in objects:
        print(key, size)
else:
    print("(none)")

print()
print("prefix_count:", len(prefixes))
print("object_count:", len(objects))
PY

echo
echo "===== OPENFOLD RODA HANDLING ====="

grep -RniE \
    'RODA|OpenProteinSet|duplicate.*chain|100% sequence|flatten_roda' \
    "$OPENFOLD_CODE_DIR/scripts" \
    2>/dev/null \
    | head -n 200 \
    || true

echo
echo "============================================================"
echo "RODA REMOTE AUDIT COMPLETE"
echo "Completed: $(date --iso-8601=seconds)"
echo "============================================================"
