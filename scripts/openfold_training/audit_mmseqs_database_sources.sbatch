#!/bin/bash -l
#SBATCH --job-name=of_db_sources
#SBATCH --partition=nodes
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:10:00

set -euo pipefail

OPENFOLD="${1:?OpenFold root required}"
OPENFOLD="$(readlink -f "$OPENFOLD")"

UNIREF_SCRIPT="$OPENFOLD/scripts/download_uniref30.sh"
ENV_SCRIPT="$OPENFOLD/scripts/download_colabfold_envdb.sh"
PREP_SCRIPT="$OPENFOLD/scripts/prep_mmseqs_dbs.sh"

for f in "$UNIREF_SCRIPT" "$ENV_SCRIPT" "$PREP_SCRIPT"; do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: required OpenFold script missing: $f"
        exit 1
    fi
done

if ! command -v curl >/dev/null 2>&1; then
    echo "ERROR: curl unavailable"
    exit 1
fi

extract_source_url() {
    local file="$1"

    sed -n \
        's/^SOURCE_URL="\([^"]*\)".*/\1/p' \
        "$file" \
        | head -n 1
}

UNIREF_URL="$(extract_source_url "$UNIREF_SCRIPT")"
ENV_URL="$(extract_source_url "$ENV_SCRIPT")"

if [[ -z "$UNIREF_URL" ]]; then
    echo "ERROR: could not derive UniRef source URL"
    exit 1
fi

if [[ -z "$ENV_URL" ]]; then
    echo "ERROR: could not derive environmental DB source URL"
    exit 1
fi

echo "============================================================"
echo "MMSEQS DATABASE SOURCE PREFLIGHT"
echo "Job:       ${SLURM_JOB_ID}"
echo "Node:      $(hostname)"
echo "OpenFold:  $OPENFOLD"
echo "Started:   $(date --iso-8601=seconds)"
echo "============================================================"

check_url() {
    local label="$1"
    local url="$2"
    local headers
    local metrics

    headers="$(mktemp)"
    metrics="$(mktemp)"

    echo
    echo "===== $label ====="
    echo "archive: $(basename "$url")"

    if ! curl \
        --silent \
        --show-error \
        --fail \
        --location \
        --head \
        --retry 3 \
        --retry-delay 2 \
        --connect-timeout 30 \
        --max-time 120 \
        --output "$headers" \
        --write-out \
'http_code=%{http_code}
url_effective=%{url_effective}
content_type=%{content_type}
' \
        "$url" \
        > "$metrics"
    then
        echo "ERROR: HEAD request failed for $label"
        rm -f "$headers" "$metrics"
        return 1
    fi

    cat "$metrics"

    local length
    length="$(
        awk '
            BEGIN { IGNORECASE=1 }
            /^content-length:/ {
                gsub("\r", "", $2)
                value=$2
            }
            END {
                if (value != "")
                    print value
            }
        ' "$headers"
    )"

    if [[ -n "$length" ]]; then
        echo "content_length_bytes=$length"

        awk \
            -v n="$length" \
            'BEGIN {
                printf "content_length_GiB=%.3f\n",
                    n / 1024 / 1024 / 1024
            }'
    else
        echo "content_length_bytes=NOT_REPORTED"
    fi

    echo
    echo "final response headers:"
    tail -n 25 "$headers"

    rm -f "$headers" "$metrics"
}

check_url \
    "UNIREF30 SOURCE" \
    "$UNIREF_URL"

check_url \
    "ENVIRONMENTAL DB SOURCE" \
    "$ENV_URL"

echo
echo "===== PREPARATION PROTOCOL ====="

grep -E \
    'tar |tsv2exprofiledb|createindex' \
    "$PREP_SCRIPT"

echo
echo "===== SCRATCH STORAGE ====="
df -h "$HOME/scratch"

echo
echo "============================================================"
echo "STATUS: PASS"
echo "Completed: $(date --iso-8601=seconds)"
echo "============================================================"
