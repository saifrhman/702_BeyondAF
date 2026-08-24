#!/bin/bash -l
#SBATCH --job-name=of_mmseqs_db_dl
#SBATCH --partition=nodes
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=1-00:00:00

set -euo pipefail

OPENFOLD="${1:?Usage: sbatch download_mmseqs_reference_dbs.sbatch <openfold_root> <db_root>}"
DBROOT="${2:?Usage: sbatch download_mmseqs_reference_dbs.sbatch <openfold_root> <db_root>}"

OPENFOLD="$(readlink -f "$OPENFOLD")"

mkdir -p "$DBROOT"
DBROOT="$(readlink -f "$DBROOT")"

ARCHIVES="$DBROOT/archives"
PROVENANCE="$DBROOT/provenance"

mkdir -p "$ARCHIVES" "$PROVENANCE"

UNIREF_SCRIPT="$OPENFOLD/scripts/download_uniref30.sh"
ENV_SCRIPT="$OPENFOLD/scripts/download_colabfold_envdb.sh"

for f in "$UNIREF_SCRIPT" "$ENV_SCRIPT"; do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: required OpenFold source script missing: $f"
        exit 1
    fi
done

for cmd in curl sha256sum stat awk df flock; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "ERROR: required command unavailable: $cmd"
        exit 1
    fi
done

# Prevent two database-install jobs writing the same files.
exec 9>"$DBROOT/.download.lock"

if ! flock -n 9; then
    echo "ERROR: another MMseqs database download is already using:"
    echo "       $DBROOT"
    exit 1
fi

extract_source_url() {
    local script="$1"

    sed -n \
        's/^SOURCE_URL="\([^"]*\)".*/\1/p' \
        "$script" \
        | head -n 1
}

UNIREF_URL="$(extract_source_url "$UNIREF_SCRIPT")"
ENV_URL="$(extract_source_url "$ENV_SCRIPT")"

if [[ -z "$UNIREF_URL" || -z "$ENV_URL" ]]; then
    echo "ERROR: failed to derive database URLs from OpenFold scripts"
    exit 1
fi

probe_url() {
    local url="$1"
    local headers="$2"
    local resolved
    local bytes

    resolved="$(
        curl \
            --silent \
            --show-error \
            --fail \
            --location \
            --head \
            --connect-timeout 30 \
            --output "$headers" \
            --write-out '%{url_effective}' \
            "$url"
    )"

    bytes="$(
        awk '
            BEGIN { IGNORECASE=1 }
            /^content-length:/ {
                gsub("\r", "", $2)
                if ($2 ~ /^[0-9]+$/)
                    value=$2
            }
            END {
                if (value != "")
                    print value
            }
        ' "$headers"
    )"

    if [[ -z "$bytes" ]]; then
        echo "ERROR: server did not report Content-Length for $url" >&2
        return 1
    fi

    printf '%s\t%s\n' "$bytes" "$resolved"
}

UNIREF_ARCHIVE="$(basename "$UNIREF_URL")"
ENV_ARCHIVE="$(basename "$ENV_URL")"

IFS=$'\t' read -r UNIREF_BYTES UNIREF_RESOLVED < <(
    probe_url \
        "$UNIREF_URL" \
        "$PROVENANCE/${UNIREF_ARCHIVE}.headers.txt"
)

IFS=$'\t' read -r ENV_BYTES ENV_RESOLVED < <(
    probe_url \
        "$ENV_URL" \
        "$PROVENANCE/${ENV_ARCHIVE}.headers.txt"
)

EXPECTED_TOTAL=$((UNIREF_BYTES + ENV_BYTES))

# Download stage only: require 10% headroom over remaining source payload.
AVAILABLE_BYTES="$(
    df --output=avail -B1 "$DBROOT" \
        | tail -n 1 \
        | tr -d ' '
)"

if [[ ! "$AVAILABLE_BYTES" =~ ^[0-9]+$ ]]; then
    echo "ERROR: could not determine available scratch space"
    exit 1
fi

REQUIRED_BYTES=$((EXPECTED_TOTAL + EXPECTED_TOTAL / 10))

echo "============================================================"
echo "MMSEQS REFERENCE DATABASE DOWNLOAD"
echo "Job:          ${SLURM_JOB_ID}"
echo "Node:         $(hostname)"
echo "OpenFold:     $OPENFOLD"
echo "Database root:$DBROOT"
echo "Started:      $(date --iso-8601=seconds)"
echo "============================================================"

echo
echo "===== SOURCE DATABASES ====="
echo "UniRef archive:       $UNIREF_ARCHIVE"
echo "UniRef bytes:         $UNIREF_BYTES"
echo "UniRef resolved URL:  $UNIREF_RESOLVED"
echo
echo "Environment archive: $ENV_ARCHIVE"
echo "Environment bytes:   $ENV_BYTES"
echo "Environment URL:     $ENV_RESOLVED"
echo
echo "Expected total bytes:$EXPECTED_TOTAL"
echo "Available bytes:     $AVAILABLE_BYTES"
echo "Required with margin:$REQUIRED_BYTES"

if (( AVAILABLE_BYTES < REQUIRED_BYTES )); then
    echo "ERROR: insufficient scratch space for source archives"
    exit 1
fi

download_one() {
    local label="$1"
    local source_url="$2"
    local resolved_url="$3"
    local expected_bytes="$4"
    local archive_name="$5"
    local source_script="$6"

    local final="$ARCHIVES/$archive_name"
    local part="$final.part"
    local actual
    local attempt
    local rc
    local sha
    local source_script_sha

    echo
    echo "============================================================"
    echo "DOWNLOADING: $label"
    echo "Archive:     $archive_name"
    echo "Expected:    $expected_bytes bytes"
    echo "============================================================"

    if [[ -f "$final" ]]; then
        actual="$(stat -c '%s' "$final")"

        if [[ "$actual" != "$expected_bytes" ]]; then
            echo "ERROR: completed archive exists with wrong size:"
            echo "       $final"
            echo "       expected=$expected_bytes actual=$actual"
            exit 1
        fi

        echo "Complete archive already present; download skipped."

    else
        if [[ -f "$part" ]]; then
            actual="$(stat -c '%s' "$part")"

            if (( actual > expected_bytes )); then
                echo "ERROR: partial archive is larger than remote object:"
                echo "       $part"
                echo "       expected=$expected_bytes actual=$actual"
                exit 1
            fi
        fi

        for attempt in $(seq 1 20); do
            if [[ -f "$part" ]]; then
                actual="$(stat -c '%s' "$part")"
            else
                actual=0
            fi

            if [[ "$actual" == "$expected_bytes" ]]; then
                break
            fi

            echo
            echo "Attempt $attempt/20"
            echo "Current bytes: $actual / $expected_bytes"

            if (( actual > 0 )); then
                if curl \
                    --fail \
                    --location \
                    --show-error \
                    --connect-timeout 30 \
                    --speed-time 600 \
                    --speed-limit 1024 \
                    --continue-at - \
                    --output "$part" \
                    "$source_url"
                then
                    rc=0
                else
                    rc=$?
                fi
            else
                if curl \
                    --fail \
                    --location \
                    --show-error \
                    --connect-timeout 30 \
                    --speed-time 600 \
                    --speed-limit 1024 \
                    --output "$part" \
                    "$source_url"
                then
                    rc=0
                else
                    rc=$?
                fi
            fi

            if [[ -f "$part" ]]; then
                actual="$(stat -c '%s' "$part")"
            else
                actual=0
            fi

            echo "Bytes after attempt: $actual / $expected_bytes"

            if [[ "$actual" == "$expected_bytes" ]]; then
                break
            fi

            echo "Transfer incomplete (curl rc=$rc)."
            echo "Waiting 15 seconds before resuming..."
            sleep 15
        done

        if [[ ! -f "$part" ]]; then
            echo "ERROR: no partial archive was created"
            exit 1
        fi

        actual="$(stat -c '%s' "$part")"

        if [[ "$actual" != "$expected_bytes" ]]; then
            echo "ERROR: archive incomplete after retry budget"
            echo "       expected=$expected_bytes actual=$actual"
            echo "       partial file retained for a future resume:"
            echo "       $part"
            exit 1
        fi

        mv "$part" "$final"
    fi

    echo
    echo "Computing SHA256..."
    sha="$(sha256sum "$final" | awk '{print $1}')"
    source_script_sha="$(sha256sum "$source_script" | awk '{print $1}')"

    printf '%s  %s\n' \
        "$sha" \
        "$final" \
        > "$PROVENANCE/${archive_name}.sha256"

    {
        printf 'label\t%s\n' "$label"
        printf 'archive\t%s\n' "$archive_name"
        printf 'source_url\t%s\n' "$source_url"
        printf 'resolved_url\t%s\n' "$resolved_url"
        printf 'expected_bytes\t%s\n' "$expected_bytes"
        printf 'local_bytes\t%s\n' "$(stat -c '%s' "$final")"
        printf 'sha256\t%s\n' "$sha"
        printf 'source_script\t%s\n' "$source_script"
        printf 'source_script_sha256\t%s\n' "$source_script_sha"
        printf 'verified_at\t%s\n' "$(date --iso-8601=seconds)"
    } > "$PROVENANCE/${archive_name}.provenance.tsv"

    echo "STATUS: $label DOWNLOAD PASS"
    echo "SHA256: $sha"
}

download_one \
    "UniRef30" \
    "$UNIREF_URL" \
    "$UNIREF_RESOLVED" \
    "$UNIREF_BYTES" \
    "$UNIREF_ARCHIVE" \
    "$UNIREF_SCRIPT"

download_one \
    "ColabFold environmental database" \
    "$ENV_URL" \
    "$ENV_RESOLVED" \
    "$ENV_BYTES" \
    "$ENV_ARCHIVE" \
    "$ENV_SCRIPT"

MANIFEST_TMP="$PROVENANCE/download_manifest.tsv.tmp"
MANIFEST="$PROVENANCE/download_manifest.tsv"

printf '%s\n' \
    $'archive\texpected_bytes\tsha256\tprovenance_file' \
    > "$MANIFEST_TMP"

for archive in "$UNIREF_ARCHIVE" "$ENV_ARCHIVE"; do
    bytes="$(
        awk -F '\t' \
            '$1 == "expected_bytes" {print $2}' \
            "$PROVENANCE/${archive}.provenance.tsv"
    )"

    sha="$(
        awk -F '\t' \
            '$1 == "sha256" {print $2}' \
            "$PROVENANCE/${archive}.provenance.tsv"
    )"

    printf '%s\t%s\t%s\t%s\n' \
        "$archive" \
        "$bytes" \
        "$sha" \
        "$PROVENANCE/${archive}.provenance.tsv" \
        >> "$MANIFEST_TMP"
done

mv "$MANIFEST_TMP" "$MANIFEST"
sha256sum "$MANIFEST" > "$PROVENANCE/download_manifest.tsv.sha256"

echo
echo "===== FINAL SOURCE ARCHIVES ====="
ls -lh "$ARCHIVES"

echo
echo "===== DOWNLOAD MANIFEST ====="
cat "$MANIFEST"

echo
echo "============================================================"
echo "STATUS: PASS"
echo "Both MMseqs reference database archives are complete,"
echo "byte-verified, SHA256-hashed, and preserved."
echo "No MSA has been generated yet."
echo "No database index has been built yet."
echo "Completed: $(date --iso-8601=seconds)"
echo "============================================================"
