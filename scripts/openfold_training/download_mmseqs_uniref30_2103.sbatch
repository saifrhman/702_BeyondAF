#!/bin/bash -l
#SBATCH --job-name=of_uniref_mm
#SBATCH --partition=nodes
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=1-00:00:00

set -euo pipefail

DBROOT="${1:?database root required}"
DBROOT="$(readlink -f "$DBROOT")"

ARCHIVES="$DBROOT/archives"
PROVENANCE="$DBROOT/provenance"

mkdir -p "$ARCHIVES" "$PROVENANCE"

URL="https://wwwuser.gwdg.de/~compbiol/colabfold/uniref30_2103.tar.gz"
NAME="uniref30_2103.tar.gz"

# Official checksum published by the ColabFold download site.
EXPECTED_MD5="2af35af92a9b1cf287a2845a80adf0c3"

FINAL="$ARCHIVES/$NAME"
PART="$FINAL.part"

for cmd in curl md5sum sha256sum stat tar awk; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "ERROR: missing required command: $cmd"
        exit 1
    fi
done

echo "============================================================"
echo "MMSEQS UNIREF30 2103 DOWNLOAD"
echo "Job:      ${SLURM_JOB_ID}"
echo "Node:     $(hostname)"
echo "Root:     $DBROOT"
echo "Started:  $(date --iso-8601=seconds)"
echo "============================================================"

echo
echo "===== REMOTE OBJECT ====="

HEADERS="$(mktemp)"
trap 'rm -f "$HEADERS"' EXIT

RESOLVED="$(
    curl \
        --silent \
        --show-error \
        --fail \
        --location \
        --head \
        --connect-timeout 30 \
        --output "$HEADERS" \
        --write-out '%{url_effective}' \
        "$URL"
)"

REMOTE_BYTES="$(
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
    ' "$HEADERS"
)"

if [[ -z "$REMOTE_BYTES" ]]; then
    echo "ERROR: remote server did not report Content-Length"
    exit 1
fi

echo "source URL:    $URL"
echo "resolved URL:  $RESOLVED"
echo "remote bytes:  $REMOTE_BYTES"
echo "expected MD5:  $EXPECTED_MD5"

echo
echo "===== DOWNLOAD ====="

if [[ -f "$FINAL" ]]; then
    LOCAL_BYTES="$(stat -c '%s' "$FINAL")"

    if [[ "$LOCAL_BYTES" != "$REMOTE_BYTES" ]]; then
        echo "ERROR: existing final archive has wrong size"
        echo "remote=$REMOTE_BYTES local=$LOCAL_BYTES"
        exit 1
    fi

    echo "Existing complete-sized archive found; transfer skipped."

else
    for attempt in $(seq 1 20); do

        if [[ -f "$PART" ]]; then
            HAVE="$(stat -c '%s' "$PART")"
        else
            HAVE=0
        fi

        echo "Attempt $attempt/20: $HAVE / $REMOTE_BYTES bytes"

        if [[ "$HAVE" == "$REMOTE_BYTES" ]]; then
            break
        fi

        if (( HAVE > REMOTE_BYTES )); then
            echo "ERROR: partial file exceeds remote object size"
            exit 1
        fi

        if (( HAVE > 0 )); then
            curl \
                --fail \
                --location \
                --show-error \
                --connect-timeout 30 \
                --speed-time 600 \
                --speed-limit 1024 \
                --continue-at - \
                --output "$PART" \
                "$URL" \
                || true
        else
            curl \
                --fail \
                --location \
                --show-error \
                --connect-timeout 30 \
                --speed-time 600 \
                --speed-limit 1024 \
                --output "$PART" \
                "$URL" \
                || true
        fi

        if [[ -f "$PART" ]]; then
            HAVE="$(stat -c '%s' "$PART")"
        else
            HAVE=0
        fi

        if [[ "$HAVE" == "$REMOTE_BYTES" ]]; then
            break
        fi

        sleep 15
    done

    HAVE="$(stat -c '%s' "$PART")"

    if [[ "$HAVE" != "$REMOTE_BYTES" ]]; then
        echo "ERROR: incomplete after retry budget"
        echo "remote=$REMOTE_BYTES local=$HAVE"
        exit 1
    fi

    mv "$PART" "$FINAL"
fi

echo
echo "===== CHECKSUM VERIFICATION ====="

ACTUAL_MD5="$(md5sum "$FINAL" | awk '{print $1}')"

echo "expected MD5: $EXPECTED_MD5"
echo "actual MD5:   $ACTUAL_MD5"

if [[ "$ACTUAL_MD5" != "$EXPECTED_MD5" ]]; then
    echo "ERROR: MD5 verification failed"
    exit 1
fi

SHA256="$(sha256sum "$FINAL" | awk '{print $1}')"

echo "SHA256:       $SHA256"

echo
echo "===== FORMAT VERIFICATION ====="

FIRST_MEMBER="$(
    tar -tzf "$FINAL" \
        | head -n 1
)"

echo "first archive member: $FIRST_MEMBER"

case "$FIRST_MEMBER" in
    *.tsv|*_h.tsv)
        ;;
    *)
        echo "ERROR: expected MMseqs/ColabFold TSV-format archive"
        echo "       but first member was: $FIRST_MEMBER"
        exit 1
        ;;
esac

{
    printf 'database\tuniref30_2103\n'
    printf 'source_url\t%s\n' "$URL"
    printf 'resolved_url\t%s\n' "$RESOLVED"
    printf 'bytes\t%s\n' "$(stat -c '%s' "$FINAL")"
    printf 'md5\t%s\n' "$ACTUAL_MD5"
    printf 'sha256\t%s\n' "$SHA256"
    printf 'first_member\t%s\n' "$FIRST_MEMBER"
    printf 'verified_at\t%s\n' "$(date --iso-8601=seconds)"
} > "$PROVENANCE/uniref30_2103.provenance.tsv"

echo
echo "============================================================"
echo "STATUS: PASS"
echo "Correct MMseqs-format UniRef30 2103 archive verified."
echo "No extraction performed."
echo "No MMseqs index created."
echo "No MSA generated."
echo "Completed: $(date --iso-8601=seconds)"
echo "============================================================"
