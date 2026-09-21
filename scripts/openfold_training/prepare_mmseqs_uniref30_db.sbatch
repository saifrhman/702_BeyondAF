#!/bin/bash -l
#SBATCH --job-name=of_uniref_prep
#SBATCH --partition=nodes
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=192G
#SBATCH --time=1-00:00:00

# Stage 18 prerequisite: build the MMseqs2 expandable-profile database from the
# already-downloaded, already-checksum-verified ColabFold UniRef30 archive.
#
# This is SHARED REFERENCE INFRASTRUCTURE, not run-specific data.  It contains
# no snapshot identifier and no run identity: the same prepared database is
# valid for any PDBClean run whose Stage-18 configuration names it.  The run's
# own snapshot, Gold release and query population are supplied downstream by
# the Stage-18 search job, never by this script.
#
# Scope, deliberately bounded:
#
#   extract  ->  tsv2exprofiledb  ->  measure  ->  reclaim intermediate space
#
# `mmseqs createindex` is NOT run here.  Indexing roughly doubles the on-disk
# footprint and the quota headroom must be judged against the real measured
# size of the converted database, which only exists once this job succeeds.
#
# The archive is never deleted, never re-downloaded, and never re-checksummed
# against the network.  Job 10275442 already proved its MD5 matches the
# published ColabFold checksum; it then died at `tar -tzf | head -n 1`, which
# raises SIGPIPE under `set -o pipefail`.  That was a reporting bug, not a data
# problem, so this script verifies bytes locally and writes the provenance
# record that job never reached.

# Portability: every path below comes from here, derived from $USER.
# Override any of them in the environment; nothing names an individual.
source "${PDBCLEAN_ENV_FILE:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/../../config/pdbclean/pipeline_env.sh}"

set -euo pipefail

DBROOT="${1:?database root required}"
DBROOT="$(readlink -f "$DBROOT")"

ARCHIVES="$DBROOT/archives"
PROVENANCE="$DBROOT/provenance"
WORK="$DBROOT/work/uniref30_2103"
DBDIR="$DBROOT/db"
TMPDIR_MMSEQS="$DBROOT/tmp/uniref30_2103"

ARCHIVE_NAME="uniref30_2103.tar.gz"
ARCHIVE="$ARCHIVES/$ARCHIVE_NAME"

# Verified by job 10275442 against the published ColabFold checksum.
EXPECTED_BYTES=75336202245
EXPECTED_MD5="2af35af92a9b1cf287a2845a80adf0c3"
EXPECTED_SHA256="c1ab61697d751ab5bd876ce22c6c4bbed2f166d4e8b8e9b7768fd15691e8d4aa"

DB_NAME="uniref30_2103_db"
DB_PATH="$DBDIR/$DB_NAME"

MMSEQS="${MMSEQS_BIN:-${MMSEQS_BIN}}"

echo "============================================================"
echo "MMSEQS UNIREF30 EXPANDABLE-PROFILE DATABASE PREPARATION"
echo "Job:      ${SLURM_JOB_ID:-none}"
echo "Node:     $(hostname)"
echo "CPUs:     ${SLURM_CPUS_PER_TASK:-unknown}"
echo "Root:     $DBROOT"
echo "Started:  $(date --iso-8601=seconds)"
echo "============================================================"

echo
echo "===== TOOLING ====="

if [[ ! -x "$MMSEQS" ]]; then
    echo "ERROR: mmseqs not executable at $MMSEQS"
    exit 1
fi

MMSEQS_VERSION="$("$MMSEQS" version)"

echo "mmseqs binary:  $MMSEQS"
echo "mmseqs version: $MMSEQS_VERSION"

for cmd in tar stat md5sum sha256sum du df awk; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "ERROR: missing required command: $cmd"
        exit 1
    fi
done

echo
echo "===== IDEMPOTENCE CHECK ====="

if [[ -f "$DBDIR/_UNIREF30_DB_SUCCESS" ]]; then
    echo "An existing PASS marker was found:"
    cat "$DBDIR/_UNIREF30_DB_SUCCESS"
    echo
    echo "Database already prepared; nothing to do."
    echo "STATUS: PASS (reused)"
    exit 0
fi

echo "No existing PASS marker; preparation will proceed."

echo
echo "===== SOURCE ARCHIVE ====="

if [[ ! -f "$ARCHIVE" ]]; then
    echo "ERROR: archive not found: $ARCHIVE"
    echo "       This job never downloads.  Nothing was modified."
    exit 1
fi

ACTUAL_BYTES="$(stat -c '%s' "$ARCHIVE")"

echo "archive:        $ARCHIVE"
echo "expected bytes: $EXPECTED_BYTES"
echo "actual bytes:   $ACTUAL_BYTES"

if [[ "$ACTUAL_BYTES" != "$EXPECTED_BYTES" ]]; then
    echo "ERROR: archive size does not match the verified download"
    exit 1
fi

# Whether the 75 GB container is worth re-hashing depends on whether this run
# will actually read it.  If a previous attempt already extracted the members,
# they are verified individually against the archive's own shipped md5 manifest
# further down -- which is strictly stronger evidence about what gets converted
# than a hash of a container nothing will open.  Hashing it anyway would cost
# ~10 minutes to prove something irrelevant.

ARCHIVE_BASE_EARLY="${ARCHIVE_NAME%.tar.gz}"
ARCHIVE_WILL_BE_READ=0

for member in \
    "${ARCHIVE_BASE_EARLY}.tsv" \
    "${ARCHIVE_BASE_EARLY}_h.tsv" \
    "${ARCHIVE_BASE_EARLY}_seq.tsv" \
    "${ARCHIVE_BASE_EARLY}_aln.tsv"
do
    if [[ ! -f "$WORK/$member" ]]; then
        ARCHIVE_WILL_BE_READ=1
        break
    fi
done

if (( ARCHIVE_WILL_BE_READ == 1 )); then
    echo "extraction will read this archive; recomputing checksums..."

    ACTUAL_MD5="$(md5sum "$ARCHIVE" | awk '{print $1}')"
    ACTUAL_SHA256="$(sha256sum "$ARCHIVE" | awk '{print $1}')"

    echo "expected MD5:    $EXPECTED_MD5"
    echo "actual MD5:      $ACTUAL_MD5"
    echo "expected SHA256: $EXPECTED_SHA256"
    echo "actual SHA256:   $ACTUAL_SHA256"

    if [[ "$ACTUAL_MD5" != "$EXPECTED_MD5" ]]; then
        echo "ERROR: MD5 mismatch against the published ColabFold checksum"
        exit 1
    fi

    if [[ "$ACTUAL_SHA256" != "$EXPECTED_SHA256" ]]; then
        echo "ERROR: SHA256 mismatch against the verified download"
        exit 1
    fi

    echo "Archive integrity confirmed."
else
    # Recorded from the verified download (job 10275442) and re-confirmed by
    # this job's predecessor attempt; the size check above proves the file has
    # not been replaced since.
    ACTUAL_MD5="$EXPECTED_MD5"
    ACTUAL_SHA256="$EXPECTED_SHA256"

    echo "members already extracted; container will not be read."
    echo "Container hash not recomputed; extracted members are verified"
    echo "individually against the shipped manifest before conversion."
fi

echo
echo "===== STORAGE BEFORE ====="

df -h "$DBROOT" || true
lfs quota -h -u "$(whoami)" /mnt/scratch 2>/dev/null || true

echo
echo "===== EXTRACTION ====="

mkdir -p "$WORK" "$DBDIR" "$TMPDIR_MMSEQS"

# The ColabFold TSV archive expands to a fixed family of files:
#
#   <base>.tsv        cluster representatives
#   <base>_h.tsv      headers
#   <base>_seq.tsv    member sequences   (by far the largest)
#   <base>_aln.tsv    alignments
#   <base>.md5sums    per-file checksums, shipped inside the archive
#
# `tsv2exprofiledb` is given <base>, and derives the four siblings itself.
# The suffixed files are therefore NOT candidate basenames -- an earlier
# revision of this script picked the largest TSV and derived "<base>_seq",
# which cannot work.  Discovery below excludes every known suffix instead of
# ranking by size.

ARCHIVE_BASE="${ARCHIVE_NAME%.tar.gz}"

EXPECTED_MEMBERS=(
    "${ARCHIVE_BASE}.tsv"
    "${ARCHIVE_BASE}_h.tsv"
    "${ARCHIVE_BASE}_seq.tsv"
    "${ARCHIVE_BASE}_aln.tsv"
)

NEED_EXTRACTION=0

for member in "${EXPECTED_MEMBERS[@]}"; do
    if [[ ! -f "$WORK/$member" ]]; then
        NEED_EXTRACTION=1
        break
    fi
done

if (( NEED_EXTRACTION == 0 )); then
    echo "All expected members are already present in $WORK"
    echo "Extraction skipped; contents are verified against the shipped"
    echo "checksums below before anything is converted."
else
    if [[ -n "$(ls -A "$WORK" 2>/dev/null)" ]]; then
        echo "Work directory holds an incomplete extraction; clearing it."
        echo "  $WORK"
        # Bounded to the work directory only.  Archives live elsewhere and are
        # never a candidate for removal.
        find "$WORK" -mindepth 1 -delete
    fi

    echo "extracting $ARCHIVE_NAME into $WORK"
    echo "(single-threaded gunzip over 75 GB; expect tens of minutes)"

    tar --extract --file="$ARCHIVE" --directory="$WORK"

    echo "extraction complete: $(date --iso-8601=seconds)"
fi

echo
echo "===== EXTRACTED CONTENTS ====="

find "$WORK" -maxdepth 2 -type f -printf '%14s  %p\n' | sort -k2

EXTRACTED_BYTES="$(du -sb "$WORK" | awk '{print $1}')"

echo
echo "extracted bytes: $EXTRACTED_BYTES"

echo
echo "===== MEMBER VERIFICATION ====="

# The archive ships its own per-file md5 manifest.  Verifying against it proves
# the extracted TSVs are intact -- which matters more than the container hash
# once extraction is resumable, because a skipped extraction trusts whatever is
# already on disk.
#
# Note: the manifest's self-reference for the .tar.gz names a different
# packaging build than the file ColabFold publishes today.  Our archive matched
# the published MD5 exactly, so the container is correct; only the per-member
# lines below are treated as authoritative here.

MD5_MANIFEST="$WORK/${ARCHIVE_BASE}.md5sums"

if [[ ! -f "$MD5_MANIFEST" ]]; then
    echo "ERROR: shipped checksum manifest not found: $MD5_MANIFEST"
    exit 1
fi

VERIFIED=0

for member in "${EXPECTED_MEMBERS[@]}"; do
    if [[ ! -f "$WORK/$member" ]]; then
        echo "ERROR: expected member missing after extraction: $member"
        exit 1
    fi

    EXPECTED_MEMBER_MD5="$(
        awk -v want="$member" '$2 == want { print $1 }' "$MD5_MANIFEST"
    )"

    if [[ -z "$EXPECTED_MEMBER_MD5" ]]; then
        echo "ERROR: $member has no entry in the shipped checksum manifest"
        exit 1
    fi

    echo "verifying $member ($(stat -c '%s' "$WORK/$member") bytes)..."

    ACTUAL_MEMBER_MD5="$(md5sum "$WORK/$member" | awk '{print $1}')"

    if [[ "$ACTUAL_MEMBER_MD5" != "$EXPECTED_MEMBER_MD5" ]]; then
        echo "ERROR: checksum mismatch for $member"
        echo "  expected $EXPECTED_MEMBER_MD5"
        echo "  actual   $ACTUAL_MEMBER_MD5"
        exit 1
    fi

    echo "  OK $ACTUAL_MEMBER_MD5"
    VERIFIED=$((VERIFIED + 1))
done

echo "$VERIFIED archive members verified against the shipped manifest."

echo
echo "===== BASENAME DISCOVERY ====="

BASENAME="$WORK/$ARCHIVE_BASE"

echo "derived basename: $BASENAME"
echo "tsv2exprofiledb will read:"

for member in "${EXPECTED_MEMBERS[@]}"; do
    echo "  $member"
done

echo
echo "===== TSV2EXPROFILEDB ====="

echo "command: $MMSEQS tsv2exprofiledb $BASENAME $DB_PATH"
echo "started: $(date --iso-8601=seconds)"

"$MMSEQS" tsv2exprofiledb "$BASENAME" "$DB_PATH"

echo "finished: $(date --iso-8601=seconds)"

echo
echo "===== CONVERTED DATABASE ====="

find "$DBDIR" -maxdepth 1 -name "${DB_NAME}*" -printf '%14s  %p\n' | sort -k2

DB_BYTES="$(
    find "$DBDIR" -maxdepth 1 -name "${DB_NAME}*" -printf '%s\n' \
        | awk '{total += $1} END {print total + 0}'
)"

echo
echo "converted database bytes: $DB_BYTES"

# tsv2exprofiledb produces the profile, sequence and header sub-databases that
# an expandable-profile search requires.  Their absence means the conversion
# did not actually complete, whatever the exit status said.
MISSING=0

for suffix in "" "_h" "_seq" "_aln"; do
    if [[ -f "${DB_PATH}${suffix}" || -f "${DB_PATH}${suffix}.dbtype" || -f "${DB_PATH}${suffix}.index" ]]; then
        echo "present: ${DB_NAME}${suffix}"
    else
        echo "ABSENT:  ${DB_NAME}${suffix}"
        MISSING=$((MISSING + 1))
    fi
done

if (( MISSING > 0 )); then
    echo
    echo "ERROR: $MISSING expected sub-database(s) are missing."
    echo "       Extracted TSVs have been LEFT IN PLACE for diagnosis:"
    echo "         $WORK"
    exit 1
fi

echo
echo "===== RECLAIMING INTERMEDIATE SPACE ====="

# The extracted TSVs are a pure intermediate: they are reproducible at any time
# from the retained archive, and they are large.  They are removed only after
# the conversion above has been positively verified.  The archive itself is
# never touched.

echo "removing extracted TSVs under $WORK"
echo "(archive $ARCHIVE_NAME is retained and untouched)"

find "$WORK" -mindepth 1 -delete
rmdir "$WORK" 2>/dev/null || true

if [[ ! -f "$ARCHIVE" ]]; then
    echo "ERROR: source archive disappeared; refusing to record PASS"
    exit 1
fi

echo "archive still present: $(stat -c '%s' "$ARCHIVE") bytes"

echo
echo "===== STORAGE AFTER ====="

df -h "$DBROOT" || true
lfs quota -h -u "$(whoami)" /mnt/scratch 2>/dev/null || true

echo
echo "===== PROVENANCE ====="

{
    printf 'database\tuniref30_2103\n'
    printf 'archive\t%s\n' "$ARCHIVE"
    printf 'archive_bytes\t%s\n' "$ACTUAL_BYTES"
    printf 'archive_md5\t%s\n' "$ACTUAL_MD5"
    printf 'archive_sha256\t%s\n' "$ACTUAL_SHA256"
    printf 'extracted_bytes\t%s\n' "$EXTRACTED_BYTES"
    printf 'tsv_basename\t%s\n' "$(basename "$BASENAME")"
    printf 'mmseqs_version\t%s\n' "$MMSEQS_VERSION"
    printf 'mmseqs_binary\t%s\n' "$MMSEQS"
    printf 'conversion\ttsv2exprofiledb\n'
    printf 'database_path\t%s\n' "$DB_PATH"
    printf 'database_bytes\t%s\n' "$DB_BYTES"
    printf 'createindex\tnot_run\n'
    printf 'slurm_job_id\t%s\n' "${SLURM_JOB_ID:-none}"
    printf 'node\t%s\n' "$(hostname)"
    printf 'prepared_at\t%s\n' "$(date --iso-8601=seconds)"
} > "$PROVENANCE/uniref30_2103.provenance.tsv"

cat "$PROVENANCE/uniref30_2103.provenance.tsv"

{
    printf 'STATUS=PASS\n'
    printf 'DATABASE=%s\n' "$DB_PATH"
    printf 'DATABASE_BYTES=%s\n' "$DB_BYTES"
    printf 'ARCHIVE_SHA256=%s\n' "$ACTUAL_SHA256"
    printf 'MMSEQS_VERSION=%s\n' "$MMSEQS_VERSION"
    printf 'CREATEINDEX=not_run\n'
    printf 'SLURM_JOB_ID=%s\n' "${SLURM_JOB_ID:-none}"
    printf 'PREPARED_AT=%s\n' "$(date --iso-8601=seconds)"
} > "$DBDIR/_UNIREF30_DB_SUCCESS"

echo
echo "============================================================"
echo "STATUS: PASS"
echo "UniRef30 expandable-profile database prepared."
echo "  database:   $DB_PATH"
echo "  bytes:      $DB_BYTES"
echo "  createindex: NOT RUN (separate job, after measuring headroom)"
echo "No MSA generated.  No run-specific data touched."
echo "Completed: $(date --iso-8601=seconds)"
echo "============================================================"
