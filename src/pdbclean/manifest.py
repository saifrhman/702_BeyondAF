"""Creation and validation of immutable PDB snapshot manifests."""

from __future__ import annotations

import hashlib
import re
from dataclasses import dataclass
from pathlib import Path

import pyarrow as pa
import pyarrow.csv as pacsv
import pyarrow.parquet as pq

from pdbclean.schemas import SOURCE_MANIFEST_SCHEMA


PDB_ID_PATTERN = re.compile(r"^[a-z0-9]{4}$")


class ManifestError(ValueError):
    """Raised when a source manifest fails validation."""


@dataclass(frozen=True)
class ManifestSummary:
    """Validated source-manifest statistics."""

    row_count: int
    total_bytes: int
    unique_pdb_ids: int
    unique_s3_keys: int
    sha256: str


def file_sha256(path: str | Path) -> str:
    """Return the SHA-256 checksum of a file."""

    file_path = Path(path)
    digest = hashlib.sha256()

    with file_path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)

    return digest.hexdigest()


def validate_manifest_table(
    table: pa.Table,
    *,
    expected_snapshot: str,
    expected_count: int | None = None,
    expected_total_bytes: int | None = None,
) -> ManifestSummary:
    """Validate manifest contents and return summary statistics."""

    required_columns = {
        "snapshot",
        "source_layout",
        "pdb_id",
        "s3_key",
        "size_bytes",
        "etag",
    }

    missing_columns = required_columns.difference(table.column_names)

    if missing_columns:
        raise ManifestError(
            "Manifest is missing required columns: "
            + ", ".join(sorted(missing_columns))
        )

    rows = table.select(sorted(required_columns)).to_pylist()

    seen_pdb_ids: set[str] = set()
    seen_s3_keys: set[str] = set()
    total_bytes = 0

    expected_prefix = (
        f"{expected_snapshot}/pub/pdb/data/structures/"
        "divided/mmCIF/"
    )

    for index, row in enumerate(rows, start=1):
        snapshot = str(row["snapshot"])
        source_layout = str(row["source_layout"])
        pdb_id = str(row["pdb_id"]).lower()
        s3_key = str(row["s3_key"])
        etag = str(row["etag"]).strip()
        size_bytes = row["size_bytes"]

        if snapshot != expected_snapshot:
            raise ManifestError(
                f"Row {index}: snapshot {snapshot!r} does not match "
                f"{expected_snapshot!r}"
            )

        if not PDB_ID_PATTERN.fullmatch(pdb_id):
            raise ManifestError(
                f"Row {index}: invalid PDB ID {pdb_id!r}"
            )

        if pdb_id in seen_pdb_ids:
            raise ManifestError(
                f"Row {index}: duplicate PDB ID {pdb_id!r}"
            )

        if s3_key in seen_s3_keys:
            raise ManifestError(
                f"Row {index}: duplicate S3 key {s3_key!r}"
            )

        if source_layout == "canonical_divided_mmcif":
            expected_suffix = (
                f"/{pdb_id[1:3]}/{pdb_id}.cif.gz"
            )

            if not s3_key.startswith(expected_prefix):
                raise ManifestError(
                    f"Row {index}: S3 key is outside snapshot prefix: "
                    f"{s3_key!r}"
                )

            if not s3_key.endswith(expected_suffix):
                raise ManifestError(
                    f"Row {index}: S3 key does not match PDB layout: "
                    f"{s3_key!r}"
                )

        elif source_layout == "recursive_coordinate_files":
            if not s3_key.startswith(f"{expected_snapshot}/"):
                raise ManifestError(
                    f"Row {index}: recursively discovered S3 key is "
                    f"outside snapshot {expected_snapshot}: {s3_key!r}"
                )

            filename = s3_key.rsplit("/", 1)[-1].lower()

            if filename != f"{pdb_id}.cif.gz":
                raise ManifestError(
                    f"Row {index}: recursive coordinate filename does "
                    f"not match PDB ID {pdb_id!r}: {s3_key!r}"
                )

        else:
            raise ManifestError(
                f"Row {index}: unsupported source_layout "
                f"{source_layout!r}"
            )

        if not isinstance(size_bytes, int) or size_bytes <= 0:
            raise ManifestError(
                f"Row {index}: size_bytes must be positive"
            )

        if not etag:
            raise ManifestError(
                f"Row {index}: ETag is empty"
            )

        seen_pdb_ids.add(pdb_id)
        seen_s3_keys.add(s3_key)
        total_bytes += size_bytes

    row_count = len(rows)

    if expected_count is not None and row_count != expected_count:
        raise ManifestError(
            f"Expected {expected_count} rows, found {row_count}"
        )

    if (
        expected_total_bytes is not None
        and total_bytes != expected_total_bytes
    ):
        raise ManifestError(
            "Expected total size "
            f"{expected_total_bytes}, found {total_bytes}"
        )

    return ManifestSummary(
        row_count=row_count,
        total_bytes=total_bytes,
        unique_pdb_ids=len(seen_pdb_ids),
        unique_s3_keys=len(seen_s3_keys),
        sha256="",
    )


def read_manifest_csv(path: str | Path) -> pa.Table:
    """Read a source-manifest CSV with explicit input types."""

    return pacsv.read_csv(
        Path(path),
        convert_options=pacsv.ConvertOptions(
            column_types={
                "snapshot": pa.string(),
                "pdb_id": pa.string(),
                "s3_key": pa.string(),
                "size_bytes": pa.int64(),
                "etag": pa.string(),
            }
        ),
    )


def write_manifest_parquet(
    table: pa.Table,
    output_path: str | Path,
) -> Path:
    """Write a manifest using the explicit Arrow schema."""

    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    missing_timestamp_columns = {
        "last_modified_utc",
        "manifest_generated_at_utc",
    }.difference(table.column_names)

    if missing_timestamp_columns:
        raise ManifestError(
            "Cannot write canonical manifest without columns: "
            + ", ".join(sorted(missing_timestamp_columns))
        )

    canonical = table.select(SOURCE_MANIFEST_SCHEMA.names).cast(
        SOURCE_MANIFEST_SCHEMA
    )

    temporary = output.with_suffix(output.suffix + ".tmp")

    pq.write_table(
        canonical,
        temporary,
        compression="zstd",
        version="2.6",
    )

    temporary.replace(output)
    return output


def manifest_partition_count(
    row_count: int,
    batch_size: int,
) -> int:
    """Return the number of zero-based manifest partitions required."""

    if (
        not isinstance(row_count, int)
        or isinstance(row_count, bool)
        or row_count < 0
    ):
        raise ManifestError(
            "Manifest row_count must be a non-negative integer"
        )

    if (
        not isinstance(batch_size, int)
        or isinstance(batch_size, bool)
        or batch_size <= 0
    ):
        raise ManifestError(
            "Manifest batch_size must be a positive integer"
        )

    if row_count == 0:
        return 0

    return (row_count + batch_size - 1) // batch_size


def select_manifest_partition(
    table: pa.Table,
    *,
    task_id: int,
    batch_size: int,
) -> pa.Table:
    """Select one deterministic zero-based manifest partition."""

    if (
        not isinstance(task_id, int)
        or isinstance(task_id, bool)
        or task_id < 0
    ):
        raise ManifestError(
            "Manifest task_id must be a non-negative integer"
        )

    partition_count = manifest_partition_count(
        table.num_rows,
        batch_size,
    )

    if partition_count == 0:
        raise ManifestError(
            "Cannot select a task from an empty manifest"
        )

    if task_id >= partition_count:
        raise ManifestError(
            f"Manifest task_id {task_id} is out of range; "
            f"valid task IDs are 0 through {partition_count - 1}"
        )

    start = task_id * batch_size
    length = min(
        batch_size,
        table.num_rows - start,
    )

    return table.slice(start, length)
