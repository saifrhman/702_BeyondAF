"""Stage 6: exact bucketing by final cleaned chain length m."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
import json
from pathlib import Path
from typing import Any

import pyarrow as pa
import pyarrow.parquet as pq


class LengthBucketError(RuntimeError):
    """Raised when Stage-6 length bucketing is invalid."""


LENGTH_BUCKET_SCHEMA_NAME = "pdbclean_stage6_length_bucket_index"
LENGTH_BUCKET_SCHEMA_VERSION = "1.0"

LENGTH_BUCKET_GLOBAL_SUMMARY_SCHEMA_NAME = (
    "pdbclean_stage6_length_bucket_global_summary"
)
LENGTH_BUCKET_GLOBAL_SUMMARY_SCHEMA_VERSION = "1.0"

LENGTH_BUCKET_SUCCESS_SCHEMA_NAME = (
    "pdbclean_stage6_length_bucket_success"
)
LENGTH_BUCKET_SUCCESS_SCHEMA_VERSION = "1.0"

BRAIN_UNDEFINED_M1_REASON = "definition_5_1_undefined_for_m1"


LENGTH_BUCKET_SCHEMA = pa.schema(
    [
        pa.field(
            "retained_residue_count",
            pa.int64(),
            nullable=False,
        ),
        pa.field(
            "brain_defined_count",
            pa.int64(),
            nullable=False,
        ),
        pa.field(
            "brain_undefined_count",
            pa.int64(),
            nullable=False,
        ),
        pa.field(
            "total_chain_count",
            pa.int64(),
            nullable=False,
        ),
    ],
    metadata={
        b"schema_name": LENGTH_BUCKET_SCHEMA_NAME.encode(),
        b"schema_version": LENGTH_BUCKET_SCHEMA_VERSION.encode(),
        b"bucket_definition": (
            b"exact final cleaned retained_residue_count m"
        ),
    },
)


@dataclass(frozen=True)
class Stage5BrainPublication:
    root: Path
    brain_path: Path
    undefined_path: Path
    processing_errors_path: Path

    snapshot: str
    cleaning_protocol: str

    input_bri_chain_count: int
    brain_chain_count: int
    undefined_chain_count: int
    processing_error_count: int

    provenance: dict[str, str]


@dataclass(frozen=True)
class LengthBucketPublication:
    bucket_index_path: Path
    global_summary_path: Path
    success_path: Path
    global_summary: dict[str, Any]


def _read_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(
            path.read_text(encoding="utf-8")
        )
    except Exception as exc:
        raise LengthBucketError(
            f"Cannot read JSON artifact {path}: {exc}"
        ) from exc

    if not isinstance(value, dict):
        raise LengthBucketError(
            f"JSON artifact is not an object: {path}"
        )

    return value


def _validated_git_commit(value: str) -> str:
    if (
        not isinstance(value, str)
        or len(value) != 40
        or any(
            character not in "0123456789abcdef"
            for character in value
        )
    ):
        raise LengthBucketError(
            "length_bucket_pipeline_git_commit must "
            "be a 40-character lowercase hexadecimal commit"
        )

    return value


def validate_stage5_brain_publication(
    brain_root: str | Path,
) -> Stage5BrainPublication:
    """Validate the canonical completed Stage-5 publication."""

    root = Path(brain_root)

    success = _read_json(root / "_SUCCESS")
    summary = _read_json(root / "global_summary.json")

    if (
        success.get("success_schema_name")
        != "pdbclean_stage5_brain_success"
        or success.get("success_schema_version") != "1.0"
    ):
        raise LengthBucketError(
            "Unexpected Stage-5 _SUCCESS schema"
        )

    if (
        summary.get("summary_schema_name")
        != "pdbclean_stage5_brain_global_summary"
        or summary.get("summary_schema_version") != "1.0"
    ):
        raise LengthBucketError(
            "Unexpected Stage-5 global-summary schema"
        )

    pointers = {
        "brain_population": "finalized/brain.parquet",
        "undefined_population": "finalized/undefined.parquet",
        "processing_errors": (
            "finalized/processing_errors.parquet"
        ),
    }

    for field, expected in pointers.items():
        if success.get(field) != expected:
            raise LengthBucketError(
                f"Unexpected Stage-5 pointer: {field}"
            )

    provenance_fields = (
        "quality_pipeline_git_commit",
        "geometric_validation_pipeline_git_commit",
        "geometric_validation_finalizer_git_commit",
        "bri_pipeline_git_commit",
        "bri_finalizer_git_commit",
        "brain_pipeline_git_commit",
        "finalizer_pipeline_git_commit",
    )

    provenance = {}

    for field in provenance_fields:
        if success.get(field) != summary.get(field):
            raise LengthBucketError(
                f"Stage-5 provenance mismatch: {field}"
            )

        provenance[field] = success[field]

    for field in (
        "snapshot",
        "cleaning_protocol",
        "task_count",
    ):
        if success.get(field) != summary.get(field):
            raise LengthBucketError(
                f"Stage-5 publication mismatch: {field}"
            )

    brain_path = root / pointers["brain_population"]
    undefined_path = root / pointers["undefined_population"]
    error_path = root / pointers["processing_errors"]

    for path in (
        brain_path,
        undefined_path,
        error_path,
    ):
        if not path.is_file():
            raise LengthBucketError(
                f"Missing Stage-5 finalized artifact: {path}"
            )

    brain_count = pq.read_metadata(brain_path).num_rows
    undefined_count = pq.read_metadata(undefined_path).num_rows
    error_count = pq.read_metadata(error_path).num_rows

    if brain_count != summary.get("brain_chain_count"):
        raise LengthBucketError(
            "Stage-5 Brain row-count mismatch"
        )

    if undefined_count != summary.get(
        "undefined_chain_count"
    ):
        raise LengthBucketError(
            "Stage-5 undefined row-count mismatch"
        )

    if error_count != summary.get(
        "processing_error_count"
    ):
        raise LengthBucketError(
            "Stage-5 processing-error row-count mismatch"
        )

    if error_count != 0:
        raise LengthBucketError(
            "Stage 6 requires zero Stage-5 processing errors"
        )

    input_count = summary.get(
        "input_bri_chain_count"
    )

    if input_count != brain_count + undefined_count:
        raise LengthBucketError(
            "Stage-5 population accounting mismatch"
        )

    return Stage5BrainPublication(
        root=root,
        brain_path=brain_path,
        undefined_path=undefined_path,
        processing_errors_path=error_path,
        snapshot=summary["snapshot"],
        cleaning_protocol=summary["cleaning_protocol"],
        input_bri_chain_count=input_count,
        brain_chain_count=brain_count,
        undefined_chain_count=undefined_count,
        processing_error_count=error_count,
        provenance=provenance,
    )


def build_length_bucket_index(
    upstream: Stage5BrainPublication,
) -> pa.Table:
    """Count every finalized Stage-5 chain in its exact-m bucket."""

    defined = Counter()
    undefined = Counter()

    brain_pf = pq.ParquetFile(
        upstream.brain_path
    )

    for batch in brain_pf.iter_batches(
        columns=["retained_residue_count"],
        batch_size=65536,
    ):
        for m in batch.column(0).to_pylist():
            if (
                not isinstance(m, int)
                or isinstance(m, bool)
                or m < 2
            ):
                raise LengthBucketError(
                    "Brain-defined population contains invalid m"
                )

            defined[m] += 1

    undefined_pf = pq.ParquetFile(
        upstream.undefined_path
    )

    for batch in undefined_pf.iter_batches(
        columns=[
            "retained_residue_count",
            "undefined_reason",
        ],
        batch_size=65536,
    ):
        for row in batch.to_pylist():
            m = row["retained_residue_count"]

            if m != 1:
                raise LengthBucketError(
                    "Brain-undefined population contains m != 1"
                )

            if (
                row["undefined_reason"]
                != BRAIN_UNDEFINED_M1_REASON
            ):
                raise LengthBucketError(
                    "Unexpected Stage-5 undefined reason"
                )

            undefined[m] += 1

    lengths = sorted(
        set(defined) | set(undefined)
    )

    rows = []

    for m in lengths:
        brain_count = defined[m]
        undefined_count = undefined[m]

        rows.append(
            {
                "retained_residue_count": m,
                "brain_defined_count": brain_count,
                "brain_undefined_count": undefined_count,
                "total_chain_count": (
                    brain_count + undefined_count
                ),
            }
        )

    table = pa.Table.from_pylist(
        rows,
        schema=LENGTH_BUCKET_SCHEMA,
    )

    if (
        sum(
            table["total_chain_count"].to_pylist()
        )
        != upstream.input_bri_chain_count
    ):
        raise LengthBucketError(
            "Stage-6 population accounting failed"
        )

    if (
        sum(
            table["brain_defined_count"].to_pylist()
        )
        != upstream.brain_chain_count
    ):
        raise LengthBucketError(
            "Stage-6 Brain-defined accounting failed"
        )

    if (
        sum(
            table["brain_undefined_count"].to_pylist()
        )
        != upstream.undefined_chain_count
    ):
        raise LengthBucketError(
            "Stage-6 undefined accounting failed"
        )

    return table


def _write_parquet_atomic(
    table: pa.Table,
    path: Path,
) -> Path:
    path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    temporary = path.with_suffix(
        path.suffix + ".tmp"
    )

    try:
        pq.write_table(
            table,
            temporary,
            compression="zstd",
            version="2.6",
        )
        temporary.replace(path)
    except Exception:
        if temporary.exists():
            temporary.unlink()
        raise

    return path


def _write_json_atomic(
    payload: dict[str, Any],
    path: Path,
) -> Path:
    path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    temporary = path.with_suffix(
        path.suffix + ".tmp"
    )

    temporary.write_text(
        json.dumps(
            payload,
            sort_keys=True,
            indent=2,
            ensure_ascii=True,
        )
        + "\n",
        encoding="utf-8",
    )

    temporary.replace(path)

    return path


def finalize_length_buckets(
    *,
    upstream: Stage5BrainPublication,
    output_root: str | Path,
    length_bucket_pipeline_git_commit: str,
) -> LengthBucketPublication:
    """Publish the canonical Stage-6 exact-length bucket index."""

    commit = _validated_git_commit(
        length_bucket_pipeline_git_commit
    )

    root = Path(output_root)
    success_path = root / "_SUCCESS"

    if success_path.exists():
        success_path.unlink()

    table = build_length_bucket_index(
        upstream
    )

    if table.num_rows == 0:
        raise LengthBucketError(
            "Stage-6 bucket index cannot be empty"
        )

    m_values = table[
        "retained_residue_count"
    ].to_pylist()

    bucket_path = _write_parquet_atomic(
        table,
        root / "finalized/bucket_index.parquet",
    )

    summary = {
        "summary_schema_name": (
            LENGTH_BUCKET_GLOBAL_SUMMARY_SCHEMA_NAME
        ),
        "summary_schema_version": (
            LENGTH_BUCKET_GLOBAL_SUMMARY_SCHEMA_VERSION
        ),
        "snapshot": upstream.snapshot,
        "cleaning_protocol": (
            upstream.cleaning_protocol
        ),
        **upstream.provenance,
        "length_bucket_pipeline_git_commit": commit,
        "input_chain_count": (
            upstream.input_bri_chain_count
        ),
        "brain_defined_chain_count": (
            upstream.brain_chain_count
        ),
        "brain_undefined_chain_count": (
            upstream.undefined_chain_count
        ),
        "distinct_length_bucket_count": (
            table.num_rows
        ),
        "minimum_retained_residue_count": (
            min(m_values)
        ),
        "maximum_retained_residue_count": (
            max(m_values)
        ),
        "m1_chain_count": int(
            table.filter(
                pa.compute.equal(
                    table["retained_residue_count"],
                    1,
                )
            )["total_chain_count"][0].as_py()
        )
        if 1 in m_values
        else 0,
        "population_accounting_valid": True,
    }

    summary_path = _write_json_atomic(
        summary,
        root / "global_summary.json",
    )

    success = {
        "success_schema_name": (
            LENGTH_BUCKET_SUCCESS_SCHEMA_NAME
        ),
        "success_schema_version": (
            LENGTH_BUCKET_SUCCESS_SCHEMA_VERSION
        ),
        "snapshot": upstream.snapshot,
        "cleaning_protocol": (
            upstream.cleaning_protocol
        ),
        **upstream.provenance,
        "length_bucket_pipeline_git_commit": commit,
        "global_summary": "global_summary.json",
        "bucket_index": (
            "finalized/bucket_index.parquet"
        ),
    }

    # Completion marker strictly last.
    success_path = _write_json_atomic(
        success,
        success_path,
    )

    return LengthBucketPublication(
        bucket_index_path=bucket_path,
        global_summary_path=summary_path,
        success_path=success_path,
        global_summary=summary,
    )
