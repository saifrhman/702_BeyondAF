"""Production orchestration for Stage-3 Definition 3.4 BRI.

This module consumes only the canonical population published by the
completed post-cleaning geometric-validation stage.

The scientific BRI calculation itself lives in ``pdbclean.bri``.
"""

from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
from typing import Any

import pyarrow as pa
import pyarrow.parquet as pq

from pdbclean.geometric_validation_finalize import (
    GEOMETRIC_VALIDATION_GLOBAL_SUMMARY_SCHEMA_NAME,
    GEOMETRIC_VALIDATION_GLOBAL_SUMMARY_SCHEMA_VERSION,
    GEOMETRIC_VALIDATION_SUCCESS_SCHEMA_NAME,
    GEOMETRIC_VALIDATION_SUCCESS_SCHEMA_VERSION,
)
from pdbclean.schemas import STAGE3_ELIGIBLE_CHAIN_SCHEMA


class BRIRunnerError(RuntimeError):
    """Raised when Stage-3 BRI orchestration cannot proceed safely."""


@dataclass(frozen=True)
class UpstreamGeometricValidation:
    """Validated completion state consumed by Stage-3 BRI."""

    geometric_validation_root: Path
    eligible_path: Path

    snapshot: str
    cleaning_protocol: str
    task_count: int
    eligible_chain_count: int

    quality_pipeline_git_commit: str
    geometric_validation_pipeline_git_commit: str
    geometric_validation_finalizer_git_commit: str


def _read_json_object(
    path: Path,
    *,
    description: str,
) -> dict[str, Any]:
    """Read one required JSON object."""

    if not path.is_file():
        raise BRIRunnerError(
            f"Required {description} does not exist: {path}"
        )

    try:
        payload = json.loads(
            path.read_text(encoding="utf-8")
        )
    except (
        OSError,
        UnicodeError,
        json.JSONDecodeError,
    ) as exc:
        raise BRIRunnerError(
            f"Cannot read {description} {path}: {exc}"
        ) from exc

    if not isinstance(payload, dict):
        raise BRIRunnerError(
            f"{description} must contain a JSON object: {path}"
        )

    return payload


def _validated_git_commit(
    value: object,
    *,
    field: str,
) -> str:
    """Validate one exact Git commit identifier."""

    if (
        not isinstance(value, str)
        or len(value) != 40
        or any(
            character not in "0123456789abcdef"
            for character in value.lower()
        )
    ):
        raise BRIRunnerError(
            f"Invalid upstream Git commit field: {field}"
        )

    return value


def _nonnegative_integer(
    payload: dict[str, Any],
    field: str,
) -> int:
    value = payload.get(field)

    if (
        not isinstance(value, int)
        or isinstance(value, bool)
        or value < 0
    ):
        raise BRIRunnerError(
            f"Upstream field {field} must be a non-negative integer"
        )

    return value


def _parquet_field_matches(
    observed: pa.Field,
    expected: pa.Field,
) -> bool:
    """Compare fields while tolerating Parquet list-child renaming."""

    if (
        observed.name != expected.name
        or observed.nullable != expected.nullable
        or observed.metadata != expected.metadata
    ):
        return False

    observed_type = observed.type
    expected_type = expected.type

    if (
        pa.types.is_list(observed_type)
        and pa.types.is_list(expected_type)
    ):
        observed_value = observed_type.value_field
        expected_value = expected_type.value_field

        return (
            observed_value.type == expected_value.type
            and observed_value.nullable == expected_value.nullable
            and observed_value.metadata == expected_value.metadata
        )

    return observed_type == expected_type


def _parquet_schema_matches(
    observed: pa.Schema,
    expected: pa.Schema,
) -> bool:
    """Compare a Parquet schema to its canonical Arrow contract."""

    if observed.metadata != expected.metadata:
        return False

    if len(observed) != len(expected):
        return False

    return all(
        _parquet_field_matches(observed_field, expected_field)
        for observed_field, expected_field in zip(
            observed,
            expected,
            strict=True,
        )
    )



def validate_upstream_geometric_validation_stage(
    geometric_validation_root: str | Path,
    *,
    expected_snapshot: str,
    expected_cleaning_protocol: str,
    expected_task_count: int | None = None,
) -> UpstreamGeometricValidation:
    """Validate the canonical completed Stage-2 input to Stage 3.

    Stage 3 trusts only a completed geometric-validation publication.
    The completion marker, global summary, canonical eligible Parquet
    population, accounting, and producer provenance must agree.

    No scientific recomputation occurs here.
    """

    root = Path(geometric_validation_root)

    success_path = root / "_SUCCESS"
    summary_path = root / "global_summary.json"

    success = _read_json_object(
        success_path,
        description="Stage-2 _SUCCESS marker",
    )
    summary = _read_json_object(
        summary_path,
        description="Stage-2 global summary",
    )

    if (
        success.get("success_schema_name")
        != GEOMETRIC_VALIDATION_SUCCESS_SCHEMA_NAME
        or success.get("success_schema_version")
        != GEOMETRIC_VALIDATION_SUCCESS_SCHEMA_VERSION
    ):
        raise BRIRunnerError(
            "Unexpected Stage-2 _SUCCESS schema"
        )

    if (
        summary.get("summary_schema_name")
        != GEOMETRIC_VALIDATION_GLOBAL_SUMMARY_SCHEMA_NAME
        or summary.get("summary_schema_version")
        != GEOMETRIC_VALIDATION_GLOBAL_SUMMARY_SCHEMA_VERSION
    ):
        raise BRIRunnerError(
            "Unexpected Stage-2 global-summary schema"
        )

    canonical_pointers = {
        "global_summary": "global_summary.json",
        "finalized_directory": "finalized",
        "eligible_population": "finalized/eligible.parquet",
        "quarantined_population": "finalized/quarantined.parquet",
    }

    for field, expected in canonical_pointers.items():
        if success.get(field) != expected:
            raise BRIRunnerError(
                f"Unexpected Stage-2 completion pointer: {field}"
            )

    for field, expected in (
        ("snapshot", expected_snapshot),
        ("cleaning_protocol", expected_cleaning_protocol),
    ):
        if success.get(field) != expected:
            raise BRIRunnerError(
                f"Stage-2 _SUCCESS {field} mismatch"
            )

        if summary.get(field) != expected:
            raise BRIRunnerError(
                f"Stage-2 global-summary {field} mismatch"
            )

    provenance_fields = (
        "quality_pipeline_git_commit",
        "geometric_validation_pipeline_git_commit",
        "finalizer_pipeline_git_commit",
    )

    provenance: dict[str, str] = {}

    for field in provenance_fields:
        success_commit = _validated_git_commit(
            success.get(field),
            field=f"_SUCCESS.{field}",
        )
        summary_commit = _validated_git_commit(
            summary.get(field),
            field=f"global_summary.{field}",
        )

        if success_commit != summary_commit:
            raise BRIRunnerError(
                f"Stage-2 provenance mismatch: {field}"
            )

        provenance[field] = success_commit

    success_task_count = _nonnegative_integer(
        success,
        "task_count",
    )
    summary_task_count = _nonnegative_integer(
        summary,
        "task_count",
    )

    if success_task_count != summary_task_count:
        raise BRIRunnerError(
            "Stage-2 task-count mismatch"
        )

    if (
        expected_task_count is not None
        and success_task_count != expected_task_count
    ):
        raise BRIRunnerError(
            "Stage-2 task count does not match the current "
            "manifest partition contract"
        )

    input_count = _nonnegative_integer(
        summary,
        "input_accepted_chain_count",
    )
    eligible_count = _nonnegative_integer(
        summary,
        "eligible_chain_count",
    )
    quarantined_count = _nonnegative_integer(
        summary,
        "quarantined_chain_count",
    )
    processing_error_count = _nonnegative_integer(
        summary,
        "processing_error_count",
    )

    if summary.get("chain_accounting_valid") is not True:
        raise BRIRunnerError(
            "Stage-2 global chain accounting is not valid"
        )

    if (
        input_count
        != eligible_count
        + quarantined_count
        + processing_error_count
    ):
        raise BRIRunnerError(
            "Stage-2 global chain counts do not reconcile"
        )

    if processing_error_count != 0:
        raise BRIRunnerError(
            "Stage-2 completed population contains processing errors"
        )

    eligible_path = (
        root / "finalized" / "eligible.parquet"
    )

    if not eligible_path.is_file():
        raise BRIRunnerError(
            "Canonical Stage-3 eligible population does not exist: "
            f"{eligible_path}"
        )

    observed_schema = pq.read_schema(eligible_path)

    if not _parquet_schema_matches(
        observed_schema,
        STAGE3_ELIGIBLE_CHAIN_SCHEMA,
    ):
        raise BRIRunnerError(
            "Canonical Stage-3 eligible population schema mismatch"
        )

    parquet = pq.ParquetFile(eligible_path)
    parquet_row_count = parquet.metadata.num_rows

    if parquet_row_count != eligible_count:
        raise BRIRunnerError(
            "Canonical Stage-3 eligible row count does not match "
            "Stage-2 global summary"
        )

    return UpstreamGeometricValidation(
        geometric_validation_root=root,
        eligible_path=eligible_path,
        snapshot=expected_snapshot,
        cleaning_protocol=expected_cleaning_protocol,
        task_count=success_task_count,
        eligible_chain_count=eligible_count,
        quality_pipeline_git_commit=provenance[
            "quality_pipeline_git_commit"
        ],
        geometric_validation_pipeline_git_commit=provenance[
            "geometric_validation_pipeline_git_commit"
        ],
        geometric_validation_finalizer_git_commit=provenance[
            "finalizer_pipeline_git_commit"
        ],
    )
