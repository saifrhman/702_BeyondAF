"""Stage-5 Brain production orchestration and upstream validation."""

from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
from typing import Any

import pyarrow as pa
import pyarrow.parquet as pq

from pdbclean.bri_finalize import (
    BRI_GLOBAL_SUMMARY_SCHEMA_NAME,
    BRI_GLOBAL_SUMMARY_SCHEMA_VERSION,
    BRI_SUCCESS_SCHEMA_NAME,
    BRI_SUCCESS_SCHEMA_VERSION,
)
from pdbclean.schemas import STAGE3_BRI_CHAIN_SCHEMA


class BrainRunnerError(RuntimeError):
    """Raised when Stage-5 Brain production cannot proceed safely."""


@dataclass(frozen=True)
class UpstreamBRI:
    """Validated completed Stage-3 publication consumed by Stage 5."""

    bri_root: Path
    bri_path: Path

    snapshot: str
    cleaning_protocol: str
    task_count: int
    bri_chain_count: int
    minimum_retained_residue_count: int
    maximum_retained_residue_count: int

    quality_pipeline_git_commit: str
    geometric_validation_pipeline_git_commit: str
    geometric_validation_finalizer_git_commit: str
    bri_pipeline_git_commit: str
    bri_finalizer_git_commit: str


def _read_json_object(
    path: Path,
    *,
    description: str,
) -> dict[str, Any]:
    if not path.is_file():
        raise BrainRunnerError(
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
        raise BrainRunnerError(
            f"Cannot read {description} {path}: {exc}"
        ) from exc

    if not isinstance(payload, dict):
        raise BrainRunnerError(
            f"{description} must contain a JSON object: {path}"
        )

    return payload


def _validated_git_commit(
    value: object,
    *,
    field: str,
) -> str:
    if (
        not isinstance(value, str)
        or len(value) != 40
        or any(
            character not in "0123456789abcdef"
            for character in value.lower()
        )
    ):
        raise BrainRunnerError(
            f"Invalid Stage-3 Git commit field: {field}"
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
        raise BrainRunnerError(
            f"Stage-3 field {field} must be a non-negative integer"
        )

    return value


def _parquet_child_field_matches(
    observed: pa.Field,
    expected: pa.Field,
) -> bool:
    if (
        observed.nullable != expected.nullable
        or observed.metadata != expected.metadata
    ):
        return False

    return _parquet_type_matches(
        observed.type,
        expected.type,
    )


def _parquet_type_matches(
    observed: pa.DataType,
    expected: pa.DataType,
) -> bool:
    if (
        pa.types.is_list(observed)
        and pa.types.is_list(expected)
    ):
        return _parquet_child_field_matches(
            observed.value_field,
            expected.value_field,
        )

    if (
        pa.types.is_fixed_size_list(observed)
        and pa.types.is_fixed_size_list(expected)
    ):
        return (
            observed.list_size == expected.list_size
            and _parquet_child_field_matches(
                observed.value_field,
                expected.value_field,
            )
        )

    return observed == expected


def _parquet_schema_matches(
    observed: pa.Schema,
    expected: pa.Schema,
) -> bool:
    if observed.metadata != expected.metadata:
        return False

    if len(observed) != len(expected):
        return False

    for observed_field, expected_field in zip(
        observed,
        expected,
        strict=True,
    ):
        if observed_field.name != expected_field.name:
            return False

        if (
            observed_field.nullable != expected_field.nullable
            or observed_field.metadata != expected_field.metadata
            or not _parquet_type_matches(
                observed_field.type,
                expected_field.type,
            )
        ):
            return False

    return True


def validate_upstream_bri_stage(
    bri_root: str | Path,
    *,
    expected_snapshot: str,
    expected_cleaning_protocol: str,
    expected_task_count: int | None = None,
) -> UpstreamBRI:
    """Validate the completed canonical Stage-3 input to Stage 5."""

    root = Path(bri_root)

    success_path = root / "_SUCCESS"
    summary_path = root / "global_summary.json"

    success = _read_json_object(
        success_path,
        description="Stage-3 BRI _SUCCESS marker",
    )
    summary = _read_json_object(
        summary_path,
        description="Stage-3 BRI global summary",
    )

    if (
        success.get("success_schema_name")
        != BRI_SUCCESS_SCHEMA_NAME
        or success.get("success_schema_version")
        != BRI_SUCCESS_SCHEMA_VERSION
    ):
        raise BrainRunnerError(
            "Unexpected Stage-3 BRI _SUCCESS schema"
        )

    if (
        summary.get("summary_schema_name")
        != BRI_GLOBAL_SUMMARY_SCHEMA_NAME
        or summary.get("summary_schema_version")
        != BRI_GLOBAL_SUMMARY_SCHEMA_VERSION
    ):
        raise BrainRunnerError(
            "Unexpected Stage-3 BRI global-summary schema"
        )

    canonical_pointers = {
        "global_summary": "global_summary.json",
        "finalized_directory": "finalized",
        "bri_population": "finalized/bri.parquet",
    }

    for field, expected in canonical_pointers.items():
        if success.get(field) != expected:
            raise BrainRunnerError(
                f"Unexpected Stage-3 BRI completion pointer: {field}"
            )

    for field, expected in (
        ("snapshot", expected_snapshot),
        ("cleaning_protocol", expected_cleaning_protocol),
    ):
        if success.get(field) != expected:
            raise BrainRunnerError(
                f"Stage-3 BRI _SUCCESS {field} mismatch"
            )

        if summary.get(field) != expected:
            raise BrainRunnerError(
                f"Stage-3 BRI global-summary {field} mismatch"
            )

    provenance_fields = (
        "quality_pipeline_git_commit",
        "geometric_validation_pipeline_git_commit",
        "geometric_validation_finalizer_git_commit",
        "bri_pipeline_git_commit",
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
            raise BrainRunnerError(
                f"Stage-3 BRI provenance mismatch: {field}"
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
        raise BrainRunnerError(
            "Stage-3 BRI task-count mismatch"
        )

    if (
        expected_task_count is not None
        and success_task_count != expected_task_count
    ):
        raise BrainRunnerError(
            "Stage-3 BRI task count does not match the "
            "current manifest partition contract"
        )

    input_count = _nonnegative_integer(
        summary,
        "input_eligible_chain_count",
    )
    bri_count = _nonnegative_integer(
        summary,
        "bri_chain_count",
    )
    processing_error_count = _nonnegative_integer(
        summary,
        "processing_error_count",
    )
    unique_eligible_count = _nonnegative_integer(
        summary,
        "unique_eligible_identity_count",
    )
    unique_bri_count = _nonnegative_integer(
        summary,
        "unique_bri_identity_count",
    )

    manifest_source_count = _nonnegative_integer(
        summary,
        "manifest_source_object_count",
    )
    relevant_source_count = _nonnegative_integer(
        summary,
        "relevant_source_object_count",
    )
    downloaded_source_count = _nonnegative_integer(
        summary,
        "downloaded_source_object_count",
    )
    parsed_source_count = _nonnegative_integer(
        summary,
        "parsed_source_object_count",
    )

    minimum_m = _nonnegative_integer(
        summary,
        "minimum_retained_residue_count",
    )
    maximum_m = _nonnegative_integer(
        summary,
        "maximum_retained_residue_count",
    )

    if summary.get("chain_accounting_valid") is not True:
        raise BrainRunnerError(
            "Stage-3 BRI global chain accounting is not valid"
        )

    if processing_error_count != 0:
        raise BrainRunnerError(
            "Stage-3 BRI completed population contains "
            "processing errors"
        )

    if input_count != bri_count + processing_error_count:
        raise BrainRunnerError(
            "Stage-3 BRI global chain counts do not reconcile"
        )

    if (
        unique_eligible_count != input_count
        or unique_bri_count != bri_count
    ):
        raise BrainRunnerError(
            "Stage-3 BRI global identity counts do not reconcile"
        )

    if bri_count < 1:
        raise BrainRunnerError(
            "Stage-3 BRI completed population is empty"
        )

    if (
        minimum_m < 1
        or maximum_m < minimum_m
    ):
        raise BrainRunnerError(
            "Stage-3 BRI retained-residue range is invalid"
        )

    if relevant_source_count > manifest_source_count:
        raise BrainRunnerError(
            "Stage-3 BRI relevant source count exceeds manifest count"
        )

    if not (
        relevant_source_count
        == downloaded_source_count
        == parsed_source_count
    ):
        raise BrainRunnerError(
            "Stage-3 BRI source accounting does not reconcile"
        )

    bri_path = (
        root
        / "finalized"
        / "bri.parquet"
    )

    if not bri_path.is_file():
        raise BrainRunnerError(
            "Canonical Stage-3 BRI population does not exist: "
            f"{bri_path}"
        )

    try:
        observed_schema = pq.read_schema(
            bri_path
        )
    except Exception as exc:
        raise BrainRunnerError(
            f"Cannot read canonical Stage-3 BRI schema: {exc}"
        ) from exc

    if not _parquet_schema_matches(
        observed_schema,
        STAGE3_BRI_CHAIN_SCHEMA,
    ):
        raise BrainRunnerError(
            "Canonical Stage-3 BRI population schema mismatch"
        )

    try:
        parquet_row_count = (
            pq.ParquetFile(
                bri_path
            ).metadata.num_rows
        )
    except Exception as exc:
        raise BrainRunnerError(
            f"Cannot read canonical Stage-3 BRI metadata: {exc}"
        ) from exc

    if parquet_row_count != bri_count:
        raise BrainRunnerError(
            "Canonical Stage-3 BRI row count does not match "
            "Stage-3 global summary"
        )

    return UpstreamBRI(
        bri_root=root,
        bri_path=bri_path,
        snapshot=expected_snapshot,
        cleaning_protocol=expected_cleaning_protocol,
        task_count=success_task_count,
        bri_chain_count=bri_count,
        minimum_retained_residue_count=minimum_m,
        maximum_retained_residue_count=maximum_m,
        quality_pipeline_git_commit=provenance[
            "quality_pipeline_git_commit"
        ],
        geometric_validation_pipeline_git_commit=provenance[
            "geometric_validation_pipeline_git_commit"
        ],
        geometric_validation_finalizer_git_commit=provenance[
            "geometric_validation_finalizer_git_commit"
        ],
        bri_pipeline_git_commit=provenance[
            "bri_pipeline_git_commit"
        ],
        bri_finalizer_git_commit=provenance[
            "finalizer_pipeline_git_commit"
        ],
    )
