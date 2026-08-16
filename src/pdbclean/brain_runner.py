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


# ---------------------------------------------------------------------------
# Stage-5 single-record Brain processing
# ---------------------------------------------------------------------------

from typing import Callable, Mapping

import numpy as np

from pdbclean.brain import compute_brain
from pdbclean.schemas import (
    STAGE5_BRAIN_CHAIN_SCHEMA,
    STAGE5_BRAIN_PROCESSING_ERROR_SCHEMA,
    STAGE5_BRAIN_UNDEFINED_CHAIN_SCHEMA,
)


BRAIN_UNDEFINED_M1_REASON = (
    "definition_5_1_undefined_for_m1"
)


@dataclass(frozen=True)
class BrainRecordResult:
    """Terminal Stage-5 outcome for one canonical Stage-3 BRI row."""

    input_bri_chain_count: int = 1

    brain_records: tuple[dict[str, Any], ...] = ()
    undefined_records: tuple[dict[str, Any], ...] = ()
    processing_errors: tuple[dict[str, Any], ...] = ()

    @property
    def chain_accounting_valid(self) -> bool:
        """Exactly one terminal outcome must exist for one input row."""

        return self.input_bri_chain_count == (
            len(self.brain_records)
            + len(self.undefined_records)
            + len(self.processing_errors)
        )


def _brain_row_identity(
    row: Mapping[str, Any],
) -> tuple[str, str, int, str]:
    """Return the canonical Stage-3 BRI chain identity."""

    try:
        return (
            row["snapshot"],
            row["pdb_id"],
            row["model_id"],
            row["label_chain_id"],
        )
    except KeyError as exc:
        raise BrainRunnerError(
            "Canonical Stage-3 BRI row is missing identity lineage"
        ) from exc


def _brain_input_lineage_issue(
    row: Mapping[str, Any],
    *,
    upstream: UpstreamBRI,
) -> str | None:
    """Return an immutable Stage-3 lineage mismatch, if any."""

    missing = [
        field.name
        for field in STAGE3_BRI_CHAIN_SCHEMA
        if field.name not in row
    ]

    if missing:
        raise BrainRunnerError(
            "Canonical Stage-3 BRI row missing required field(s): "
            + ", ".join(missing)
        )

    if row["snapshot"] != upstream.snapshot:
        return (
            "BRI row snapshot does not match completed "
            "Stage-3 publication"
        )

    if row["cleaning_protocol"] != upstream.cleaning_protocol:
        return (
            "BRI row cleaning_protocol does not match completed "
            "Stage-3 publication"
        )

    provenance_checks = (
        (
            "quality_pipeline_git_commit",
            upstream.quality_pipeline_git_commit,
        ),
        (
            "geometric_validation_pipeline_git_commit",
            upstream.geometric_validation_pipeline_git_commit,
        ),
        (
            "geometric_validation_finalizer_git_commit",
            upstream.geometric_validation_finalizer_git_commit,
        ),
        (
            "bri_pipeline_git_commit",
            upstream.bri_pipeline_git_commit,
        ),
    )

    for field, expected in provenance_checks:
        if row[field] != expected:
            return (
                f"BRI row {field} does not match completed "
                "Stage-3 publication"
            )

    model_id = row["model_id"]

    if (
        not isinstance(model_id, int)
        or isinstance(model_id, bool)
        or model_id <= 0
    ):
        return "BRI row model_id must be a positive integer"

    label_chain_id = row["label_chain_id"]

    if (
        not isinstance(label_chain_id, str)
        or not label_chain_id
    ):
        return "BRI row label_chain_id must be non-empty"

    retained_count = row["retained_residue_count"]
    retained_ids = row["retained_label_seq_ids"]

    if (
        not isinstance(retained_count, int)
        or isinstance(retained_count, bool)
        or retained_count < 1
    ):
        return (
            "BRI row retained_residue_count must be "
            "a positive integer"
        )

    if not isinstance(retained_ids, (list, tuple)):
        return (
            "BRI row retained_label_seq_ids must be a sequence"
        )

    if len(retained_ids) != retained_count:
        return (
            "BRI row retained_residue_count does not match "
            "retained_label_seq_ids"
        )

    return None


def _validated_stage3_bri_matrix(
    value: object,
    *,
    expected_residue_count: int,
) -> np.ndarray:
    """Validate one canonical Stage-3 BRI matrix for Brain input."""

    try:
        matrix = np.asarray(
            value,
            dtype=np.float64,
        )
    except (TypeError, ValueError) as exc:
        raise ValueError(
            "Stage-3 BRI payload cannot be converted to float64"
        ) from exc

    if matrix.shape != (
        expected_residue_count,
        9,
    ):
        raise ValueError(
            "Stage-3 BRI payload has an unexpected matrix shape"
        )

    if not np.isfinite(matrix).all():
        raise ValueError(
            "Stage-3 BRI payload contains non-finite coordinates"
        )

    if not np.array_equal(
        matrix,
        np.around(matrix, 3),
    ):
        raise ValueError(
            "Stage-3 BRI payload is not canonical at 3 decimals"
        )

    return matrix


def _validated_brain_vector(
    value: object,
) -> np.ndarray:
    """Validate the scientific output of compute_brain()."""

    if not isinstance(value, np.ndarray):
        raise ValueError(
            "compute_brain() did not return a NumPy array"
        )

    if value.dtype != np.float64:
        raise ValueError(
            "compute_brain() did not return float64 coordinates"
        )

    if value.shape != (9,):
        raise ValueError(
            "compute_brain() returned an unexpected vector shape"
        )

    if not np.isfinite(value).all():
        raise ValueError(
            "compute_brain() returned non-finite coordinates"
        )

    return value


def _stage5_upstream_record_fields(
    row: Mapping[str, Any],
) -> dict[str, Any]:
    """Copy exact Stage-3 lineage without duplicating the BRI matrix."""

    return {
        field.name: row[field.name]
        for field in STAGE3_BRI_CHAIN_SCHEMA
        if field.name != "bri"
    }


def _brain_defined_record(
    row: Mapping[str, Any],
    *,
    brain: np.ndarray,
    upstream: UpstreamBRI,
    brain_pipeline_git_commit: str,
) -> dict[str, Any]:
    """Materialize one Definition-5.1 Brain-defined record."""

    record = _stage5_upstream_record_fields(
        row
    )

    record.update(
        {
            "bri_finalizer_git_commit": (
                upstream.bri_finalizer_git_commit
            ),
            "brain_pipeline_git_commit": (
                brain_pipeline_git_commit
            ),
            "brain": brain.tolist(),
        }
    )

    return record


def _brain_undefined_record(
    row: Mapping[str, Any],
    *,
    upstream: UpstreamBRI,
    brain_pipeline_git_commit: str,
) -> dict[str, Any]:
    """Materialize the legitimate m=1 Brain-bypass outcome."""

    record = _stage5_upstream_record_fields(
        row
    )

    record.update(
        {
            "bri_finalizer_git_commit": (
                upstream.bri_finalizer_git_commit
            ),
            "brain_pipeline_git_commit": (
                brain_pipeline_git_commit
            ),
            "undefined_reason": (
                BRAIN_UNDEFINED_M1_REASON
            ),
        }
    )

    return record


def _brain_processing_error(
    row: Mapping[str, Any],
    *,
    processing_stage: str,
    error_type: str,
    error_message: str,
    upstream: UpstreamBRI,
    brain_pipeline_git_commit: str,
) -> dict[str, Any]:
    """Materialize one genuine Stage-5 processing failure."""

    record = {
        field.name: row[field.name]
        for field in STAGE5_BRAIN_PROCESSING_ERROR_SCHEMA
        if field.name
        not in {
            "bri_finalizer_git_commit",
            "brain_pipeline_git_commit",
            "processing_stage",
            "error_type",
            "error_message",
        }
    }

    record.update(
        {
            "bri_finalizer_git_commit": (
                upstream.bri_finalizer_git_commit
            ),
            "brain_pipeline_git_commit": (
                brain_pipeline_git_commit
            ),
            "processing_stage": processing_stage,
            "error_type": error_type,
            "error_message": error_message,
        }
    )

    return record


def process_brain_record(
    row: Mapping[str, Any],
    *,
    upstream: UpstreamBRI,
    brain_pipeline_git_commit: str,
    brain_computer: Callable[
        [np.ndarray],
        np.ndarray,
    ] = compute_brain,
) -> BrainRecordResult:
    """Produce exactly one terminal Stage-5 outcome for one BRI row."""

    brain_pipeline_git_commit = _validated_git_commit(
        brain_pipeline_git_commit,
        field="brain_pipeline_git_commit",
    )

    _brain_row_identity(row)

    lineage_issue = _brain_input_lineage_issue(
        row,
        upstream=upstream,
    )

    if lineage_issue is not None:
        result = BrainRecordResult(
            processing_errors=(
                _brain_processing_error(
                    row,
                    processing_stage=(
                        "bri_lineage_validation"
                    ),
                    error_type="BRILineageError",
                    error_message=lineage_issue,
                    upstream=upstream,
                    brain_pipeline_git_commit=(
                        brain_pipeline_git_commit
                    ),
                ),
            ),
        )

        if not result.chain_accounting_valid:
            raise BrainRunnerError(
                "Stage-5 single-record accounting failed"
            )

        return result

    retained_count = row[
        "retained_residue_count"
    ]

    try:
        bri = _validated_stage3_bri_matrix(
            row["bri"],
            expected_residue_count=retained_count,
        )
    except ValueError as exc:
        result = BrainRecordResult(
            processing_errors=(
                _brain_processing_error(
                    row,
                    processing_stage="bri_input_validation",
                    error_type=type(exc).__name__,
                    error_message=str(exc),
                    upstream=upstream,
                    brain_pipeline_git_commit=(
                        brain_pipeline_git_commit
                    ),
                ),
            ),
        )

        if not result.chain_accounting_valid:
            raise BrainRunnerError(
                "Stage-5 single-record accounting failed"
            )

        return result

    if retained_count == 1:
        result = BrainRecordResult(
            undefined_records=(
                _brain_undefined_record(
                    row,
                    upstream=upstream,
                    brain_pipeline_git_commit=(
                        brain_pipeline_git_commit
                    ),
                ),
            ),
        )

        if not result.chain_accounting_valid:
            raise BrainRunnerError(
                "Stage-5 single-record accounting failed"
            )

        return result

    try:
        brain = brain_computer(
            bri
        )
        brain = _validated_brain_vector(
            brain
        )
    except ValueError as exc:
        result = BrainRecordResult(
            processing_errors=(
                _brain_processing_error(
                    row,
                    processing_stage="brain_computation",
                    error_type=type(exc).__name__,
                    error_message=str(exc),
                    upstream=upstream,
                    brain_pipeline_git_commit=(
                        brain_pipeline_git_commit
                    ),
                ),
            ),
        )

        if not result.chain_accounting_valid:
            raise BrainRunnerError(
                "Stage-5 single-record accounting failed"
            )

        return result

    result = BrainRecordResult(
        brain_records=(
            _brain_defined_record(
                row,
                brain=brain,
                upstream=upstream,
                brain_pipeline_git_commit=(
                    brain_pipeline_git_commit
                ),
            ),
        ),
    )

    if not result.chain_accounting_valid:
        raise BrainRunnerError(
            "Stage-5 single-record accounting failed"
        )

    return result


# ---------------------------------------------------------------------------
# Stage-5 deterministic physical task partitioning
# ---------------------------------------------------------------------------

DEFAULT_BRAIN_ROW_GROUPS_PER_TASK = 64


@dataclass(frozen=True)
class BrainTaskPartition:
    """One contiguous physical partition of canonical Stage-3 BRI."""

    task_id: int
    task_count: int

    start_row_group: int
    stop_row_group: int  # exclusive

    row_group_count: int
    input_bri_chain_count: int


def _positive_partition_integer(
    value: object,
    *,
    field: str,
) -> int:
    if (
        not isinstance(value, int)
        or isinstance(value, bool)
        or value < 1
    ):
        raise BrainRunnerError(
            f"{field} must be a positive integer"
        )

    return value


def brain_task_count(
    row_group_count: int,
    *,
    row_groups_per_task: int = (
        DEFAULT_BRAIN_ROW_GROUPS_PER_TASK
    ),
) -> int:
    """Return deterministic Stage-5 task count."""

    row_group_count = _positive_partition_integer(
        row_group_count,
        field="row_group_count",
    )
    row_groups_per_task = _positive_partition_integer(
        row_groups_per_task,
        field="row_groups_per_task",
    )

    return (
        row_group_count
        + row_groups_per_task
        - 1
    ) // row_groups_per_task


def brain_task_partition(
    row_group_row_counts: tuple[int, ...],
    *,
    task_id: int,
    row_groups_per_task: int = (
        DEFAULT_BRAIN_ROW_GROUPS_PER_TASK
    ),
) -> BrainTaskPartition:
    """Describe one contiguous canonical-BRI row-group partition."""

    counts = tuple(
        row_group_row_counts
    )

    if not counts:
        raise BrainRunnerError(
            "Canonical Stage-3 BRI has no Parquet row groups"
        )

    for index, count in enumerate(counts):
        if (
            not isinstance(count, int)
            or isinstance(count, bool)
            or count < 1
        ):
            raise BrainRunnerError(
                "Invalid canonical Stage-3 BRI row-group "
                f"row count at index {index}: {count!r}"
            )

    row_groups_per_task = _positive_partition_integer(
        row_groups_per_task,
        field="row_groups_per_task",
    )

    if (
        not isinstance(task_id, int)
        or isinstance(task_id, bool)
        or task_id < 0
    ):
        raise BrainRunnerError(
            "task_id must be a non-negative integer"
        )

    task_count = brain_task_count(
        len(counts),
        row_groups_per_task=row_groups_per_task,
    )

    if task_id >= task_count:
        raise BrainRunnerError(
            "Stage-5 Brain task_id is outside the "
            "physical partition range"
        )

    start = task_id * row_groups_per_task
    stop = min(
        start + row_groups_per_task,
        len(counts),
    )

    selected = counts[
        start:stop
    ]

    return BrainTaskPartition(
        task_id=task_id,
        task_count=task_count,
        start_row_group=start,
        stop_row_group=stop,
        row_group_count=len(selected),
        input_bri_chain_count=sum(selected),
    )


# ---------------------------------------------------------------------------
# Stage-5 task-level Brain processing
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class BrainTables:
    """Schema-enforced terminal tables for one Stage-5 task."""

    chains: pa.Table
    undefined: pa.Table
    processing_errors: pa.Table


@dataclass(frozen=True)
class BrainBatchResult:
    """Terminal Stage-5 outcomes for one physical BRI partition."""

    partition: BrainTaskPartition
    tables: BrainTables

    @property
    def input_bri_chain_count(self) -> int:
        return self.partition.input_bri_chain_count

    @property
    def brain_chain_count(self) -> int:
        return self.tables.chains.num_rows

    @property
    def undefined_chain_count(self) -> int:
        return self.tables.undefined.num_rows

    @property
    def processing_error_count(self) -> int:
        return self.tables.processing_errors.num_rows

    @property
    def chain_accounting_valid(self) -> bool:
        return self.input_bri_chain_count == (
            self.brain_chain_count
            + self.undefined_chain_count
            + self.processing_error_count
        )


def _materialize_brain_tables(
    brain_records: list[Mapping[str, Any]],
    undefined_records: list[Mapping[str, Any]],
    processing_errors: list[Mapping[str, Any]],
) -> BrainTables:
    """Materialize Stage-5 terminal records against frozen schemas."""

    try:
        chains = pa.Table.from_pylist(
            list(brain_records),
            schema=STAGE5_BRAIN_CHAIN_SCHEMA,
        )
        undefined = pa.Table.from_pylist(
            list(undefined_records),
            schema=STAGE5_BRAIN_UNDEFINED_CHAIN_SCHEMA,
        )
        errors = pa.Table.from_pylist(
            list(processing_errors),
            schema=STAGE5_BRAIN_PROCESSING_ERROR_SCHEMA,
        )
    except (
        pa.ArrowException,
        TypeError,
        ValueError,
    ) as exc:
        raise BrainRunnerError(
            f"Cannot materialize Stage-5 Brain tables: {exc}"
        ) from exc

    expected = (
        (
            "Brain-defined",
            chains,
            STAGE5_BRAIN_CHAIN_SCHEMA,
        ),
        (
            "Brain-undefined",
            undefined,
            STAGE5_BRAIN_UNDEFINED_CHAIN_SCHEMA,
        ),
        (
            "Brain processing-error",
            errors,
            STAGE5_BRAIN_PROCESSING_ERROR_SCHEMA,
        ),
    )

    for name, table, schema in expected:
        if not table.schema.equals(
            schema,
            check_metadata=True,
        ):
            raise BrainRunnerError(
                f"Materialized {name} table schema mismatch"
            )

    return BrainTables(
        chains=chains,
        undefined=undefined,
        processing_errors=errors,
    )


def process_brain_partition(
    *,
    upstream: UpstreamBRI,
    partition: BrainTaskPartition,
    brain_pipeline_git_commit: str,
    record_processor: Callable[
        ...,
        BrainRecordResult,
    ] = process_brain_record,
) -> BrainBatchResult:
    """Process exactly one physical partition of canonical Stage-3 BRI."""

    brain_pipeline_git_commit = _validated_git_commit(
        brain_pipeline_git_commit,
        field="brain_pipeline_git_commit",
    )

    try:
        parquet = pq.ParquetFile(
            upstream.bri_path
        )
    except Exception as exc:
        raise BrainRunnerError(
            f"Cannot open canonical Stage-3 BRI population: {exc}"
        ) from exc

    total_row_groups = parquet.metadata.num_row_groups

    if (
        partition.start_row_group < 0
        or partition.stop_row_group
        > total_row_groups
        or partition.start_row_group
        >= partition.stop_row_group
    ):
        raise BrainRunnerError(
            "Stage-5 Brain partition is incompatible with "
            "canonical Stage-3 BRI row groups"
        )

    observed_partition_rows = sum(
        parquet.metadata.row_group(index).num_rows
        for index in range(
            partition.start_row_group,
            partition.stop_row_group,
        )
    )

    if (
        observed_partition_rows
        != partition.input_bri_chain_count
    ):
        raise BrainRunnerError(
            "Stage-5 Brain partition row accounting does not match "
            "canonical Stage-3 BRI metadata"
        )

    brain_records: list[Mapping[str, Any]] = []
    undefined_records: list[Mapping[str, Any]] = []
    processing_errors: list[Mapping[str, Any]] = []

    seen_identities: set[
        tuple[str, str, int, str]
    ] = set()

    input_row_count = 0

    try:
        batches = parquet.iter_batches(
            batch_size=64,
            row_groups=list(
                range(
                    partition.start_row_group,
                    partition.stop_row_group,
                )
            ),
        )

        for batch in batches:
            for row in batch.to_pylist():
                input_row_count += 1

                identity = _brain_row_identity(
                    row
                )

                if identity in seen_identities:
                    raise BrainRunnerError(
                        "Stage-5 Brain task contains duplicate "
                        f"canonical BRI identity: {identity!r}"
                    )

                seen_identities.add(
                    identity
                )

                result = record_processor(
                    row,
                    upstream=upstream,
                    brain_pipeline_git_commit=(
                        brain_pipeline_git_commit
                    ),
                )

                if (
                    result.input_bri_chain_count != 1
                    or not result.chain_accounting_valid
                ):
                    raise BrainRunnerError(
                        "Stage-5 per-record Brain accounting failed"
                    )

                brain_records.extend(
                    result.brain_records
                )
                undefined_records.extend(
                    result.undefined_records
                )
                processing_errors.extend(
                    result.processing_errors
                )

    except BrainRunnerError:
        raise
    except Exception as exc:
        raise BrainRunnerError(
            f"Cannot read Stage-5 canonical BRI partition: {exc}"
        ) from exc

    if input_row_count != partition.input_bri_chain_count:
        raise BrainRunnerError(
            "Stage-5 Brain task consumed an unexpected number "
            "of canonical BRI rows"
        )

    tables = _materialize_brain_tables(
        brain_records,
        undefined_records,
        processing_errors,
    )

    batch_result = BrainBatchResult(
        partition=partition,
        tables=tables,
    )

    if not batch_result.chain_accounting_valid:
        raise BrainRunnerError(
            "Stage-5 Brain task chain accounting failed"
        )

    return batch_result
