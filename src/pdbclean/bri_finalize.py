"""Validate task-level artifacts from distributed Stage-3 BRI production."""

from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
import re
from typing import Any, Iterable, Mapping

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

from pdbclean.bri_runner import (
    BRI_TASK_SUMMARY_SCHEMA_NAME,
    BRI_TASK_SUMMARY_SCHEMA_VERSION,
)
from pdbclean.schemas import (
    STAGE3_BRI_CHAIN_SCHEMA,
    STAGE3_ELIGIBLE_CHAIN_SCHEMA,
    STAGE3_BRI_PROCESSING_ERROR_SCHEMA,
)


class BRIFinalizeError(RuntimeError):
    """Raised when Stage-3 BRI finalization cannot proceed safely."""


_TASK_SUMMARY_RE = re.compile(r"^task_(\d+)\.json$")
_TASK_SHARD_RE = re.compile(r"^task_(\d+)\.parquet$")


@dataclass(frozen=True)
class BRITaskArtifacts:
    """Validated terminal artifacts for one Stage-3 logical task."""

    task_id: int
    chains_path: Path
    processing_errors_path: Path
    summary_path: Path
    summary: dict[str, Any]


def _parse_task_id(
    path: Path,
    pattern: re.Pattern[str],
    *,
    artifact_type: str,
) -> int:
    match = pattern.fullmatch(path.name)

    if match is None:
        raise BRIFinalizeError(
            f"Unexpected {artifact_type} filename: {path.name!r}"
        )

    return int(match.group(1))


def _discover_task_ids(
    directory: Path,
    *,
    suffix: str,
    pattern: re.Pattern[str],
    artifact_type: str,
) -> set[int]:
    if not directory.is_dir():
        return set()

    task_ids: set[int] = set()

    for path in sorted(directory.glob(f"*{suffix}")):
        task_id = _parse_task_id(
            path,
            pattern,
            artifact_type=artifact_type,
        )

        if task_id in task_ids:
            raise BRIFinalizeError(
                f"Duplicate {artifact_type} task ID: {task_id}"
            )

        task_ids.add(task_id)

    return task_ids


def _validate_exact_task_ids(
    observed: set[int],
    expected: set[int],
    *,
    artifact_type: str,
) -> None:
    missing = sorted(expected - observed)
    unexpected = sorted(observed - expected)

    if not missing and not unexpected:
        return

    pieces = []

    if missing:
        pieces.append(f"missing={missing}")

    if unexpected:
        pieces.append(f"unexpected={unexpected}")

    raise BRIFinalizeError(
        f"Invalid {artifact_type} task set: "
        + ", ".join(pieces)
    )


def _parquet_child_field_matches(
    observed: pa.Field,
    expected: pa.Field,
) -> bool:
    """Compare nested Parquet fields while ignoring child-field names."""

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
    """Compare Arrow types while tolerating Parquet list-child renaming."""

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


def _parquet_field_matches(
    observed: pa.Field,
    expected: pa.Field,
) -> bool:
    if (
        observed.name != expected.name
        or observed.nullable != expected.nullable
        or observed.metadata != expected.metadata
    ):
        return False

    return _parquet_type_matches(
        observed.type,
        expected.type,
    )


def _validate_shard_schema(
    path: Path,
    expected_schema: pa.Schema,
) -> None:
    try:
        observed_schema = pq.read_schema(path)
    except Exception as exc:
        raise BRIFinalizeError(
            f"Cannot read Parquet schema from {path}: {exc}"
        ) from exc

    if observed_schema.metadata != expected_schema.metadata:
        raise BRIFinalizeError(
            f"Unexpected Parquet schema metadata for {path}"
        )

    if len(observed_schema) != len(expected_schema):
        raise BRIFinalizeError(
            f"Unexpected Parquet field count for {path}"
        )

    if any(
        not _parquet_field_matches(observed, expected)
        for observed, expected in zip(
            observed_schema,
            expected_schema,
            strict=True,
        )
    ):
        raise BRIFinalizeError(
            f"Unexpected Parquet schema for {path}"
        )


def _nonnegative_int(
    summary: dict[str, Any],
    field: str,
    *,
    task_id: int,
) -> int:
    value = summary.get(field)

    if (
        not isinstance(value, int)
        or isinstance(value, bool)
        or value < 0
    ):
        raise BRIFinalizeError(
            f"Invalid {field} in Stage-3 BRI task "
            f"{task_id} summary"
        )

    return value


def _load_task_summary(
    path: Path,
    *,
    expected_task_id: int,
    expected_snapshot: str,
    expected_cleaning_protocol: str,
    expected_quality_pipeline_git_commit: str,
    expected_geometric_validation_pipeline_git_commit: str,
    expected_geometric_validation_finalizer_git_commit: str,
    expected_bri_pipeline_git_commit: str,
) -> dict[str, Any]:
    try:
        summary = json.loads(
            path.read_text(encoding="utf-8")
        )
    except (
        OSError,
        UnicodeError,
        json.JSONDecodeError,
    ) as exc:
        raise BRIFinalizeError(
            f"Cannot read Stage-3 BRI task summary "
            f"{path}: {exc}"
        ) from exc

    if not isinstance(summary, dict):
        raise BRIFinalizeError(
            "Stage-3 BRI task summary must contain a "
            f"JSON object: {path}"
        )

    if (
        summary.get("summary_schema_name")
        != BRI_TASK_SUMMARY_SCHEMA_NAME
    ):
        raise BRIFinalizeError(
            f"Unexpected summary schema name in {path}"
        )

    if (
        summary.get("summary_schema_version")
        != BRI_TASK_SUMMARY_SCHEMA_VERSION
    ):
        raise BRIFinalizeError(
            f"Unexpected summary schema version in {path}"
        )

    if str(summary.get("task_id", "")) != str(
        expected_task_id
    ):
        raise BRIFinalizeError(
            "Summary task ID does not match filename for "
            f"Stage-3 BRI task {expected_task_id}"
        )

    if summary.get("snapshot") != expected_snapshot:
        raise BRIFinalizeError(
            "Snapshot mismatch in Stage-3 BRI task "
            f"{expected_task_id}"
        )

    if (
        summary.get("cleaning_protocol")
        != expected_cleaning_protocol
    ):
        raise BRIFinalizeError(
            "Cleaning protocol mismatch in Stage-3 BRI task "
            f"{expected_task_id}"
        )

    provenance = (
        (
            "quality_pipeline_git_commit",
            expected_quality_pipeline_git_commit,
            "Quality producer",
        ),
        (
            "geometric_validation_pipeline_git_commit",
            expected_geometric_validation_pipeline_git_commit,
            "Geometry producer",
        ),
        (
            "geometric_validation_finalizer_git_commit",
            expected_geometric_validation_finalizer_git_commit,
            "Geometry finalizer",
        ),
        (
            "bri_pipeline_git_commit",
            expected_bri_pipeline_git_commit,
            "BRI producer",
        ),
    )

    for field, expected, label in provenance:
        if summary.get(field) != expected:
            raise BRIFinalizeError(
                f"{label} Git commit mismatch in "
                f"Stage-3 BRI task {expected_task_id}"
            )

    if summary.get("chain_accounting_valid") is not True:
        raise BRIFinalizeError(
            "Chain accounting is invalid for Stage-3 BRI "
            f"task {expected_task_id}"
        )

    return summary


def discover_bri_task_artifacts(
    bri_root: str | Path,
    *,
    expected_task_ids: Iterable[int],
    expected_snapshot: str,
    expected_cleaning_protocol: str,
    expected_quality_pipeline_git_commit: str,
    expected_geometric_validation_pipeline_git_commit: str,
    expected_geometric_validation_finalizer_git_commit: str,
    expected_bri_pipeline_git_commit: str,
) -> tuple[BRITaskArtifacts, ...]:
    """Discover and validate the complete Stage-3 task artifact set."""

    root = Path(bri_root)

    expected_ids = tuple(
        sorted(set(expected_task_ids))
    )
    expected_set = set(expected_ids)

    if expected_ids != tuple(range(len(expected_ids))):
        raise BRIFinalizeError(
            "Expected task IDs must be contiguous and zero-based"
        )

    temporary_files = sorted(
        root.rglob("*.tmp")
    )

    if temporary_files:
        raise BRIFinalizeError(
            "Temporary Stage-3 BRI files remain: "
            + ", ".join(
                str(path)
                for path in temporary_files
            )
        )

    artifact_specs = (
        (
            root / "chains",
            ".parquet",
            _TASK_SHARD_RE,
            "BRI chain shard",
        ),
        (
            root / "processing_errors",
            ".parquet",
            _TASK_SHARD_RE,
            "BRI processing-error shard",
        ),
        (
            root / "summaries",
            ".json",
            _TASK_SUMMARY_RE,
            "BRI task summary",
        ),
    )

    for directory, suffix, pattern, artifact_type in artifact_specs:
        observed = _discover_task_ids(
            directory,
            suffix=suffix,
            pattern=pattern,
            artifact_type=artifact_type,
        )

        _validate_exact_task_ids(
            observed,
            expected_set,
            artifact_type=artifact_type,
        )

    artifacts: list[BRITaskArtifacts] = []

    for task_id in expected_ids:
        chains_path = (
            root
            / "chains"
            / f"task_{task_id}.parquet"
        )

        errors_path = (
            root
            / "processing_errors"
            / f"task_{task_id}.parquet"
        )

        summary_path = (
            root
            / "summaries"
            / f"task_{task_id}.json"
        )

        _validate_shard_schema(
            chains_path,
            STAGE3_BRI_CHAIN_SCHEMA,
        )

        _validate_shard_schema(
            errors_path,
            STAGE3_BRI_PROCESSING_ERROR_SCHEMA,
        )

        summary = _load_task_summary(
            summary_path,
            expected_task_id=task_id,
            expected_snapshot=expected_snapshot,
            expected_cleaning_protocol=(
                expected_cleaning_protocol
            ),
            expected_quality_pipeline_git_commit=(
                expected_quality_pipeline_git_commit
            ),
            expected_geometric_validation_pipeline_git_commit=(
                expected_geometric_validation_pipeline_git_commit
            ),
            expected_geometric_validation_finalizer_git_commit=(
                expected_geometric_validation_finalizer_git_commit
            ),
            expected_bri_pipeline_git_commit=(
                expected_bri_pipeline_git_commit
            ),
        )

        artifacts.append(
            BRITaskArtifacts(
                task_id=task_id,
                chains_path=chains_path,
                processing_errors_path=errors_path,
                summary_path=summary_path,
                summary=summary,
            )
        )

    return tuple(artifacts)


def validate_bri_task_accounting(
    artifacts: Iterable[BRITaskArtifacts],
) -> None:
    """Recompute task-level Stage-3 terminal row accounting."""

    for artifact in artifacts:
        task_id = artifact.task_id
        summary = artifact.summary

        input_count = _nonnegative_int(
            summary,
            "input_eligible_chain_count",
            task_id=task_id,
        )

        bri_count = _nonnegative_int(
            summary,
            "bri_chain_count",
            task_id=task_id,
        )

        error_count = _nonnegative_int(
            summary,
            "processing_error_count",
            task_id=task_id,
        )

        chains_rows = pq.read_metadata(
            artifact.chains_path
        ).num_rows

        error_rows = pq.read_metadata(
            artifact.processing_errors_path
        ).num_rows

        if chains_rows != bri_count:
            raise BRIFinalizeError(
                "BRI chain row-count mismatch for Stage-3 "
                f"task {task_id}: parquet={chains_rows}, "
                f"summary={bri_count}"
            )

        if error_rows != error_count:
            raise BRIFinalizeError(
                "Processing-error row-count mismatch for "
                f"Stage-3 task {task_id}: parquet={error_rows}, "
                f"summary={error_count}"
            )

        if input_count != chains_rows + error_rows:
            raise BRIFinalizeError(
                "Terminal chain accounting mismatch for "
                f"Stage-3 BRI task {task_id}"
            )


@dataclass(frozen=True)
class BRIGlobalValidation:
    """Globally recomputed Stage-3 BRI identity and scientific state."""

    eligible_chain_count: int
    bri_chain_count: int
    processing_error_count: int

    unique_eligible_identity_count: int
    unique_bri_identity_count: int

    minimum_retained_residue_count: int
    maximum_retained_residue_count: int


def _chain_identity(
    row: Mapping[str, Any],
) -> tuple[str, str, int, str]:
    """Return the canonical Stage-3 chain identity."""

    try:
        return (
            row["snapshot"],
            row["pdb_id"],
            row["model_id"],
            row["label_chain_id"],
        )
    except KeyError as exc:
        raise BRIFinalizeError(
            "Stage-3 row is missing canonical identity lineage"
        ) from exc


def validate_bri_global_state(
    artifacts: Iterable[BRITaskArtifacts],
    *,
    eligible_path: str | Path,
    expected_quality_pipeline_git_commit: str,
    expected_geometric_validation_pipeline_git_commit: str,
    expected_geometric_validation_finalizer_git_commit: str,
    expected_bri_pipeline_git_commit: str,
) -> BRIGlobalValidation:
    """Validate global Stage-2-to-BRI identity, lineage, and BRI state."""

    artifacts = tuple(
        sorted(
            artifacts,
            key=lambda artifact: artifact.task_id,
        )
    )

    validate_bri_task_accounting(artifacts)

    eligible_path = Path(eligible_path)

    _validate_shard_schema(
        eligible_path,
        STAGE3_ELIGIBLE_CHAIN_SCHEMA,
    )

    eligible = pq.read_table(
        eligible_path,
        columns=[
            "snapshot",
            "pdb_id",
            "model_id",
            "label_chain_id",
            "retained_residue_count",
            "retained_label_seq_ids",
            "source_mmcif_key",
            "source_etag",
            "cleaning_protocol",
            "pipeline_git_commit",
        ],
    )

    eligible_identity_to_index: dict[
        tuple[str, str, int, str],
        int,
    ] = {}

    absolute_index = 0

    for batch in eligible.select(
        [
            "snapshot",
            "pdb_id",
            "model_id",
            "label_chain_id",
        ]
    ).to_batches(max_chunksize=65_536):
        for row in batch.to_pylist():
            identity = _chain_identity(row)

            if identity in eligible_identity_to_index:
                raise BRIFinalizeError(
                    "Duplicate canonical Stage-2 eligible identity: "
                    f"{identity!r}"
                )

            eligible_identity_to_index[
                identity
            ] = absolute_index

            absolute_index += 1

    if absolute_index != eligible.num_rows:
        raise BRIFinalizeError(
            "Canonical eligible identity indexing failed"
        )

    eligible_m = eligible["retained_residue_count"]
    eligible_ids = eligible["retained_label_seq_ids"]
    eligible_source = eligible["source_mmcif_key"]
    eligible_etag = eligible["source_etag"]
    eligible_protocol = eligible["cleaning_protocol"]
    eligible_quality_commit = eligible["pipeline_git_commit"]

    seen_bri: set[
        tuple[str, str, int, str]
    ] = set()

    bri_count = 0
    processing_error_count = 0

    minimum_m: int | None = None
    maximum_m: int | None = None

    for artifact in artifacts:
        task_id = artifact.task_id

        task_error_count = pq.read_metadata(
            artifact.processing_errors_path
        ).num_rows

        processing_error_count += task_error_count

        if task_error_count != 0:
            raise BRIFinalizeError(
                "Cannot globally validate Stage-3 BRI while "
                f"processing errors remain in task {task_id}: "
                f"{task_error_count}"
            )

        parquet_file = pq.ParquetFile(
            artifact.chains_path
        )

        for batch in parquet_file.iter_batches(
            batch_size=64,
        ):
            for row in batch.to_pylist():
                identity = _chain_identity(row)

                if identity in seen_bri:
                    raise BRIFinalizeError(
                        "Duplicate BRI-chain identity across "
                        f"Stage-3 tasks: {identity!r}"
                    )

                seen_bri.add(identity)

                eligible_index = (
                    eligible_identity_to_index.get(
                        identity
                    )
                )

                if eligible_index is None:
                    raise BRIFinalizeError(
                        "Unexpected BRI-chain identity not present "
                        "in canonical Stage-2 eligible population: "
                        f"{identity!r}"
                    )

                expected_m = eligible_m[
                    eligible_index
                ].as_py()

                if row["retained_residue_count"] != expected_m:
                    raise BRIFinalizeError(
                        "Retained residue-count lineage mismatch "
                        f"for {identity!r}"
                    )

                if expected_m < 1:
                    raise BRIFinalizeError(
                        "Canonical retained residue count must be "
                        f"positive for {identity!r}"
                    )

                expected_ids = eligible_ids[
                    eligible_index
                ].as_py()

                if row["retained_label_seq_ids"] != expected_ids:
                    raise BRIFinalizeError(
                        "Retained residue-ID lineage mismatch "
                        f"for {identity!r}"
                    )

                lineage_comparisons = (
                    (
                        "source_mmcif_key",
                        eligible_source[
                            eligible_index
                        ].as_py(),
                    ),
                    (
                        "source_etag",
                        eligible_etag[
                            eligible_index
                        ].as_py(),
                    ),
                    (
                        "cleaning_protocol",
                        eligible_protocol[
                            eligible_index
                        ].as_py(),
                    ),
                )

                for field, expected in lineage_comparisons:
                    if row[field] != expected:
                        raise BRIFinalizeError(
                            "Stage-2/Stage-3 lineage mismatch for "
                            f"{identity!r}: {field}"
                        )

                upstream_quality_commit = (
                    eligible_quality_commit[
                        eligible_index
                    ].as_py()
                )

                if (
                    upstream_quality_commit
                    != expected_quality_pipeline_git_commit
                ):
                    raise BRIFinalizeError(
                        "Canonical Stage-2 eligible quality "
                        f"provenance mismatch for {identity!r}"
                    )

                provenance_comparisons = (
                    (
                        "quality_pipeline_git_commit",
                        expected_quality_pipeline_git_commit,
                    ),
                    (
                        "geometric_validation_pipeline_git_commit",
                        expected_geometric_validation_pipeline_git_commit,
                    ),
                    (
                        "geometric_validation_finalizer_git_commit",
                        expected_geometric_validation_finalizer_git_commit,
                    ),
                    (
                        "bri_pipeline_git_commit",
                        expected_bri_pipeline_git_commit,
                    ),
                )

                for field, expected in provenance_comparisons:
                    if row[field] != expected:
                        raise BRIFinalizeError(
                            "Stage-3 BRI row provenance mismatch "
                            f"for {identity!r}: {field}"
                        )

                matrix = np.asarray(
                    row["bri"],
                    dtype=np.float64,
                )

                if matrix.shape != (expected_m, 9):
                    raise BRIFinalizeError(
                        "Definition 3.4 BRI shape mismatch for "
                        f"{identity!r}: observed={matrix.shape!r}, "
                        f"expected={(expected_m, 9)!r}"
                    )

                if not np.isfinite(matrix).all():
                    raise BRIFinalizeError(
                        "Non-finite Definition 3.4 BRI coordinate "
                        f"for {identity!r}"
                    )

                if not np.array_equal(
                    matrix,
                    np.around(matrix, 3),
                ):
                    raise BRIFinalizeError(
                        "Non-canonical Definition 3.4 BRI "
                        f"precision for {identity!r}"
                    )

                # Definition 3.4 row 1 is represented in its own
                # residue basis. Only x(N), x(C), and y(C) can be
                # non-zero. Therefore these six coordinates are
                # exactly zero:
                #
                # y(N), z(N), x(A), y(A), z(A), z(C)
                first_row_zero_indices = (
                    1,
                    2,
                    3,
                    4,
                    5,
                    8,
                )

                if not np.array_equal(
                    matrix[
                        0,
                        first_row_zero_indices,
                    ],
                    np.zeros(
                        len(first_row_zero_indices),
                        dtype=np.float64,
                    ),
                ):
                    raise BRIFinalizeError(
                        "Definition 3.4 first-row zero structure "
                        f"mismatch for {identity!r}"
                    )

                bri_count += 1

                minimum_m = (
                    expected_m
                    if minimum_m is None
                    else min(minimum_m, expected_m)
                )

                maximum_m = (
                    expected_m
                    if maximum_m is None
                    else max(maximum_m, expected_m)
                )

    eligible_count = eligible.num_rows

    if bri_count != len(seen_bri):
        raise BRIFinalizeError(
            "Global BRI identity uniqueness accounting failed"
        )

    if len(seen_bri) != eligible_count:
        missing = []

        for identity in eligible_identity_to_index:
            if identity not in seen_bri:
                missing.append(identity)

                if len(missing) == 10:
                    break

        raise BRIFinalizeError(
            "Canonical Stage-2 eligible / Stage-3 BRI identity "
            "population mismatch: "
            f"eligible={eligible_count}, "
            f"BRI={len(seen_bri)}, "
            f"first_missing={missing!r}"
        )

    if processing_error_count != 0:
        raise BRIFinalizeError(
            "Cannot globally validate Stage-3 BRI with "
            "processing errors"
        )

    if minimum_m is None or maximum_m is None:
        raise BRIFinalizeError(
            "Stage-3 BRI population is unexpectedly empty"
        )

    return BRIGlobalValidation(
        eligible_chain_count=eligible_count,
        bri_chain_count=bri_count,
        processing_error_count=processing_error_count,
        unique_eligible_identity_count=len(
            eligible_identity_to_index
        ),
        unique_bri_identity_count=len(seen_bri),
        minimum_retained_residue_count=minimum_m,
        maximum_retained_residue_count=maximum_m,
    )
