"""Finalize a complete distributed post-cleaning geometric-validation stage."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
import json
from pathlib import Path
import re
from typing import Any, Iterable

import pyarrow as pa
import pyarrow.parquet as pq

from pdbclean.geometric_validation_runner import (
    GEOMETRIC_VALIDATION_TASK_SUMMARY_SCHEMA_NAME,
    GEOMETRIC_VALIDATION_TASK_SUMMARY_SCHEMA_VERSION,
)
from pdbclean.schemas import (
    GEOMETRIC_VALIDATION_AUDIT_SCHEMA,
    GEOMETRIC_VALIDATION_PROCESSING_ERROR_SCHEMA,
    GOLD_ACCEPTED_CHAIN_SCHEMA,
    STAGE3_ELIGIBLE_CHAIN_SCHEMA,
    STAGE3_QUARANTINED_CHAIN_SCHEMA,
)


class GeometricValidationFinalizeError(RuntimeError):
    """Raised when Step-2 finalization cannot proceed safely."""


_TASK_SUMMARY_RE = re.compile(r"^task_(\d+)\.json$")
_TASK_SHARD_RE = re.compile(r"^task_(\d+)\.parquet$")


@dataclass(frozen=True)
class GeometricValidationTaskArtifacts:
    """Validated filesystem artifacts for one Step-2 logical task."""

    task_id: int
    accepted_path: Path
    audit_path: Path
    error_path: Path
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
        raise GeometricValidationFinalizeError(
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
            raise GeometricValidationFinalizeError(
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

    raise GeometricValidationFinalizeError(
        f"Invalid {artifact_type} task set: "
        + ", ".join(pieces)
    )


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

    if pa.types.is_list(observed_type) and pa.types.is_list(expected_type):
        observed_value = observed_type.value_field
        expected_value = expected_type.value_field

        return (
            observed_value.type == expected_value.type
            and observed_value.nullable == expected_value.nullable
            and observed_value.metadata == expected_value.metadata
        )

    return observed_type == expected_type


def _validate_shard_schema(
    path: Path,
    expected_schema: pa.Schema,
) -> None:
    try:
        observed_schema = pq.read_schema(path)
    except Exception as exc:
        raise GeometricValidationFinalizeError(
            f"Cannot read Parquet schema from {path}: {exc}"
        ) from exc

    if observed_schema.metadata != expected_schema.metadata:
        raise GeometricValidationFinalizeError(
            f"Unexpected Parquet schema metadata for {path}"
        )

    if len(observed_schema) != len(expected_schema):
        raise GeometricValidationFinalizeError(
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
        raise GeometricValidationFinalizeError(
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
        raise GeometricValidationFinalizeError(
            f"Invalid {field} in Step-2 task {task_id} summary"
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
    expected_minimum_backbone_distance_angstrom: float,
    expected_minimum_triangle_angle_degrees: float,
) -> dict[str, Any]:
    try:
        summary = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise GeometricValidationFinalizeError(
            f"Cannot read Step-2 task summary {path}: {exc}"
        ) from exc

    if not isinstance(summary, dict):
        raise GeometricValidationFinalizeError(
            f"Step-2 task summary must contain a JSON object: {path}"
        )

    if (
        summary.get("summary_schema_name")
        != GEOMETRIC_VALIDATION_TASK_SUMMARY_SCHEMA_NAME
    ):
        raise GeometricValidationFinalizeError(
            f"Unexpected summary schema name in {path}"
        )

    if (
        summary.get("summary_schema_version")
        != GEOMETRIC_VALIDATION_TASK_SUMMARY_SCHEMA_VERSION
    ):
        raise GeometricValidationFinalizeError(
            f"Unexpected summary schema version in {path}"
        )

    if str(summary.get("task_id", "")) != str(expected_task_id):
        raise GeometricValidationFinalizeError(
            f"Summary task ID does not match filename for task "
            f"{expected_task_id}"
        )

    if summary.get("snapshot") != expected_snapshot:
        raise GeometricValidationFinalizeError(
            f"Snapshot mismatch in Step-2 task {expected_task_id}"
        )

    if summary.get("cleaning_protocol") != expected_cleaning_protocol:
        raise GeometricValidationFinalizeError(
            f"Cleaning protocol mismatch in Step-2 task "
            f"{expected_task_id}"
        )

    if (
        summary.get("quality_pipeline_git_commit")
        != expected_quality_pipeline_git_commit
    ):
        raise GeometricValidationFinalizeError(
            f"Quality producer Git commit mismatch in Step-2 task "
            f"{expected_task_id}"
        )

    if (
        summary.get("geometric_validation_pipeline_git_commit")
        != expected_geometric_validation_pipeline_git_commit
    ):
        raise GeometricValidationFinalizeError(
            f"Geometry producer Git commit mismatch in Step-2 task "
            f"{expected_task_id}"
        )

    if (
        summary.get(
            "configured_minimum_backbone_distance_angstrom"
        )
        != expected_minimum_backbone_distance_angstrom
    ):
        raise GeometricValidationFinalizeError(
            f"Backbone-distance threshold mismatch in Step-2 task "
            f"{expected_task_id}"
        )

    if (
        summary.get(
            "configured_minimum_triangle_angle_degrees"
        )
        != expected_minimum_triangle_angle_degrees
    ):
        raise GeometricValidationFinalizeError(
            f"Triangle-angle threshold mismatch in Step-2 task "
            f"{expected_task_id}"
        )

    if summary.get("chain_accounting_valid") is not True:
        raise GeometricValidationFinalizeError(
            f"Chain accounting is invalid for Step-2 task "
            f"{expected_task_id}"
        )

    return summary


def discover_geometric_validation_task_artifacts(
    quality_root: str | Path,
    geometric_validation_root: str | Path,
    *,
    expected_task_ids: Iterable[int],
    expected_snapshot: str,
    expected_cleaning_protocol: str,
    expected_quality_pipeline_git_commit: str,
    expected_geometric_validation_pipeline_git_commit: str,
    expected_minimum_backbone_distance_angstrom: float,
    expected_minimum_triangle_angle_degrees: float,
) -> tuple[GeometricValidationTaskArtifacts, ...]:
    """Discover and validate the complete task-level Stage-2 artifact set."""

    quality_root = Path(quality_root)
    geometry_root = Path(geometric_validation_root)

    expected_ids = tuple(sorted(set(expected_task_ids)))
    expected_set = set(expected_ids)

    if expected_ids != tuple(range(len(expected_ids))):
        raise GeometricValidationFinalizeError(
            "Expected task IDs must be contiguous and zero-based"
        )

    temporary_files = sorted(
        list(geometry_root.rglob("*.tmp"))
        + list((quality_root / "accepted").glob("*.tmp"))
    )

    if temporary_files:
        raise GeometricValidationFinalizeError(
            "Temporary upstream/finalization input files remain: "
            + ", ".join(str(path) for path in temporary_files)
        )

    artifact_specs = (
        (
            quality_root / "accepted",
            ".parquet",
            _TASK_SHARD_RE,
            "quality accepted shard",
        ),
        (
            geometry_root / "audit",
            ".parquet",
            _TASK_SHARD_RE,
            "geometric-validation audit shard",
        ),
        (
            geometry_root / "errors",
            ".parquet",
            _TASK_SHARD_RE,
            "geometric-validation error shard",
        ),
        (
            geometry_root / "summaries",
            ".json",
            _TASK_SUMMARY_RE,
            "geometric-validation summary",
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

    artifacts = []

    for task_id in expected_ids:
        accepted_path = (
            quality_root
            / "accepted"
            / f"task_{task_id}.parquet"
        )
        audit_path = (
            geometry_root
            / "audit"
            / f"task_{task_id}.parquet"
        )
        error_path = (
            geometry_root
            / "errors"
            / f"task_{task_id}.parquet"
        )
        summary_path = (
            geometry_root
            / "summaries"
            / f"task_{task_id}.json"
        )

        _validate_shard_schema(
            accepted_path,
            GOLD_ACCEPTED_CHAIN_SCHEMA,
        )
        _validate_shard_schema(
            audit_path,
            GEOMETRIC_VALIDATION_AUDIT_SCHEMA,
        )
        _validate_shard_schema(
            error_path,
            GEOMETRIC_VALIDATION_PROCESSING_ERROR_SCHEMA,
        )

        summary = _load_task_summary(
            summary_path,
            expected_task_id=task_id,
            expected_snapshot=expected_snapshot,
            expected_cleaning_protocol=expected_cleaning_protocol,
            expected_quality_pipeline_git_commit=(
                expected_quality_pipeline_git_commit
            ),
            expected_geometric_validation_pipeline_git_commit=(
                expected_geometric_validation_pipeline_git_commit
            ),
            expected_minimum_backbone_distance_angstrom=(
                expected_minimum_backbone_distance_angstrom
            ),
            expected_minimum_triangle_angle_degrees=(
                expected_minimum_triangle_angle_degrees
            ),
        )

        artifacts.append(
            GeometricValidationTaskArtifacts(
                task_id=task_id,
                accepted_path=accepted_path,
                audit_path=audit_path,
                error_path=error_path,
                summary_path=summary_path,
                summary=summary,
            )
        )

    return tuple(artifacts)


def validate_geometric_validation_task_accounting(
    artifacts: Iterable[GeometricValidationTaskArtifacts],
) -> None:
    """Recompute task-level Step-2 row/accounting invariants."""

    for artifact in artifacts:
        task_id = artifact.task_id
        summary = artifact.summary

        input_count = _nonnegative_int(
            summary,
            "input_accepted_chain_count",
            task_id=task_id,
        )
        audit_count = _nonnegative_int(
            summary,
            "audit_chain_count",
            task_id=task_id,
        )
        passed_count = _nonnegative_int(
            summary,
            "geometric_passed_chain_count",
            task_id=task_id,
        )
        violated_count = _nonnegative_int(
            summary,
            "geometric_violated_chain_count",
            task_id=task_id,
        )
        error_count = _nonnegative_int(
            summary,
            "processing_error_count",
            task_id=task_id,
        )

        accepted_rows = pq.read_metadata(
            artifact.accepted_path
        ).num_rows
        audit_rows = pq.read_metadata(
            artifact.audit_path
        ).num_rows
        error_rows = pq.read_metadata(
            artifact.error_path
        ).num_rows

        if accepted_rows != input_count:
            raise GeometricValidationFinalizeError(
                f"Accepted row-count mismatch for Step-2 task "
                f"{task_id}: parquet={accepted_rows}, "
                f"summary={input_count}"
            )

        if audit_rows != audit_count:
            raise GeometricValidationFinalizeError(
                f"Audit row-count mismatch for Step-2 task "
                f"{task_id}: parquet={audit_rows}, "
                f"summary={audit_count}"
            )

        if error_rows != error_count:
            raise GeometricValidationFinalizeError(
                f"Error row-count mismatch for Step-2 task "
                f"{task_id}: parquet={error_rows}, "
                f"summary={error_count}"
            )

        if input_count != audit_count + error_count:
            raise GeometricValidationFinalizeError(
                f"Terminal chain accounting mismatch for Step-2 task "
                f"{task_id}"
            )

        if audit_count != passed_count + violated_count:
            raise GeometricValidationFinalizeError(
                f"Audit pass/violation accounting mismatch for "
                f"Step-2 task {task_id}"
            )


GEOMETRIC_VALIDATION_GLOBAL_SUMMARY_SCHEMA_NAME = (
    "pdbclean_geometric_validation_global_summary"
)
GEOMETRIC_VALIDATION_GLOBAL_SUMMARY_SCHEMA_VERSION = "1.0"

GEOMETRIC_VALIDATION_SUCCESS_SCHEMA_NAME = (
    "pdbclean_geometric_validation_success"
)
GEOMETRIC_VALIDATION_SUCCESS_SCHEMA_VERSION = "1.0"


@dataclass(frozen=True)
class GeometricValidationGlobalValidation:
    """Globally recomputed Stage-2 scientific/accounting state."""

    input_accepted_chain_count: int
    audit_chain_count: int
    eligible_chain_count: int
    quarantined_chain_count: int
    processing_error_count: int
    unique_chain_identity_count: int
    violation_event_count: int
    violations_by_type: dict[str, int]


@dataclass(frozen=True)
class GeometricValidationFinalizePublication:
    """Canonical Stage-3 population and completion metadata."""

    eligible_path: Path
    quarantined_path: Path
    global_summary_path: Path
    success_path: Path
    global_summary: dict[str, Any]


def _chain_identity(
    row: dict[str, Any],
) -> tuple[str, str, int, str]:
    return (
        row["snapshot"],
        row["pdb_id"],
        row["model_id"],
        row["label_chain_id"],
    )


def _validated_task_partition_rows(
    artifact: GeometricValidationTaskArtifacts,
) -> tuple[
    list[dict[str, Any]],
    list[dict[str, Any]],
    Counter[str],
    int,
]:
    """Join one accepted shard exactly to its Step-2 terminal audit."""

    summary = artifact.summary
    task_id = artifact.task_id

    error_rows = pq.read_table(
        artifact.error_path
    ).to_pylist()

    if error_rows:
        raise GeometricValidationFinalizeError(
            "Cannot finalize Step-2 while processing errors remain "
            f"in task {task_id}: {len(error_rows)}"
        )

    accepted_rows = pq.read_table(
        artifact.accepted_path
    ).to_pylist()

    audit_rows = pq.read_table(
        artifact.audit_path
    ).to_pylist()

    accepted_map: dict[
        tuple[str, str, int, str],
        dict[str, Any],
    ] = {}

    for row in accepted_rows:
        identity = _chain_identity(row)

        if identity in accepted_map:
            raise GeometricValidationFinalizeError(
                f"Duplicate accepted-chain identity in task "
                f"{task_id}: {identity!r}"
            )

        accepted_map[identity] = row

    audit_map: dict[
        tuple[str, str, int, str],
        dict[str, Any],
    ] = {}

    for row in audit_rows:
        identity = _chain_identity(row)

        if identity in audit_map:
            raise GeometricValidationFinalizeError(
                f"Duplicate audit-chain identity in task "
                f"{task_id}: {identity!r}"
            )

        audit_map[identity] = row

    accepted_ids = set(accepted_map)
    audit_ids = set(audit_map)

    if accepted_ids != audit_ids:
        missing_audit = sorted(accepted_ids - audit_ids)
        unexpected_audit = sorted(audit_ids - accepted_ids)

        raise GeometricValidationFinalizeError(
            f"Accepted/audit identity-set mismatch in task {task_id}: "
            f"missing_audit={missing_audit[:10]!r}, "
            f"unexpected_audit={unexpected_audit[:10]!r}"
        )

    expected_snapshot = summary["snapshot"]
    expected_protocol = summary["cleaning_protocol"]
    expected_quality_commit = summary[
        "quality_pipeline_git_commit"
    ]
    expected_geometry_commit = summary[
        "geometric_validation_pipeline_git_commit"
    ]
    expected_distance = summary[
        "configured_minimum_backbone_distance_angstrom"
    ]
    expected_angle = summary[
        "configured_minimum_triangle_angle_degrees"
    ]

    eligible: list[dict[str, Any]] = []
    quarantined: list[dict[str, Any]] = []

    violations = Counter()
    violation_event_count = 0

    for identity in accepted_map:
        upstream = accepted_map[identity]
        audit = audit_map[identity]

        # Step-1 row itself must agree with the completed task summary.
        if upstream["snapshot"] != expected_snapshot:
            raise GeometricValidationFinalizeError(
                f"Accepted snapshot mismatch for {identity!r}"
            )

        if upstream["cleaning_protocol"] != expected_protocol:
            raise GeometricValidationFinalizeError(
                f"Accepted cleaning protocol mismatch for "
                f"{identity!r}"
            )

        if upstream["pipeline_git_commit"] != expected_quality_commit:
            raise GeometricValidationFinalizeError(
                f"Accepted quality producer mismatch for "
                f"{identity!r}"
            )

        # Exact retained-chain lineage must survive Step 2 unchanged.
        comparisons = (
            (
                "retained_residue_count",
                upstream["retained_residue_count"],
                audit["retained_residue_count"],
            ),
            (
                "retained_label_seq_ids",
                upstream["retained_label_seq_ids"],
                audit["retained_label_seq_ids"],
            ),
            (
                "source_mmcif_key",
                upstream["source_mmcif_key"],
                audit["source_mmcif_key"],
            ),
            (
                "source_etag",
                upstream["source_etag"],
                audit["source_etag"],
            ),
            (
                "cleaning_protocol",
                upstream["cleaning_protocol"],
                audit["cleaning_protocol"],
            ),
            (
                "quality_pipeline_git_commit",
                upstream["pipeline_git_commit"],
                audit["quality_pipeline_git_commit"],
            ),
        )

        for field, expected, observed in comparisons:
            if observed != expected:
                raise GeometricValidationFinalizeError(
                    f"Step-1/Step-2 lineage mismatch for "
                    f"{identity!r}: {field}"
                )

        if (
            audit["geometric_validation_pipeline_git_commit"]
            != expected_geometry_commit
        ):
            raise GeometricValidationFinalizeError(
                f"Geometry producer mismatch for {identity!r}"
            )

        if (
            audit[
                "configured_minimum_backbone_distance_angstrom"
            ]
            != expected_distance
        ):
            raise GeometricValidationFinalizeError(
                f"Backbone-distance threshold mismatch for "
                f"{identity!r}"
            )

        if (
            audit[
                "configured_minimum_triangle_angle_degrees"
            ]
            != expected_angle
        ):
            raise GeometricValidationFinalizeError(
                f"Triangle-angle threshold mismatch for "
                f"{identity!r}"
            )

        violation_count = audit["violation_count"]
        violation_types = audit["violation_types"]
        violation_residue_ids = audit[
            "violation_residue_ids"
        ]
        violation_details = audit["violation_details"]

        if not (
            len(violation_types)
            == len(violation_residue_ids)
            == len(violation_details)
            == violation_count
        ):
            raise GeometricValidationFinalizeError(
                f"Violation evidence length mismatch for "
                f"{identity!r}"
            )

        if audit["passed"]:
            if violation_count != 0:
                raise GeometricValidationFinalizeError(
                    f"Passing chain has violations: {identity!r}"
                )

            eligible.append(upstream)

        else:
            if violation_count < 1:
                raise GeometricValidationFinalizeError(
                    f"Failing chain has no violation evidence: "
                    f"{identity!r}"
                )

            quarantined.append(upstream)
            violations.update(violation_types)
            violation_event_count += violation_count

    return (
        eligible,
        quarantined,
        violations,
        violation_event_count,
    )


def validate_geometric_validation_global_state(
    artifacts: Iterable[GeometricValidationTaskArtifacts],
) -> GeometricValidationGlobalValidation:
    """Recompute global identity, lineage, and eligibility state."""

    artifacts = tuple(
        sorted(
            artifacts,
            key=lambda artifact: artifact.task_id,
        )
    )

    seen: set[tuple[str, str, int, str]] = set()

    input_count = 0
    audit_count = 0
    eligible_count = 0
    quarantined_count = 0
    processing_error_count = 0
    violation_event_count = 0
    violations = Counter()

    for artifact in artifacts:
        accepted_rows = pq.read_table(
            artifact.accepted_path,
            columns=[
                "snapshot",
                "pdb_id",
                "model_id",
                "label_chain_id",
            ],
        ).to_pylist()

        for row in accepted_rows:
            identity = _chain_identity(row)

            if identity in seen:
                raise GeometricValidationFinalizeError(
                    f"Duplicate accepted-chain identity across "
                    f"tasks: {identity!r}"
                )

            seen.add(identity)

        (
            eligible,
            quarantined,
            task_violations,
            task_violation_events,
        ) = _validated_task_partition_rows(artifact)

        input_count += len(eligible) + len(quarantined)
        audit_count += (
            pq.read_metadata(artifact.audit_path).num_rows
        )
        processing_error_count += (
            pq.read_metadata(artifact.error_path).num_rows
        )

        eligible_count += len(eligible)
        quarantined_count += len(quarantined)
        violations.update(task_violations)
        violation_event_count += task_violation_events

    if processing_error_count != 0:
        raise GeometricValidationFinalizeError(
            "Cannot finalize Step-2 with processing errors"
        )

    if input_count != audit_count:
        raise GeometricValidationFinalizeError(
            "Global accepted/audit accounting mismatch"
        )

    if input_count != eligible_count + quarantined_count:
        raise GeometricValidationFinalizeError(
            "Global eligibility accounting mismatch"
        )

    if len(seen) != input_count:
        raise GeometricValidationFinalizeError(
            "Global accepted identity uniqueness failed"
        )

    return GeometricValidationGlobalValidation(
        input_accepted_chain_count=input_count,
        audit_chain_count=audit_count,
        eligible_chain_count=eligible_count,
        quarantined_chain_count=quarantined_count,
        processing_error_count=processing_error_count,
        unique_chain_identity_count=len(seen),
        violation_event_count=violation_event_count,
        violations_by_type={
            key: violations[key]
            for key in sorted(violations)
        },
    )


def _write_json_atomic(
    payload: dict[str, Any],
    output_path: str | Path,
) -> Path:
    """Write deterministic JSON atomically."""

    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    temporary = output.with_suffix(
        output.suffix + ".tmp"
    )

    text = json.dumps(
        payload,
        sort_keys=True,
        indent=2,
        ensure_ascii=True,
    )

    temporary.write_text(
        text + "\n",
        encoding="utf-8",
    )
    temporary.replace(output)

    return output


def _write_finalized_populations_atomic(
    artifacts: Iterable[GeometricValidationTaskArtifacts],
    *,
    finalized_root: str | Path,
) -> tuple[Path, Path]:
    """Stream canonical eligible/quarantined populations atomically."""

    artifacts = tuple(
        sorted(
            artifacts,
            key=lambda artifact: artifact.task_id,
        )
    )

    root = Path(finalized_root)
    root.mkdir(parents=True, exist_ok=True)

    eligible_path = root / "eligible.parquet"
    quarantined_path = root / "quarantined.parquet"

    eligible_tmp = eligible_path.with_suffix(
        eligible_path.suffix + ".tmp"
    )
    quarantined_tmp = quarantined_path.with_suffix(
        quarantined_path.suffix + ".tmp"
    )

    for temporary in (eligible_tmp, quarantined_tmp):
        if temporary.exists():
            temporary.unlink()

    eligible_writer: pq.ParquetWriter | None = None
    quarantined_writer: pq.ParquetWriter | None = None

    try:
        eligible_writer = pq.ParquetWriter(
            eligible_tmp,
            STAGE3_ELIGIBLE_CHAIN_SCHEMA,
            compression="zstd",
            version="2.6",
        )

        quarantined_writer = pq.ParquetWriter(
            quarantined_tmp,
            STAGE3_QUARANTINED_CHAIN_SCHEMA,
            compression="zstd",
            version="2.6",
        )

        for artifact in artifacts:
            (
                eligible,
                quarantined,
                _,
                _,
            ) = _validated_task_partition_rows(
                artifact
            )

            if eligible:
                eligible_writer.write_table(
                    pa.Table.from_pylist(
                        eligible,
                        schema=STAGE3_ELIGIBLE_CHAIN_SCHEMA,
                    )
                )

            if quarantined:
                quarantined_writer.write_table(
                    pa.Table.from_pylist(
                        quarantined,
                        schema=STAGE3_QUARANTINED_CHAIN_SCHEMA,
                    )
                )

    except Exception:
        if eligible_writer is not None:
            eligible_writer.close()

        if quarantined_writer is not None:
            quarantined_writer.close()

        for temporary in (eligible_tmp, quarantined_tmp):
            if temporary.exists():
                temporary.unlink()

        raise

    else:
        assert eligible_writer is not None
        assert quarantined_writer is not None

        eligible_writer.close()
        quarantined_writer.close()

        eligible_tmp.replace(eligible_path)
        quarantined_tmp.replace(quarantined_path)

    return eligible_path, quarantined_path


def build_geometric_validation_global_summary(
    artifacts: Iterable[GeometricValidationTaskArtifacts],
    *,
    snapshot: str,
    cleaning_protocol: str,
    quality_pipeline_git_commit: str,
    geometric_validation_pipeline_git_commit: str,
    finalizer_pipeline_git_commit: str,
    minimum_backbone_distance_angstrom: float,
    minimum_triangle_angle_degrees: float,
    global_validation: GeometricValidationGlobalValidation,
) -> dict[str, Any]:
    """Build deterministic Stage-2 finalization accounting."""

    artifacts = tuple(artifacts)

    relevant_sources = sum(
        _nonnegative_int(
            artifact.summary,
            "relevant_source_object_count",
            task_id=artifact.task_id,
        )
        for artifact in artifacts
    )

    downloaded_sources = sum(
        _nonnegative_int(
            artifact.summary,
            "downloaded_source_object_count",
            task_id=artifact.task_id,
        )
        for artifact in artifacts
    )

    parsed_sources = sum(
        _nonnegative_int(
            artifact.summary,
            "parsed_source_object_count",
            task_id=artifact.task_id,
        )
        for artifact in artifacts
    )

    if not (
        relevant_sources
        == downloaded_sources
        == parsed_sources
    ):
        raise GeometricValidationFinalizeError(
            "Global relevant/downloaded/parsed source accounting "
            "does not match"
        )

    result = global_validation

    return {
        "summary_schema_name": (
            GEOMETRIC_VALIDATION_GLOBAL_SUMMARY_SCHEMA_NAME
        ),
        "summary_schema_version": (
            GEOMETRIC_VALIDATION_GLOBAL_SUMMARY_SCHEMA_VERSION
        ),
        "snapshot": snapshot,
        "cleaning_protocol": cleaning_protocol,
        "quality_pipeline_git_commit": (
            quality_pipeline_git_commit
        ),
        "geometric_validation_pipeline_git_commit": (
            geometric_validation_pipeline_git_commit
        ),
        "finalizer_pipeline_git_commit": (
            finalizer_pipeline_git_commit
        ),
        "configured_minimum_backbone_distance_angstrom": (
            minimum_backbone_distance_angstrom
        ),
        "configured_minimum_triangle_angle_degrees": (
            minimum_triangle_angle_degrees
        ),
        "task_count": len(artifacts),
        "input_accepted_chain_count": (
            result.input_accepted_chain_count
        ),
        "audit_chain_count": result.audit_chain_count,
        "eligible_chain_count": result.eligible_chain_count,
        "quarantined_chain_count": (
            result.quarantined_chain_count
        ),
        "processing_error_count": (
            result.processing_error_count
        ),
        "unique_chain_identity_count": (
            result.unique_chain_identity_count
        ),
        "violation_event_count": (
            result.violation_event_count
        ),
        "violations_by_type": result.violations_by_type,
        "relevant_source_object_count": relevant_sources,
        "downloaded_source_object_count": downloaded_sources,
        "parsed_source_object_count": parsed_sources,
        "chain_accounting_valid": (
            result.input_accepted_chain_count
            == result.eligible_chain_count
            + result.quarantined_chain_count
            + result.processing_error_count
        ),
    }


def publish_geometric_validation_finalization(
    artifacts: Iterable[GeometricValidationTaskArtifacts],
    *,
    geometric_validation_root: str | Path,
    snapshot: str,
    cleaning_protocol: str,
    quality_pipeline_git_commit: str,
    geometric_validation_pipeline_git_commit: str,
    finalizer_pipeline_git_commit: str,
    minimum_backbone_distance_angstrom: float,
    minimum_triangle_angle_degrees: float,
    global_validation: GeometricValidationGlobalValidation,
) -> GeometricValidationFinalizePublication:
    """Publish canonical Stage-3 population and _SUCCESS strictly last."""

    artifacts = tuple(
        sorted(
            artifacts,
            key=lambda artifact: artifact.task_id,
        )
    )

    root = Path(geometric_validation_root)

    # Any attempted republish invalidates previous completion.
    success_path = root / "_SUCCESS"

    if success_path.exists():
        success_path.unlink()

    eligible_path, quarantined_path = (
        _write_finalized_populations_atomic(
            artifacts,
            finalized_root=root / "finalized",
        )
    )

    summary = build_geometric_validation_global_summary(
        artifacts,
        snapshot=snapshot,
        cleaning_protocol=cleaning_protocol,
        quality_pipeline_git_commit=(
            quality_pipeline_git_commit
        ),
        geometric_validation_pipeline_git_commit=(
            geometric_validation_pipeline_git_commit
        ),
        finalizer_pipeline_git_commit=(
            finalizer_pipeline_git_commit
        ),
        minimum_backbone_distance_angstrom=(
            minimum_backbone_distance_angstrom
        ),
        minimum_triangle_angle_degrees=(
            minimum_triangle_angle_degrees
        ),
        global_validation=global_validation,
    )

    if summary["chain_accounting_valid"] is not True:
        raise GeometricValidationFinalizeError(
            "Global Step-2 finalization accounting failed"
        )

    global_summary_path = _write_json_atomic(
        summary,
        root / "global_summary.json",
    )

    success_payload = {
        "success_schema_name": (
            GEOMETRIC_VALIDATION_SUCCESS_SCHEMA_NAME
        ),
        "success_schema_version": (
            GEOMETRIC_VALIDATION_SUCCESS_SCHEMA_VERSION
        ),
        "snapshot": snapshot,
        "cleaning_protocol": cleaning_protocol,
        "quality_pipeline_git_commit": (
            quality_pipeline_git_commit
        ),
        "geometric_validation_pipeline_git_commit": (
            geometric_validation_pipeline_git_commit
        ),
        "finalizer_pipeline_git_commit": (
            finalizer_pipeline_git_commit
        ),
        "task_count": len(artifacts),
        "global_summary": "global_summary.json",
        "finalized_directory": "finalized",
        "eligible_population": (
            "finalized/eligible.parquet"
        ),
        "quarantined_population": (
            "finalized/quarantined.parquet"
        ),
    }

    # Presence of _SUCCESS means all validation and publication succeeded.
    success_path = _write_json_atomic(
        success_payload,
        success_path,
    )

    return GeometricValidationFinalizePublication(
        eligible_path=eligible_path,
        quarantined_path=quarantined_path,
        global_summary_path=global_summary_path,
        success_path=success_path,
        global_summary=summary,
    )


def finalize_geometric_validation_stage(
    *,
    quality_root: str | Path,
    geometric_validation_root: str | Path,
    manifest_row_count: int,
    batch_size: int,
    snapshot: str,
    cleaning_protocol: str,
    quality_pipeline_git_commit: str,
    geometric_validation_pipeline_git_commit: str,
    finalizer_pipeline_git_commit: str,
    minimum_backbone_distance_angstrom: float,
    minimum_triangle_angle_degrees: float,
) -> GeometricValidationFinalizePublication:
    """Validate and publish one complete distributed Step-2 stage."""

    if (
        not isinstance(manifest_row_count, int)
        or isinstance(manifest_row_count, bool)
        or manifest_row_count < 1
    ):
        raise GeometricValidationFinalizeError(
            "manifest_row_count must be a positive integer"
        )

    if (
        not isinstance(batch_size, int)
        or isinstance(batch_size, bool)
        or batch_size < 1
    ):
        raise GeometricValidationFinalizeError(
            "batch_size must be a positive integer"
        )

    root = Path(geometric_validation_root)
    success_path = root / "_SUCCESS"

    if success_path.exists():
        success_path.unlink()

    task_count = (
        manifest_row_count + batch_size - 1
    ) // batch_size

    expected_task_ids = tuple(range(task_count))

    artifacts = discover_geometric_validation_task_artifacts(
        quality_root,
        root,
        expected_task_ids=expected_task_ids,
        expected_snapshot=snapshot,
        expected_cleaning_protocol=cleaning_protocol,
        expected_quality_pipeline_git_commit=(
            quality_pipeline_git_commit
        ),
        expected_geometric_validation_pipeline_git_commit=(
            geometric_validation_pipeline_git_commit
        ),
        expected_minimum_backbone_distance_angstrom=(
            minimum_backbone_distance_angstrom
        ),
        expected_minimum_triangle_angle_degrees=(
            minimum_triangle_angle_degrees
        ),
    )

    validate_geometric_validation_task_accounting(
        artifacts
    )

    global_validation = (
        validate_geometric_validation_global_state(
            artifacts
        )
    )

    return publish_geometric_validation_finalization(
        artifacts,
        geometric_validation_root=root,
        snapshot=snapshot,
        cleaning_protocol=cleaning_protocol,
        quality_pipeline_git_commit=(
            quality_pipeline_git_commit
        ),
        geometric_validation_pipeline_git_commit=(
            geometric_validation_pipeline_git_commit
        ),
        finalizer_pipeline_git_commit=(
            finalizer_pipeline_git_commit
        ),
        minimum_backbone_distance_angstrom=(
            minimum_backbone_distance_angstrom
        ),
        minimum_triangle_angle_degrees=(
            minimum_triangle_angle_degrees
        ),
        global_validation=global_validation,
    )
