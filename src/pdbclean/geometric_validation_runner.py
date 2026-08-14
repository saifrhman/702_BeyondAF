"""Source-level orchestration for post-cleaning geometric validation."""

from __future__ import annotations

from collections import Counter
from collections.abc import Callable, Iterable, Mapping
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, replace
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import resource
import time
from typing import Any

import pyarrow as pa
import pyarrow.parquet as pq

from pdbclean.geometric_validation import (
    GeometricValidationConfig,
    reconstruct_retained_backbone_chain,
    validate_post_cleaning_geometry,
)
from pdbclean.gold import (
    QUALITY_TASK_SUMMARY_SCHEMA_NAME,
    QUALITY_TASK_SUMMARY_SCHEMA_VERSION,
)
from pdbclean.mmcif_parser import (
    ChainObservation,
    MMCIFParseError,
    parse_coordinate_mmcif_bytes,
)
from pdbclean.schemas import (
    GEOMETRIC_VALIDATION_AUDIT_SCHEMA,
    GEOMETRIC_VALIDATION_PROCESSING_ERROR_SCHEMA,
)
from pdbclean.snapshot import (
    SnapshotError,
    SnapshotTransportError,
    download_verified_s3_object_bytes,
)


class GeometricValidationRunnerError(RuntimeError):
    """Raised when Step-2 orchestration cannot proceed safely."""


@dataclass(frozen=True)
class GeometricValidationTables:
    """Schema-enforced Step-2 tables for one processing batch."""

    audit: pa.Table
    processing_errors: pa.Table


def geometric_validation_records_to_tables(
    audit_records: Iterable[Mapping[str, Any]],
    processing_errors: Iterable[Mapping[str, Any]],
) -> GeometricValidationTables:
    """Convert Step-2 records into explicit Arrow schemas."""

    return GeometricValidationTables(
        audit=pa.Table.from_pylist(
            list(audit_records),
            schema=GEOMETRIC_VALIDATION_AUDIT_SCHEMA,
        ),
        processing_errors=pa.Table.from_pylist(
            list(processing_errors),
            schema=GEOMETRIC_VALIDATION_PROCESSING_ERROR_SCHEMA,
        ),
    )


def _write_geometric_validation_table_atomic(
    table: pa.Table,
    schema: pa.Schema,
    output_path: str | Path,
) -> Path:
    """Atomically write one schema-enforced Step-2 Parquet table."""

    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    canonical = table.select(schema.names).cast(schema)
    temporary = output.with_suffix(output.suffix + ".tmp")

    pq.write_table(
        canonical,
        temporary,
        compression="zstd",
        version="2.6",
    )

    temporary.replace(output)

    return output


def write_geometric_validation_shards(
    tables: GeometricValidationTables,
    output_root: str | Path,
    task_id: str | int,
) -> dict[str, Path]:
    """Write one task's Step-2 Parquet shards atomically."""

    task_id_text = str(task_id)

    if (
        not task_id_text
        or "/" in task_id_text
        or "\\" in task_id_text
        or task_id_text in {".", ".."}
    ):
        raise GeometricValidationRunnerError(
            f"Unsafe geometric-validation task_id: {task_id_text!r}"
        )

    root = Path(output_root)
    filename = f"task_{task_id_text}.parquet"

    outputs = {
        "audit": root / "audit" / filename,
        "errors": root / "errors" / filename,
    }

    _write_geometric_validation_table_atomic(
        tables.audit,
        GEOMETRIC_VALIDATION_AUDIT_SCHEMA,
        outputs["audit"],
    )

    _write_geometric_validation_table_atomic(
        tables.processing_errors,
        GEOMETRIC_VALIDATION_PROCESSING_ERROR_SCHEMA,
        outputs["errors"],
    )

    return outputs


@dataclass(frozen=True)
class SourceGeometricValidationResult:
    """Step-2 outcome for accepted chains from one source object."""

    pdb_id: str
    input_accepted_chain_count: int

    audit_records: tuple[dict[str, Any], ...] = ()
    processing_errors: tuple[dict[str, Any], ...] = ()

    source_downloaded: bool = False
    source_parsed: bool = False

    @property
    def chain_accounting_valid(self) -> bool:
        """Whether every accepted input has exactly one terminal outcome."""

        return self.input_accepted_chain_count == (
            len(self.audit_records)
            + len(self.processing_errors)
        )




@dataclass(frozen=True)
class GeometricValidationBatchResult:
    """Aggregated Step-2 outcome for one logical task."""

    manifest_source_object_count: int
    relevant_source_object_count: int
    input_accepted_chain_count: int

    downloaded_source_object_count: int
    parsed_source_object_count: int

    tables: GeometricValidationTables

    @property
    def audit_chain_count(self) -> int:
        return self.tables.audit.num_rows

    @property
    def processing_error_count(self) -> int:
        return self.tables.processing_errors.num_rows

    @property
    def chain_accounting_valid(self) -> bool:
        return self.input_accepted_chain_count == (
            self.audit_chain_count
            + self.processing_error_count
        )



GEOMETRIC_VALIDATION_TASK_SUMMARY_SCHEMA_NAME = (
    "pdbclean_geometric_validation_task_summary"
)
GEOMETRIC_VALIDATION_TASK_SUMMARY_SCHEMA_VERSION = "1.0"


def _utc_now_text() -> str:
    """Return a timezone-explicit UTC timestamp."""

    return (
        datetime.now(timezone.utc)
        .isoformat(timespec="microseconds")
        .replace("+00:00", "Z")
    )


def _linux_process_peak_memory_bytes() -> int | None:
    """Return this Linux process's peak resident memory in bytes."""

    peak_kib = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss

    if peak_kib < 0:
        return None

    return int(peak_kib) * 1024


def _slurm_environment(
    environ: Mapping[str, str] | None = None,
) -> tuple[str | None, str | None]:
    """Read optional SLURM identifiers without inventing defaults."""

    source = os.environ if environ is None else environ

    return (
        source.get("SLURM_JOB_ID"),
        source.get("SLURM_ARRAY_TASK_ID"),
    )


@dataclass(frozen=True)
class GeometricValidationTaskContext:
    """Execution/provenance context for one Step-2 task.

    Completion time, runtime, and peak memory are populated only after
    the terminal Parquet shards have been published.
    """

    task_id: str
    snapshot: str
    cleaning_protocol: str
    quality_pipeline_git_commit: str
    geometric_validation_pipeline_git_commit: str

    configured_minimum_backbone_distance_angstrom: float
    configured_minimum_triangle_angle_degrees: float

    started_at_utc: str
    completed_at_utc: str | None = None
    runtime_seconds: float | None = None

    slurm_job_id: str | None = None
    slurm_array_task_id: str | None = None
    peak_memory_bytes: int | None = None


@dataclass(frozen=True)
class GeometricValidationTaskPublication:
    """Published Step-2 shards and task completion summary."""

    shard_paths: dict[str, Path]
    summary_path: Path
    summary: dict[str, Any]

@dataclass(frozen=True)
class UpstreamQualityTaskValidation:
    """Validated Step-1 completion metadata consumed by Step 2."""

    accepted_rows: tuple[Mapping[str, Any], ...]
    quality_pipeline_git_commit: str


def validate_upstream_quality_task(
    summary: Mapping[str, Any],
    accepted_rows: Iterable[Mapping[str, Any]],
    *,
    expected_task_id: str | int,
    expected_snapshot: str,
    expected_cleaning_protocol: str,
    expected_manifest_source_object_count: int,
) -> UpstreamQualityTaskValidation:
    """Validate one completed Step-1 task before Step-2 processing.

    The Step-1 JSON summary is the task completion marker. Its accepted
    count and producer provenance must agree with the direct accepted
    Parquet shard consumed by Step 2.
    """

    task_id = str(expected_task_id)
    rows = tuple(accepted_rows)

    if (
        summary.get("summary_schema_name")
        != QUALITY_TASK_SUMMARY_SCHEMA_NAME
    ):
        raise GeometricValidationRunnerError(
            "Unexpected upstream quality-task summary schema name"
        )

    if (
        summary.get("summary_schema_version")
        != QUALITY_TASK_SUMMARY_SCHEMA_VERSION
    ):
        raise GeometricValidationRunnerError(
            "Unexpected upstream quality-task summary schema version"
        )

    if str(summary.get("task_id", "")) != task_id:
        raise GeometricValidationRunnerError(
            "Upstream quality-task summary task_id mismatch"
        )

    if summary.get("snapshot") != expected_snapshot:
        raise GeometricValidationRunnerError(
            "Upstream quality-task summary snapshot mismatch"
        )

    if (
        summary.get("cleaning_protocol")
        != expected_cleaning_protocol
    ):
        raise GeometricValidationRunnerError(
            "Upstream quality-task summary cleaning protocol mismatch"
        )

    if summary.get("source_object_accounting_valid") is not True:
        raise GeometricValidationRunnerError(
            "Upstream quality-task source accounting is not valid"
        )

    if summary.get("selected_chain_accounting_valid") is not True:
        raise GeometricValidationRunnerError(
            "Upstream quality-task selected-chain accounting is not valid"
        )

    input_source_count = summary.get(
        "input_source_object_count"
    )

    if (
        not isinstance(input_source_count, int)
        or isinstance(input_source_count, bool)
        or input_source_count
        != expected_manifest_source_object_count
    ):
        raise GeometricValidationRunnerError(
            "Upstream quality-task input source count does not match "
            "the manifest partition"
        )

    accepted_count = summary.get("accepted_chain_count")

    if (
        not isinstance(accepted_count, int)
        or isinstance(accepted_count, bool)
        or accepted_count != len(rows)
    ):
        raise GeometricValidationRunnerError(
            "Upstream quality-task accepted-chain count does not match "
            "the accepted shard"
        )

    quality_pipeline_git_commit = summary.get(
        "pipeline_git_commit"
    )

    if (
        not isinstance(quality_pipeline_git_commit, str)
        or not quality_pipeline_git_commit
    ):
        raise GeometricValidationRunnerError(
            "Upstream quality-task pipeline_git_commit must be "
            "a non-empty string"
        )

    for row in rows:
        if row.get("snapshot") != expected_snapshot:
            raise GeometricValidationRunnerError(
                "Accepted shard row snapshot does not match "
                "the upstream quality task"
            )

        if (
            row.get("cleaning_protocol")
            != expected_cleaning_protocol
        ):
            raise GeometricValidationRunnerError(
                "Accepted shard row cleaning protocol does not match "
                "the upstream quality task"
            )

        if (
            row.get("pipeline_git_commit")
            != quality_pipeline_git_commit
        ):
            raise GeometricValidationRunnerError(
                "Accepted shard row quality producer Git commit does "
                "not match the upstream quality-task summary"
            )

    return UpstreamQualityTaskValidation(
        accepted_rows=rows,
        quality_pipeline_git_commit=quality_pipeline_git_commit,
    )


def _normalized_etag(value: str) -> str:
    """Normalize quoting used around HTTP/S3 ETags."""

    return value.strip().strip('"')


def _validate_manifest_row(
    manifest_row: Mapping[str, Any],
) -> tuple[str, str, str, int, str]:
    """Validate and return immutable source-manifest lineage."""

    required = (
        "snapshot",
        "pdb_id",
        "s3_key",
        "size_bytes",
        "etag",
    )

    missing = [
        field
        for field in required
        if field not in manifest_row
    ]

    if missing:
        raise GeometricValidationRunnerError(
            "Manifest row missing required field(s): "
            + ", ".join(missing)
        )

    snapshot = manifest_row["snapshot"]
    pdb_id = manifest_row["pdb_id"]
    s3_key = manifest_row["s3_key"]
    size_bytes = manifest_row["size_bytes"]
    etag = manifest_row["etag"]

    if not isinstance(snapshot, str) or not snapshot:
        raise GeometricValidationRunnerError(
            "Manifest snapshot must be a non-empty string"
        )

    if not isinstance(pdb_id, str) or not pdb_id:
        raise GeometricValidationRunnerError(
            "Manifest pdb_id must be a non-empty string"
        )

    if not isinstance(s3_key, str) or not s3_key:
        raise GeometricValidationRunnerError(
            "Manifest s3_key must be a non-empty string"
        )

    if (
        not isinstance(size_bytes, int)
        or isinstance(size_bytes, bool)
        or size_bytes <= 0
    ):
        raise GeometricValidationRunnerError(
            "Manifest size_bytes must be a positive integer"
        )

    if (
        not isinstance(etag, str)
        or not _normalized_etag(etag)
    ):
        raise GeometricValidationRunnerError(
            "Manifest etag must be a non-empty string"
        )

    return (
        snapshot,
        pdb_id.lower(),
        s3_key,
        size_bytes,
        _normalized_etag(etag),
    )


def _accepted_lineage_issue(
    row: Mapping[str, Any],
    *,
    manifest_snapshot: str,
    manifest_pdb_id: str,
    manifest_s3_key: str,
    manifest_etag: str,
) -> str | None:
    """Return a deterministic explanation for invalid Gold lineage."""

    required = (
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
    )

    missing = [
        field
        for field in required
        if field not in row
    ]

    if missing:
        return (
            "Accepted row missing required field(s): "
            + ", ".join(missing)
        )

    if row["snapshot"] != manifest_snapshot:
        return (
            "Accepted snapshot does not match manifest: "
            f"{row['snapshot']!r} != {manifest_snapshot!r}"
        )

    pdb_id = row["pdb_id"]

    if (
        not isinstance(pdb_id, str)
        or pdb_id.lower() != manifest_pdb_id
    ):
        return (
            "Accepted pdb_id does not match manifest: "
            f"{pdb_id!r} != {manifest_pdb_id!r}"
        )

    model_id = row["model_id"]

    if (
        not isinstance(model_id, int)
        or isinstance(model_id, bool)
        or model_id <= 0
    ):
        return "Accepted model_id must be a positive integer"

    label_chain_id = row["label_chain_id"]

    if (
        not isinstance(label_chain_id, str)
        or not label_chain_id
    ):
        return "Accepted label_chain_id must be a non-empty string"

    retained_count = row["retained_residue_count"]
    retained_ids = row["retained_label_seq_ids"]

    if (
        not isinstance(retained_count, int)
        or isinstance(retained_count, bool)
        or retained_count <= 0
    ):
        return (
            "Accepted retained_residue_count must be a positive integer"
        )

    if not isinstance(retained_ids, (list, tuple)):
        return "Accepted retained_label_seq_ids must be a list or tuple"

    if not retained_ids:
        return "Accepted retained_label_seq_ids must not be empty"

    if any(
        not isinstance(residue_id, int)
        or isinstance(residue_id, bool)
        for residue_id in retained_ids
    ):
        return (
            "Accepted retained_label_seq_ids must contain only integers"
        )

    if retained_count != len(retained_ids):
        return (
            "Accepted retained_residue_count does not match "
            "retained_label_seq_ids length"
        )

    if row["source_mmcif_key"] != manifest_s3_key:
        return (
            "Accepted source_mmcif_key does not match manifest"
        )

    source_etag = row["source_etag"]

    if (
        not isinstance(source_etag, str)
        or _normalized_etag(source_etag) != manifest_etag
    ):
        return "Accepted source_etag does not match manifest"

    if (
        not isinstance(row["cleaning_protocol"], str)
        or not row["cleaning_protocol"]
    ):
        return "Accepted cleaning_protocol must be a non-empty string"

    if (
        not isinstance(row["pipeline_git_commit"], str)
        or not row["pipeline_git_commit"]
    ):
        return (
            "Accepted pipeline_git_commit must be a non-empty string"
        )

    return None


def _processing_error(
    row: Mapping[str, Any],
    *,
    processing_stage: str,
    error_type: str,
    error_message: str,
    geometric_validation_pipeline_git_commit: str,
) -> dict[str, Any]:
    """Materialize one chain-level Step-2 processing error."""

    return {
        "snapshot": row.get("snapshot"),
        "pdb_id": (
            row["pdb_id"].lower()
            if isinstance(row.get("pdb_id"), str)
            else row.get("pdb_id")
        ),
        "model_id": row.get("model_id"),
        "label_chain_id": row.get("label_chain_id"),
        "processing_stage": processing_stage,
        "error_type": error_type,
        "error_message": error_message,
        "source_mmcif_key": row.get("source_mmcif_key"),
        "source_etag": row.get("source_etag"),
        "cleaning_protocol": row.get("cleaning_protocol"),
        "quality_pipeline_git_commit": row.get(
            "pipeline_git_commit"
        ),
        "geometric_validation_pipeline_git_commit": (
            geometric_validation_pipeline_git_commit
        ),
    }


def _audit_record(
    row: Mapping[str, Any],
    *,
    validation_result: Any,
    config: GeometricValidationConfig,
    geometric_validation_pipeline_git_commit: str,
) -> dict[str, Any]:
    """Materialize one successfully reconstructed Step-2 audit row."""

    violations = validation_result.violations

    return {
        "snapshot": row["snapshot"],
        "pdb_id": row["pdb_id"].lower(),
        "model_id": row["model_id"],
        "label_chain_id": row["label_chain_id"],
        "retained_residue_count": row["retained_residue_count"],
        "retained_label_seq_ids": list(
            row["retained_label_seq_ids"]
        ),
        "source_mmcif_key": row["source_mmcif_key"],
        "source_etag": row["source_etag"],
        "cleaning_protocol": row["cleaning_protocol"],
        "quality_pipeline_git_commit": row[
            "pipeline_git_commit"
        ],
        "geometric_validation_pipeline_git_commit": (
            geometric_validation_pipeline_git_commit
        ),
        "configured_minimum_backbone_distance_angstrom": (
            config.minimum_backbone_distance_angstrom
        ),
        "configured_minimum_triangle_angle_degrees": (
            config.minimum_triangle_angle_degrees
        ),
        "passed": validation_result.passed,
        "minimum_observed_backbone_distance_angstrom": (
            validation_result
            .minimum_observed_backbone_distance_angstrom
        ),
        "minimum_observed_triangle_angle_degrees": (
            validation_result
            .minimum_observed_triangle_angle_degrees
        ),
        "minimum_observed_basis_h_norm_angstrom": (
            validation_result
            .minimum_observed_basis_h_norm_angstrom
        ),
        "violation_count": len(violations),
        "violation_types": [
            violation.violation_type
            for violation in violations
        ],
        "violation_residue_ids": [
            violation.residue_id
            for violation in violations
        ],
        "violation_details": [
            violation.details
            for violation in violations
        ],
    }


def process_geometric_validation_source(
    manifest_row: Mapping[str, Any],
    accepted_rows: Iterable[Mapping[str, Any]],
    *,
    bucket_url: str,
    config: GeometricValidationConfig,
    geometric_validation_pipeline_git_commit: str,
    timeout_seconds: int = 60,
    max_retries: int = 0,
    downloader: Callable[..., bytes] = (
        download_verified_s3_object_bytes
    ),
    parser: Callable[..., list[ChainObservation]] = (
        parse_coordinate_mmcif_bytes
    ),
) -> SourceGeometricValidationResult:
    """Validate all accepted Gold chains from one immutable source.

    The source object is downloaded and parsed at most once. Gold
    accepted-chain lineage determines exactly which chains and residues
    are reconstructed; Protocol 3.2 cleaning is never rerun here.
    """

    if (
        not isinstance(max_retries, int)
        or isinstance(max_retries, bool)
        or max_retries < 0
    ):
        raise GeometricValidationRunnerError(
            "max_retries must be a non-negative integer"
        )

    (
        manifest_snapshot,
        manifest_pdb_id,
        manifest_s3_key,
        manifest_size_bytes,
        manifest_etag,
    ) = _validate_manifest_row(manifest_row)

    rows = list(accepted_rows)

    if not rows:
        raise GeometricValidationRunnerError(
            "Step-2 source processing requires at least one "
            "accepted Gold chain"
        )

    audit_records: list[dict[str, Any]] = []
    processing_errors: list[dict[str, Any]] = []
    valid_rows: list[Mapping[str, Any]] = []

    for row in rows:
        issue = _accepted_lineage_issue(
            row,
            manifest_snapshot=manifest_snapshot,
            manifest_pdb_id=manifest_pdb_id,
            manifest_s3_key=manifest_s3_key,
            manifest_etag=manifest_etag,
        )

        if issue is None:
            valid_rows.append(row)
            continue

        processing_errors.append(
            _processing_error(
                row,
                processing_stage="accepted_lineage_validation",
                error_type="AcceptedLineageError",
                error_message=issue,
                geometric_validation_pipeline_git_commit=(
                    geometric_validation_pipeline_git_commit
                ),
            )
        )

    if not valid_rows:
        result = SourceGeometricValidationResult(
            pdb_id=manifest_pdb_id,
            input_accepted_chain_count=len(rows),
            audit_records=tuple(audit_records),
            processing_errors=tuple(processing_errors),
            source_downloaded=False,
            source_parsed=False,
        )

        if not result.chain_accounting_valid:
            raise GeometricValidationRunnerError(
                "Step-2 chain accounting failed"
            )

        return result

    compressed_bytes: bytes | None = None
    source_error: SnapshotError | None = None

    for attempt in range(max_retries + 1):
        try:
            compressed_bytes = downloader(
                bucket_url=bucket_url,
                s3_key=manifest_s3_key,
                expected_size_bytes=manifest_size_bytes,
                expected_etag=manifest_etag,
                timeout_seconds=timeout_seconds,
            )
            source_error = None
            break
        except SnapshotTransportError as exc:
            source_error = exc

            if attempt == max_retries:
                break
        except SnapshotError as exc:
            source_error = exc
            break

    if source_error is not None:
        for row in valid_rows:
            processing_errors.append(
                _processing_error(
                    row,
                    processing_stage="source_download",
                    error_type=type(source_error).__name__,
                    error_message=str(source_error),
                    geometric_validation_pipeline_git_commit=(
                        geometric_validation_pipeline_git_commit
                    ),
                )
            )

        result = SourceGeometricValidationResult(
            pdb_id=manifest_pdb_id,
            input_accepted_chain_count=len(rows),
            processing_errors=tuple(processing_errors),
            source_downloaded=False,
            source_parsed=False,
        )

        if not result.chain_accounting_valid:
            raise GeometricValidationRunnerError(
                "Step-2 chain accounting failed"
            )

        return result

    if compressed_bytes is None:
        raise GeometricValidationRunnerError(
            "Source download completed without bytes or an error"
        )

    try:
        parsed_chains = parser(
            compressed_bytes,
            pdb_id=manifest_pdb_id,
        )
    except MMCIFParseError as exc:
        for row in valid_rows:
            processing_errors.append(
                _processing_error(
                    row,
                    processing_stage="mmcif_parse",
                    error_type=type(exc).__name__,
                    error_message=str(exc),
                    geometric_validation_pipeline_git_commit=(
                        geometric_validation_pipeline_git_commit
                    ),
                )
            )

        result = SourceGeometricValidationResult(
            pdb_id=manifest_pdb_id,
            input_accepted_chain_count=len(rows),
            processing_errors=tuple(processing_errors),
            source_downloaded=True,
            source_parsed=False,
        )

        if not result.chain_accounting_valid:
            raise GeometricValidationRunnerError(
                "Step-2 chain accounting failed"
            )

        return result

    chain_index: dict[
        tuple[int, str],
        ChainObservation,
    ] = {}

    for chain in parsed_chains:
        key = (chain.model_id, chain.label_chain_id)

        if key in chain_index:
            raise GeometricValidationRunnerError(
                "Parsed source contains duplicate chain identity "
                f"{key!r}"
            )

        chain_index[key] = chain

    for row in valid_rows:
        key = (
            row["model_id"],
            row["label_chain_id"],
        )

        chain = chain_index.get(key)

        if chain is None:
            processing_errors.append(
                _processing_error(
                    row,
                    processing_stage="source_chain_lookup",
                    error_type="SourceChainNotFoundError",
                    error_message=(
                        "Accepted Gold chain identity was not found "
                        f"in parsed source: {key!r}"
                    ),
                    geometric_validation_pipeline_git_commit=(
                        geometric_validation_pipeline_git_commit
                    ),
                )
            )
            continue

        try:
            retained = reconstruct_retained_backbone_chain(
                chain,
                row["retained_label_seq_ids"],
            )
        except ValueError as exc:
            processing_errors.append(
                _processing_error(
                    row,
                    processing_stage="retained_chain_reconstruction",
                    error_type=type(exc).__name__,
                    error_message=str(exc),
                    geometric_validation_pipeline_git_commit=(
                        geometric_validation_pipeline_git_commit
                    ),
                )
            )
            continue

        try:
            validation_result = validate_post_cleaning_geometry(
                retained,
                config=config,
            )
        except ValueError as exc:
            processing_errors.append(
                _processing_error(
                    row,
                    processing_stage="geometric_validation",
                    error_type=type(exc).__name__,
                    error_message=str(exc),
                    geometric_validation_pipeline_git_commit=(
                        geometric_validation_pipeline_git_commit
                    ),
                )
            )
            continue

        audit_records.append(
            _audit_record(
                row,
                validation_result=validation_result,
                config=config,
                geometric_validation_pipeline_git_commit=(
                    geometric_validation_pipeline_git_commit
                ),
            )
        )

    result = SourceGeometricValidationResult(
        pdb_id=manifest_pdb_id,
        input_accepted_chain_count=len(rows),
        audit_records=tuple(audit_records),
        processing_errors=tuple(processing_errors),
        source_downloaded=True,
        source_parsed=True,
    )

    if not result.chain_accounting_valid:
        raise GeometricValidationRunnerError(
            "Step-2 chain accounting failed"
        )

    return result


def process_geometric_validation_batch(
    manifest_rows: Iterable[Mapping[str, Any]],
    accepted_rows: Iterable[Mapping[str, Any]],
    *,
    bucket_url: str,
    config: GeometricValidationConfig,
    geometric_validation_pipeline_git_commit: str,
    timeout_seconds: int = 60,
    max_retries: int = 0,
    download_concurrency: int = 1,
    source_processor: Callable[..., SourceGeometricValidationResult] = (
        process_geometric_validation_source
    ),
) -> GeometricValidationBatchResult:
    """Process one manifest partition and its matching accepted shard.

    Only manifest sources that have accepted Gold chains are downloaded.
    Manifest order defines deterministic source aggregation order.
    """

    if (
        not isinstance(download_concurrency, int)
        or isinstance(download_concurrency, bool)
        or download_concurrency <= 0
    ):
        raise GeometricValidationRunnerError(
            "download_concurrency must be a positive integer"
        )

    manifest = list(manifest_rows)
    accepted = list(accepted_rows)

    manifest_by_pdb: dict[str, Mapping[str, Any]] = {}
    manifest_order: list[str] = []

    for row in manifest:
        (
            _snapshot,
            pdb_id,
            _s3_key,
            _size_bytes,
            _etag,
        ) = _validate_manifest_row(row)

        if pdb_id in manifest_by_pdb:
            raise GeometricValidationRunnerError(
                "Manifest partition contains duplicate pdb_id: "
                f"{pdb_id!r}"
            )

        manifest_by_pdb[pdb_id] = row
        manifest_order.append(pdb_id)

    accepted_by_pdb: dict[
        str,
        list[Mapping[str, Any]],
    ] = {}

    seen_chain_identities: set[tuple[str, int, str]] = set()

    for row in accepted:
        pdb_id_raw = row.get("pdb_id")

        if not isinstance(pdb_id_raw, str) or not pdb_id_raw:
            raise GeometricValidationRunnerError(
                "Accepted task row requires a non-empty pdb_id"
            )

        pdb_id = pdb_id_raw.lower()

        if pdb_id not in manifest_by_pdb:
            raise GeometricValidationRunnerError(
                "Accepted Gold source is outside the matching "
                f"manifest partition: {pdb_id!r}"
            )

        model_id = row.get("model_id")
        label_chain_id = row.get("label_chain_id")

        if (
            not isinstance(model_id, int)
            or isinstance(model_id, bool)
            or not isinstance(label_chain_id, str)
            or not label_chain_id
        ):
            raise GeometricValidationRunnerError(
                "Accepted task row has invalid chain identity"
            )

        identity = (
            pdb_id,
            model_id,
            label_chain_id,
        )

        if identity in seen_chain_identities:
            raise GeometricValidationRunnerError(
                "Accepted task shard contains duplicate chain identity: "
                f"{identity!r}"
            )

        seen_chain_identities.add(identity)

        accepted_by_pdb.setdefault(
            pdb_id,
            [],
        ).append(row)

    relevant_pdb_ids = [
        pdb_id
        for pdb_id in manifest_order
        if pdb_id in accepted_by_pdb
    ]

    audit_records: list[dict[str, Any]] = []
    processing_errors: list[dict[str, Any]] = []

    downloaded_source_object_count = 0
    parsed_source_object_count = 0

    def process_pdb(
        pdb_id: str,
    ) -> SourceGeometricValidationResult:
        return source_processor(
            manifest_by_pdb[pdb_id],
            accepted_by_pdb[pdb_id],
            bucket_url=bucket_url,
            config=config,
            geometric_validation_pipeline_git_commit=(
                geometric_validation_pipeline_git_commit
            ),
            timeout_seconds=timeout_seconds,
            max_retries=max_retries,
        )

    # executor.map preserves input ordering, so task outputs remain
    # deterministic even if source downloads finish out of order.
    with ThreadPoolExecutor(
        max_workers=download_concurrency
    ) as executor:
        results = executor.map(
            process_pdb,
            relevant_pdb_ids,
        )

        for result in results:
            if not result.chain_accounting_valid:
                raise GeometricValidationRunnerError(
                    "Source-level Step-2 chain accounting failed for "
                    f"{result.pdb_id!r}"
                )

            downloaded_source_object_count += int(
                result.source_downloaded
            )
            parsed_source_object_count += int(
                result.source_parsed
            )

            audit_records.extend(result.audit_records)
            processing_errors.extend(result.processing_errors)

    tables = geometric_validation_records_to_tables(
        audit_records,
        processing_errors,
    )

    batch = GeometricValidationBatchResult(
        manifest_source_object_count=len(manifest),
        relevant_source_object_count=len(relevant_pdb_ids),
        input_accepted_chain_count=len(accepted),
        downloaded_source_object_count=(
            downloaded_source_object_count
        ),
        parsed_source_object_count=(
            parsed_source_object_count
        ),
        tables=tables,
    )

    if not batch.chain_accounting_valid:
        raise GeometricValidationRunnerError(
            "Task-level Step-2 chain accounting failed"
        )

    return batch


def _sorted_counts(values: Iterable[str]) -> dict[str, int]:
    """Return deterministic lexical-key counts."""

    counts = Counter(values)

    return {
        key: counts[key]
        for key in sorted(counts)
    }


def build_geometric_validation_task_summary(
    batch: GeometricValidationBatchResult,
    context: GeometricValidationTaskContext,
) -> dict[str, Any]:
    """Build one deterministic Step-2 task summary."""

    if context.completed_at_utc is None:
        raise GeometricValidationRunnerError(
            "Step-2 task summary requires completed_at_utc"
        )

    if context.runtime_seconds is None:
        raise GeometricValidationRunnerError(
            "Step-2 task summary requires runtime_seconds"
        )

    if context.runtime_seconds < 0:
        raise GeometricValidationRunnerError(
            "Step-2 task runtime cannot be negative"
        )

    audit_rows = batch.tables.audit.to_pylist()
    error_rows = batch.tables.processing_errors.to_pylist()

    passed_count = sum(
        bool(row["passed"])
        for row in audit_rows
    )
    violated_count = len(audit_rows) - passed_count

    violation_types = [
        violation_type
        for row in audit_rows
        for violation_type in row["violation_types"]
    ]

    chain_accounting_valid = (
        batch.input_accepted_chain_count
        == len(audit_rows) + len(error_rows)
    )

    return {
        "summary_schema_name": (
            GEOMETRIC_VALIDATION_TASK_SUMMARY_SCHEMA_NAME
        ),
        "summary_schema_version": (
            GEOMETRIC_VALIDATION_TASK_SUMMARY_SCHEMA_VERSION
        ),
        "task_id": context.task_id,
        "snapshot": context.snapshot,
        "cleaning_protocol": context.cleaning_protocol,
        "quality_pipeline_git_commit": (
            context.quality_pipeline_git_commit
        ),
        "geometric_validation_pipeline_git_commit": (
            context.geometric_validation_pipeline_git_commit
        ),
        "configured_minimum_backbone_distance_angstrom": (
            context.configured_minimum_backbone_distance_angstrom
        ),
        "configured_minimum_triangle_angle_degrees": (
            context.configured_minimum_triangle_angle_degrees
        ),
        "started_at_utc": context.started_at_utc,
        "completed_at_utc": context.completed_at_utc,
        "runtime_seconds": context.runtime_seconds,
        "slurm_job_id": context.slurm_job_id,
        "slurm_array_task_id": context.slurm_array_task_id,
        "peak_memory_bytes": context.peak_memory_bytes,
        "manifest_source_object_count": (
            batch.manifest_source_object_count
        ),
        "relevant_source_object_count": (
            batch.relevant_source_object_count
        ),
        "downloaded_source_object_count": (
            batch.downloaded_source_object_count
        ),
        "parsed_source_object_count": (
            batch.parsed_source_object_count
        ),
        "input_accepted_chain_count": (
            batch.input_accepted_chain_count
        ),
        "audit_chain_count": len(audit_rows),
        "geometric_passed_chain_count": passed_count,
        "geometric_violated_chain_count": violated_count,
        "processing_error_count": len(error_rows),
        "violations_by_type": _sorted_counts(
            violation_types
        ),
        "processing_errors_by_stage": _sorted_counts(
            row["processing_stage"]
            for row in error_rows
        ),
        "processing_errors_by_type": _sorted_counts(
            row["error_type"]
            for row in error_rows
        ),
        "chain_accounting_valid": chain_accounting_valid,
    }


def write_geometric_validation_task_summary_atomic(
    summary: dict[str, Any],
    output_root: str | Path,
) -> Path:
    """Validate and atomically write one Step-2 task summary."""

    if (
        summary.get("summary_schema_name")
        != GEOMETRIC_VALIDATION_TASK_SUMMARY_SCHEMA_NAME
    ):
        raise ValueError(
            "Unexpected geometric-validation task-summary schema name"
        )

    if (
        summary.get("summary_schema_version")
        != GEOMETRIC_VALIDATION_TASK_SUMMARY_SCHEMA_VERSION
    ):
        raise ValueError(
            "Unexpected geometric-validation task-summary schema version"
        )

    if summary.get("chain_accounting_valid") is not True:
        raise ValueError(
            "Cannot publish geometric-validation task summary: "
            "chain accounting failed"
        )

    task_id = str(summary.get("task_id", ""))

    if (
        not task_id
        or "/" in task_id
        or "\\" in task_id
        or task_id in {".", ".."}
    ):
        raise ValueError(
            f"Unsafe geometric-validation task_id: {task_id!r}"
        )

    output = (
        Path(output_root)
        / "summaries"
        / f"task_{task_id}.json"
    )

    output.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    temporary = output.with_suffix(
        output.suffix + ".tmp"
    )

    payload = json.dumps(
        summary,
        sort_keys=True,
        indent=2,
        ensure_ascii=True,
    )

    temporary.write_text(
        payload + "\n",
        encoding="utf-8",
    )

    temporary.replace(output)

    return output


def publish_geometric_validation_batch(
    batch: GeometricValidationBatchResult,
    *,
    output_root: str | Path,
    context: GeometricValidationTaskContext,
    started_perf_counter: float,
    shard_writer: Callable[..., dict[str, Path]] = (
        write_geometric_validation_shards
    ),
    summary_writer: Callable[..., Path] = (
        write_geometric_validation_task_summary_atomic
    ),
    utc_now: Callable[[], str] = _utc_now_text,
    perf_counter: Callable[[], float] = time.perf_counter,
    peak_memory_reader: Callable[
        [], int | None
    ] = _linux_process_peak_memory_bytes,
) -> GeometricValidationTaskPublication:
    """Publish Step-2 shards, capture completion, then publish summary."""

    if not batch.chain_accounting_valid:
        raise GeometricValidationRunnerError(
            "Cannot publish Step-2 task: chain accounting failed"
        )

    if (
        not isinstance(started_perf_counter, (int, float))
        or isinstance(started_perf_counter, bool)
    ):
        raise GeometricValidationRunnerError(
            "started_perf_counter must be numeric"
        )

    # Terminal Parquet outcomes are published before completion metadata
    # is captured.
    shard_paths = shard_writer(
        batch.tables,
        output_root,
        context.task_id,
    )

    completed_at_utc = utc_now()
    runtime_seconds = perf_counter() - started_perf_counter

    if runtime_seconds < 0:
        raise GeometricValidationRunnerError(
            "Monotonic Step-2 task runtime cannot be negative"
        )

    completed_context = replace(
        context,
        completed_at_utc=completed_at_utc,
        runtime_seconds=runtime_seconds,
        peak_memory_bytes=peak_memory_reader(),
    )

    summary = build_geometric_validation_task_summary(
        batch,
        completed_context,
    )

    # JSON summary is deliberately published last and acts as the
    # task-level completion marker.
    summary_path = summary_writer(
        summary,
        output_root,
    )

    return GeometricValidationTaskPublication(
        shard_paths=shard_paths,
        summary_path=summary_path,
        summary=summary,
    )


def execute_geometric_validation_task(
    manifest_rows: Iterable[Mapping[str, Any]],
    accepted_rows: Iterable[Mapping[str, Any]],
    *,
    output_root: str | Path,
    task_id: str | int,
    snapshot: str,
    bucket_url: str,
    config: GeometricValidationConfig,
    cleaning_protocol: str,
    quality_pipeline_git_commit: str,
    geometric_validation_pipeline_git_commit: str,
    timeout_seconds: int = 60,
    max_retries: int = 0,
    download_concurrency: int = 1,
    environ: Mapping[str, str] | None = None,
    batch_processor: Callable[
        ...,
        GeometricValidationBatchResult,
    ] = process_geometric_validation_batch,
    publisher: Callable[
        ...,
        GeometricValidationTaskPublication,
    ] = publish_geometric_validation_batch,
    utc_now: Callable[[], str] = _utc_now_text,
    perf_counter: Callable[[], float] = time.perf_counter,
) -> GeometricValidationTaskPublication:
    """Execute one complete post-cleaning geometric-validation task."""

    task_id_text = str(task_id)

    if (
        not task_id_text
        or "/" in task_id_text
        or "\\" in task_id_text
        or task_id_text in {".", ".."}
    ):
        raise GeometricValidationRunnerError(
            f"Unsafe geometric-validation task_id: {task_id_text!r}"
        )

    slurm_job_id, slurm_array_task_id = _slurm_environment(
        environ
    )

    # Capture task start immediately before scientific source processing.
    started_at_utc = utc_now()
    started_perf_counter = perf_counter()

    batch = batch_processor(
        manifest_rows,
        accepted_rows,
        bucket_url=bucket_url,
        config=config,
        geometric_validation_pipeline_git_commit=(
            geometric_validation_pipeline_git_commit
        ),
        timeout_seconds=timeout_seconds,
        max_retries=max_retries,
        download_concurrency=download_concurrency,
    )

    context = GeometricValidationTaskContext(
        task_id=task_id_text,
        snapshot=snapshot,
        cleaning_protocol=cleaning_protocol,
        quality_pipeline_git_commit=quality_pipeline_git_commit,
        geometric_validation_pipeline_git_commit=(
            geometric_validation_pipeline_git_commit
        ),
        configured_minimum_backbone_distance_angstrom=(
            config.minimum_backbone_distance_angstrom
        ),
        configured_minimum_triangle_angle_degrees=(
            config.minimum_triangle_angle_degrees
        ),
        started_at_utc=started_at_utc,
        slurm_job_id=slurm_job_id,
        slurm_array_task_id=slurm_array_task_id,
    )

    # The same clocks are forwarded so publication captures completion
    # only after the terminal Parquet shards have been written.
    return publisher(
        batch,
        output_root=output_root,
        context=context,
        started_perf_counter=started_perf_counter,
        utc_now=utc_now,
        perf_counter=perf_counter,
    )
