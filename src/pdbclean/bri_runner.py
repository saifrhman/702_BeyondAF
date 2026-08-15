"""Production orchestration for Stage-3 Definition 3.4 BRI.

This module consumes only the canonical population published by the
completed post-cleaning geometric-validation stage.

The scientific BRI calculation itself lives in ``pdbclean.bri``.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, replace
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import resource
import time
from typing import Any, Callable, Iterable, Mapping

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

from pdbclean.bri import compute_bri
from pdbclean.geometric_validation import (
    reconstruct_retained_backbone_chain,
)
from pdbclean.geometric_validation_finalize import (
    GEOMETRIC_VALIDATION_GLOBAL_SUMMARY_SCHEMA_NAME,
    GEOMETRIC_VALIDATION_GLOBAL_SUMMARY_SCHEMA_VERSION,
    GEOMETRIC_VALIDATION_SUCCESS_SCHEMA_NAME,
    GEOMETRIC_VALIDATION_SUCCESS_SCHEMA_VERSION,
)
from pdbclean.mmcif_parser import (
    ChainObservation,
    MMCIFParseError,
    parse_coordinate_mmcif_bytes,
)
from pdbclean.schemas import (
    STAGE3_BRI_CHAIN_SCHEMA,
    STAGE3_BRI_PROCESSING_ERROR_SCHEMA,
    STAGE3_ELIGIBLE_CHAIN_SCHEMA,
)
from pdbclean.snapshot import (
    SnapshotError,
    SnapshotTransportError,
    download_verified_s3_object_bytes,
)


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


@dataclass(frozen=True)
class SourceBRIResult:
    """Terminal Stage-3 outcomes for one immutable source object."""

    pdb_id: str
    input_eligible_chain_count: int

    bri_records: tuple[dict[str, Any], ...] = ()
    processing_errors: tuple[dict[str, Any], ...] = ()

    source_downloaded: bool = False
    source_parsed: bool = False

    @property
    def chain_accounting_valid(self) -> bool:
        """Every eligible chain must have exactly one terminal outcome."""

        return self.input_eligible_chain_count == (
            len(self.bri_records)
            + len(self.processing_errors)
        )


def _normalized_etag(value: str) -> str:
    """Normalize optional HTTP/S3 ETag quoting."""

    return value.strip().strip('"')


def _validate_manifest_source_row(
    manifest_row: Mapping[str, Any],
    *,
    expected_snapshot: str,
) -> tuple[str, str, int, str]:
    """Validate immutable source lineage required for BRI processing."""

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
        raise BRIRunnerError(
            "Manifest row missing required field(s): "
            + ", ".join(missing)
        )

    snapshot = manifest_row["snapshot"]
    pdb_id = manifest_row["pdb_id"]
    s3_key = manifest_row["s3_key"]
    size_bytes = manifest_row["size_bytes"]
    etag = manifest_row["etag"]

    if snapshot != expected_snapshot:
        raise BRIRunnerError(
            "Manifest source snapshot does not match validated "
            "Stage-2 completion"
        )

    if not isinstance(pdb_id, str) or not pdb_id:
        raise BRIRunnerError(
            "Manifest pdb_id must be a non-empty string"
        )

    if not isinstance(s3_key, str) or not s3_key:
        raise BRIRunnerError(
            "Manifest s3_key must be a non-empty string"
        )

    if (
        not isinstance(size_bytes, int)
        or isinstance(size_bytes, bool)
        or size_bytes <= 0
    ):
        raise BRIRunnerError(
            "Manifest size_bytes must be a positive integer"
        )

    if not isinstance(etag, str) or not _normalized_etag(etag):
        raise BRIRunnerError(
            "Manifest etag must be a non-empty string"
        )

    return (
        pdb_id.lower(),
        s3_key,
        size_bytes,
        _normalized_etag(etag),
    )


def _eligible_lineage_issue(
    row: Mapping[str, Any],
    *,
    upstream: UpstreamGeometricValidation,
    manifest_pdb_id: str,
    manifest_s3_key: str,
    manifest_etag: str,
) -> str | None:
    """Return an immutable-lineage problem for one Stage-3 input row."""

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
        raise BRIRunnerError(
            "Canonical eligible row missing required field(s): "
            + ", ".join(missing)
        )

    if row["snapshot"] != upstream.snapshot:
        return "Eligible row snapshot does not match Stage-2 completion"

    pdb_id = row["pdb_id"]

    if (
        not isinstance(pdb_id, str)
        or pdb_id.lower() != manifest_pdb_id
    ):
        return "Eligible row pdb_id does not match manifest source"

    if row["source_mmcif_key"] != manifest_s3_key:
        return (
            "Eligible row source_mmcif_key does not match "
            "manifest source"
        )

    source_etag = row["source_etag"]

    if (
        not isinstance(source_etag, str)
        or _normalized_etag(source_etag) != manifest_etag
    ):
        return "Eligible row source_etag does not match manifest source"

    if row["cleaning_protocol"] != upstream.cleaning_protocol:
        return (
            "Eligible row cleaning_protocol does not match "
            "Stage-2 completion"
        )

    if (
        row["pipeline_git_commit"]
        != upstream.quality_pipeline_git_commit
    ):
        return (
            "Eligible row quality producer does not match "
            "Stage-2 completion"
        )

    model_id = row["model_id"]

    if (
        not isinstance(model_id, int)
        or isinstance(model_id, bool)
        or model_id <= 0
    ):
        return "Eligible row model_id must be a positive integer"

    label_chain_id = row["label_chain_id"]

    if not isinstance(label_chain_id, str) or not label_chain_id:
        return "Eligible row label_chain_id must be non-empty"

    retained_count = row["retained_residue_count"]
    retained_ids = row["retained_label_seq_ids"]

    if (
        not isinstance(retained_count, int)
        or isinstance(retained_count, bool)
        or retained_count < 1
    ):
        return (
            "Eligible row retained_residue_count must be "
            "a positive integer"
        )

    if not isinstance(retained_ids, (list, tuple)):
        return "Eligible row retained_label_seq_ids must be a sequence"

    if len(retained_ids) != retained_count:
        return (
            "Eligible row retained_residue_count does not match "
            "retained_label_seq_ids"
        )

    return None


def _bri_processing_error(
    row: Mapping[str, Any],
    *,
    processing_stage: str,
    error_type: str,
    error_message: str,
    upstream: UpstreamGeometricValidation,
    bri_pipeline_git_commit: str,
) -> dict[str, Any]:
    """Materialize one terminal Stage-3 processing error."""

    return {
        "snapshot": row["snapshot"],
        "pdb_id": row["pdb_id"],
        "model_id": row["model_id"],
        "label_chain_id": row["label_chain_id"],
        "retained_residue_count": row[
            "retained_residue_count"
        ],
        "retained_label_seq_ids": list(
            row["retained_label_seq_ids"]
        ),
        "processing_stage": processing_stage,
        "error_type": error_type,
        "error_message": error_message,
        "source_mmcif_key": row["source_mmcif_key"],
        "source_etag": row["source_etag"],
        "cleaning_protocol": row["cleaning_protocol"],
        "quality_pipeline_git_commit": (
            upstream.quality_pipeline_git_commit
        ),
        "geometric_validation_pipeline_git_commit": (
            upstream.geometric_validation_pipeline_git_commit
        ),
        "geometric_validation_finalizer_git_commit": (
            upstream.geometric_validation_finalizer_git_commit
        ),
        "bri_pipeline_git_commit": bri_pipeline_git_commit,
    }


def _validated_bri_matrix(
    value: object,
    *,
    expected_residue_count: int,
) -> np.ndarray:
    """Enforce the frozen Definition-3.4 production output contract."""

    if not isinstance(value, np.ndarray):
        raise BRIRunnerError(
            "compute_bri() did not return a NumPy array"
        )

    if value.dtype != np.float64:
        raise BRIRunnerError(
            "compute_bri() did not return float64 coordinates"
        )

    if value.shape != (expected_residue_count, 9):
        raise BRIRunnerError(
            "compute_bri() returned an unexpected matrix shape"
        )

    if not np.isfinite(value).all():
        raise BRIRunnerError(
            "compute_bri() returned non-finite coordinates"
        )

    if not np.array_equal(value, np.around(value, 3)):
        raise BRIRunnerError(
            "compute_bri() returned non-canonical coordinates"
        )

    return value


def _bri_record(
    row: Mapping[str, Any],
    *,
    bri: np.ndarray,
    upstream: UpstreamGeometricValidation,
    bri_pipeline_git_commit: str,
) -> dict[str, Any]:
    """Attach canonical BRI and explicit producer provenance."""

    record = {
        field.name: row[field.name]
        for field in STAGE3_ELIGIBLE_CHAIN_SCHEMA
        if field.name != "pipeline_git_commit"
    }

    record.update(
        {
            "quality_pipeline_git_commit": (
                upstream.quality_pipeline_git_commit
            ),
            "geometric_validation_pipeline_git_commit": (
                upstream.geometric_validation_pipeline_git_commit
            ),
            "geometric_validation_finalizer_git_commit": (
                upstream.geometric_validation_finalizer_git_commit
            ),
            "bri_pipeline_git_commit": (
                bri_pipeline_git_commit
            ),
            "bri": bri.tolist(),
        }
    )

    return record


def process_bri_source(
    manifest_row: Mapping[str, Any],
    eligible_rows: Iterable[Mapping[str, Any]],
    *,
    upstream: UpstreamGeometricValidation,
    bucket_url: str,
    bri_pipeline_git_commit: str,
    timeout_seconds: int = 60,
    max_retries: int = 0,
    downloader: Callable[..., bytes] = (
        download_verified_s3_object_bytes
    ),
    parser: Callable[..., list[ChainObservation]] = (
        parse_coordinate_mmcif_bytes
    ),
    bri_computer: Callable[[ChainObservation], np.ndarray] = compute_bri,
) -> SourceBRIResult:
    """Compute BRI for eligible chains belonging to one source object."""

    if (
        not isinstance(max_retries, int)
        or isinstance(max_retries, bool)
        or max_retries < 0
    ):
        raise BRIRunnerError(
            "max_retries must be a non-negative integer"
        )

    bri_pipeline_git_commit = _validated_git_commit(
        bri_pipeline_git_commit,
        field="bri_pipeline_git_commit",
    )

    (
        manifest_pdb_id,
        manifest_s3_key,
        manifest_size_bytes,
        manifest_etag,
    ) = _validate_manifest_source_row(
        manifest_row,
        expected_snapshot=upstream.snapshot,
    )

    rows = tuple(eligible_rows)

    if not rows:
        return SourceBRIResult(
            pdb_id=manifest_pdb_id,
            input_eligible_chain_count=0,
        )

    bri_records: list[dict[str, Any]] = []
    processing_errors: list[dict[str, Any]] = []
    valid_rows: list[Mapping[str, Any]] = []

    for row in rows:
        issue = _eligible_lineage_issue(
            row,
            upstream=upstream,
            manifest_pdb_id=manifest_pdb_id,
            manifest_s3_key=manifest_s3_key,
            manifest_etag=manifest_etag,
        )

        if issue is None:
            valid_rows.append(row)
            continue

        processing_errors.append(
            _bri_processing_error(
                row,
                processing_stage="eligible_lineage_validation",
                error_type="EligibleLineageError",
                error_message=issue,
                upstream=upstream,
                bri_pipeline_git_commit=bri_pipeline_git_commit,
            )
        )

    if not valid_rows:
        result = SourceBRIResult(
            pdb_id=manifest_pdb_id,
            input_eligible_chain_count=len(rows),
            bri_records=tuple(bri_records),
            processing_errors=tuple(processing_errors),
        )

        if not result.chain_accounting_valid:
            raise BRIRunnerError(
                "Stage-3 source chain accounting failed"
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
                _bri_processing_error(
                    row,
                    processing_stage="source_download",
                    error_type=type(source_error).__name__,
                    error_message=str(source_error),
                    upstream=upstream,
                    bri_pipeline_git_commit=bri_pipeline_git_commit,
                )
            )

        result = SourceBRIResult(
            pdb_id=manifest_pdb_id,
            input_eligible_chain_count=len(rows),
            processing_errors=tuple(processing_errors),
            source_downloaded=False,
            source_parsed=False,
        )

        if not result.chain_accounting_valid:
            raise BRIRunnerError(
                "Stage-3 source chain accounting failed"
            )

        return result

    if compressed_bytes is None:
        raise BRIRunnerError(
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
                _bri_processing_error(
                    row,
                    processing_stage="mmcif_parse",
                    error_type=type(exc).__name__,
                    error_message=str(exc),
                    upstream=upstream,
                    bri_pipeline_git_commit=bri_pipeline_git_commit,
                )
            )

        result = SourceBRIResult(
            pdb_id=manifest_pdb_id,
            input_eligible_chain_count=len(rows),
            processing_errors=tuple(processing_errors),
            source_downloaded=True,
            source_parsed=False,
        )

        if not result.chain_accounting_valid:
            raise BRIRunnerError(
                "Stage-3 source chain accounting failed"
            )

        return result

    chain_index: dict[
        tuple[int, str],
        ChainObservation,
    ] = {}

    for chain in parsed_chains:
        key = (
            chain.model_id,
            chain.label_chain_id,
        )

        if key in chain_index:
            raise BRIRunnerError(
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
                _bri_processing_error(
                    row,
                    processing_stage="source_chain_lookup",
                    error_type="SourceChainNotFoundError",
                    error_message=(
                        "Canonical Stage-3 chain identity was not found "
                        f"in parsed source: {key!r}"
                    ),
                    upstream=upstream,
                    bri_pipeline_git_commit=bri_pipeline_git_commit,
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
                _bri_processing_error(
                    row,
                    processing_stage="retained_chain_reconstruction",
                    error_type=type(exc).__name__,
                    error_message=str(exc),
                    upstream=upstream,
                    bri_pipeline_git_commit=bri_pipeline_git_commit,
                )
            )
            continue

        try:
            bri = bri_computer(retained)
        except ValueError as exc:
            processing_errors.append(
                _bri_processing_error(
                    row,
                    processing_stage="bri_computation",
                    error_type=type(exc).__name__,
                    error_message=str(exc),
                    upstream=upstream,
                    bri_pipeline_git_commit=bri_pipeline_git_commit,
                )
            )
            continue

        bri = _validated_bri_matrix(
            bri,
            expected_residue_count=row[
                "retained_residue_count"
            ],
        )

        bri_records.append(
            _bri_record(
                row,
                bri=bri,
                upstream=upstream,
                bri_pipeline_git_commit=bri_pipeline_git_commit,
            )
        )

    result = SourceBRIResult(
        pdb_id=manifest_pdb_id,
        input_eligible_chain_count=len(rows),
        bri_records=tuple(bri_records),
        processing_errors=tuple(processing_errors),
        source_downloaded=True,
        source_parsed=True,
    )

    if not result.chain_accounting_valid:
        raise BRIRunnerError(
            "Stage-3 source chain accounting failed"
        )

    return result


BRI_TASK_SUMMARY_SCHEMA_NAME = "pdbclean_stage3_bri_task_summary"
BRI_TASK_SUMMARY_SCHEMA_VERSION = "1.0"


@dataclass(frozen=True)
class BRITables:
    """Schema-enforced terminal Stage-3 task tables."""

    chains: pa.Table
    processing_errors: pa.Table


@dataclass(frozen=True)
class BRIBatchResult:
    """Aggregated outcome for one logical Stage-3 manifest task."""

    manifest_source_object_count: int
    relevant_source_object_count: int

    input_eligible_chain_count: int

    downloaded_source_object_count: int
    parsed_source_object_count: int

    tables: BRITables

    @property
    def chain_accounting_valid(self) -> bool:
        """Every canonical eligible chain has one terminal outcome."""

        return self.input_eligible_chain_count == (
            self.tables.chains.num_rows
            + self.tables.processing_errors.num_rows
        )


@dataclass(frozen=True)
class BRITaskContext:
    """Execution and provenance context for one Stage-3 task."""

    task_id: str
    snapshot: str
    cleaning_protocol: str

    quality_pipeline_git_commit: str
    geometric_validation_pipeline_git_commit: str
    geometric_validation_finalizer_git_commit: str
    bri_pipeline_git_commit: str

    started_at_utc: str
    completed_at_utc: str | None = None
    runtime_seconds: float | None = None

    slurm_job_id: str | None = None
    slurm_array_task_id: str | None = None
    peak_memory_bytes: int | None = None


@dataclass(frozen=True)
class BRITaskPublication:
    """Published Stage-3 task shards and completion summary."""

    shard_paths: dict[str, Path]
    summary_path: Path
    summary: dict[str, Any]


def _safe_bri_task_id(task_id: str | int) -> str:
    """Validate a task identifier before filesystem use."""

    value = str(task_id)

    if (
        not value
        or "/" in value
        or "\\" in value
        or value in {".", ".."}
    ):
        raise BRIRunnerError(
            f"Unsafe Stage-3 BRI task_id: {value!r}"
        )

    return value


def _utc_now_text() -> str:
    """Return an explicit UTC timestamp for task provenance."""

    return (
        datetime.now(timezone.utc)
        .isoformat(timespec="microseconds")
        .replace("+00:00", "Z")
    )


def _linux_process_peak_memory_bytes() -> int | None:
    """Return process peak resident memory on Linux."""

    peak_kib = resource.getrusage(
        resource.RUSAGE_SELF
    ).ru_maxrss

    if peak_kib < 0:
        return None

    return int(peak_kib) * 1024


def _slurm_environment(
    environ: Mapping[str, str] | None = None,
) -> tuple[str | None, str | None]:
    """Read optional Slurm identifiers without inventing defaults."""

    source = os.environ if environ is None else environ

    return (
        source.get("SLURM_JOB_ID"),
        source.get("SLURM_ARRAY_TASK_ID"),
    )


def _eligible_chain_identity(
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
        raise BRIRunnerError(
            "Canonical eligible row is missing identity lineage"
        ) from exc


def _materialize_bri_tables(
    bri_records: Iterable[Mapping[str, Any]],
    processing_errors: Iterable[Mapping[str, Any]],
) -> BRITables:
    """Materialize terminal records against frozen Stage-3 schemas."""

    try:
        chains = pa.Table.from_pylist(
            list(bri_records),
            schema=STAGE3_BRI_CHAIN_SCHEMA,
        )

        errors = pa.Table.from_pylist(
            list(processing_errors),
            schema=STAGE3_BRI_PROCESSING_ERROR_SCHEMA,
        )
    except (pa.ArrowException, TypeError, ValueError) as exc:
        raise BRIRunnerError(
            f"Cannot materialize Stage-3 BRI tables: {exc}"
        ) from exc

    if not chains.schema.equals(
        STAGE3_BRI_CHAIN_SCHEMA,
        check_metadata=True,
    ):
        raise BRIRunnerError(
            "Materialized BRI-chain table schema mismatch"
        )

    if not errors.schema.equals(
        STAGE3_BRI_PROCESSING_ERROR_SCHEMA,
        check_metadata=True,
    ):
        raise BRIRunnerError(
            "Materialized BRI processing-error table schema mismatch"
        )

    return BRITables(
        chains=chains,
        processing_errors=errors,
    )


def process_bri_batch(
    manifest_rows: Iterable[Mapping[str, Any]],
    eligible_rows: Iterable[Mapping[str, Any]],
    *,
    upstream: UpstreamGeometricValidation,
    bucket_url: str,
    bri_pipeline_git_commit: str,
    timeout_seconds: int = 60,
    max_retries: int = 0,
    source_processor: Callable[..., SourceBRIResult] = process_bri_source,
) -> BRIBatchResult:
    """Process one deterministic manifest partition for Stage-3 BRI."""

    manifest = tuple(manifest_rows)
    eligible = tuple(eligible_rows)

    bri_pipeline_git_commit = _validated_git_commit(
        bri_pipeline_git_commit,
        field="bri_pipeline_git_commit",
    )

    manifest_by_source: dict[
        str,
        Mapping[str, Any],
    ] = {}

    manifest_source_order: list[str] = []

    for manifest_row in manifest:
        (
            _,
            source_key,
            _,
            _,
        ) = _validate_manifest_source_row(
            manifest_row,
            expected_snapshot=upstream.snapshot,
        )

        if source_key in manifest_by_source:
            raise BRIRunnerError(
                "Manifest task contains duplicate source_mmcif key: "
                f"{source_key!r}"
            )

        manifest_by_source[source_key] = manifest_row
        manifest_source_order.append(source_key)

    eligible_by_source: dict[
        str,
        list[Mapping[str, Any]],
    ] = {}

    seen_identities: set[
        tuple[str, str, int, str]
    ] = set()

    for row in eligible:
        identity = _eligible_chain_identity(row)

        if identity in seen_identities:
            raise BRIRunnerError(
                "Stage-3 task contains duplicate eligible chain "
                f"identity: {identity!r}"
            )

        seen_identities.add(identity)

        source_key = row.get("source_mmcif_key")

        if (
            not isinstance(source_key, str)
            or not source_key
        ):
            raise BRIRunnerError(
                "Canonical eligible row has invalid source_mmcif_key"
            )

        if source_key not in manifest_by_source:
            raise BRIRunnerError(
                "Canonical eligible row belongs to a source outside "
                "the current manifest partition"
            )

        eligible_by_source.setdefault(
            source_key,
            [],
        ).append(row)

    bri_records: list[dict[str, Any]] = []
    processing_errors: list[dict[str, Any]] = []

    relevant_source_object_count = 0
    downloaded_source_object_count = 0
    parsed_source_object_count = 0

    for source_key in manifest_source_order:
        source_rows = eligible_by_source.get(source_key)

        if not source_rows:
            continue

        relevant_source_object_count += 1

        result = source_processor(
            manifest_by_source[source_key],
            source_rows,
            upstream=upstream,
            bucket_url=bucket_url,
            bri_pipeline_git_commit=bri_pipeline_git_commit,
            timeout_seconds=timeout_seconds,
            max_retries=max_retries,
        )

        if (
            result.input_eligible_chain_count
            != len(source_rows)
        ):
            raise BRIRunnerError(
                "Per-source BRI input count does not match "
                "task-selected eligible rows"
            )

        if not result.chain_accounting_valid:
            raise BRIRunnerError(
                "Per-source BRI chain accounting failed"
            )

        if result.source_downloaded:
            downloaded_source_object_count += 1

        if result.source_parsed:
            parsed_source_object_count += 1

        bri_records.extend(result.bri_records)
        processing_errors.extend(
            result.processing_errors
        )

    tables = _materialize_bri_tables(
        bri_records,
        processing_errors,
    )

    batch = BRIBatchResult(
        manifest_source_object_count=len(manifest),
        relevant_source_object_count=(
            relevant_source_object_count
        ),
        input_eligible_chain_count=len(eligible),
        downloaded_source_object_count=(
            downloaded_source_object_count
        ),
        parsed_source_object_count=(
            parsed_source_object_count
        ),
        tables=tables,
    )

    if not batch.chain_accounting_valid:
        raise BRIRunnerError(
            "Task-level Stage-3 BRI chain accounting failed"
        )

    return batch


def _write_bri_parquet_atomic(
    table: pa.Table,
    output_path: str | Path,
) -> Path:
    """Write one Stage-3 Parquet artifact atomically."""

    output = Path(output_path)
    output.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    temporary = output.with_suffix(
        output.suffix + ".tmp"
    )

    try:
        pq.write_table(
            table,
            temporary,
            compression="zstd",
            version="2.6",
        )
        temporary.replace(output)
    except Exception:
        if temporary.exists():
            temporary.unlink()

        raise

    return output


def write_bri_shards(
    tables: BRITables,
    output_root: str | Path,
    task_id: str | int,
) -> dict[str, Path]:
    """Write deterministic schema-enforced Stage-3 task shards."""

    task_id_text = _safe_bri_task_id(task_id)
    root = Path(output_root)

    if not tables.chains.schema.equals(
        STAGE3_BRI_CHAIN_SCHEMA,
        check_metadata=True,
    ):
        raise BRIRunnerError(
            "Cannot publish BRI chains with unexpected schema"
        )

    if not tables.processing_errors.schema.equals(
        STAGE3_BRI_PROCESSING_ERROR_SCHEMA,
        check_metadata=True,
    ):
        raise BRIRunnerError(
            "Cannot publish BRI errors with unexpected schema"
        )

    paths = {
        "chains": (
            root
            / "chains"
            / f"task_{task_id_text}.parquet"
        ),
        "processing_errors": (
            root
            / "processing_errors"
            / f"task_{task_id_text}.parquet"
        ),
    }

    _write_bri_parquet_atomic(
        tables.chains,
        paths["chains"],
    )
    _write_bri_parquet_atomic(
        tables.processing_errors,
        paths["processing_errors"],
    )

    return paths


def _sorted_counts(
    values: Iterable[str],
) -> dict[str, int]:
    """Return deterministic lexical-key counts."""

    counts = Counter(values)

    return {
        key: counts[key]
        for key in sorted(counts)
    }


def build_bri_task_summary(
    batch: BRIBatchResult,
    context: BRITaskContext,
) -> dict[str, Any]:
    """Build deterministic Stage-3 BRI task completion metadata."""

    if context.completed_at_utc is None:
        raise BRIRunnerError(
            "BRI task summary requires completed_at_utc"
        )

    if context.runtime_seconds is None:
        raise BRIRunnerError(
            "BRI task summary requires runtime_seconds"
        )

    if context.runtime_seconds < 0:
        raise BRIRunnerError(
            "BRI task runtime cannot be negative"
        )

    error_rows = (
        batch.tables.processing_errors.to_pylist()
    )

    chain_accounting_valid = (
        batch.input_eligible_chain_count
        == batch.tables.chains.num_rows
        + len(error_rows)
    )

    return {
        "summary_schema_name": (
            BRI_TASK_SUMMARY_SCHEMA_NAME
        ),
        "summary_schema_version": (
            BRI_TASK_SUMMARY_SCHEMA_VERSION
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
        "geometric_validation_finalizer_git_commit": (
            context.geometric_validation_finalizer_git_commit
        ),
        "bri_pipeline_git_commit": (
            context.bri_pipeline_git_commit
        ),
        "started_at_utc": context.started_at_utc,
        "completed_at_utc": context.completed_at_utc,
        "runtime_seconds": context.runtime_seconds,
        "slurm_job_id": context.slurm_job_id,
        "slurm_array_task_id": (
            context.slurm_array_task_id
        ),
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
        "input_eligible_chain_count": (
            batch.input_eligible_chain_count
        ),
        "bri_chain_count": (
            batch.tables.chains.num_rows
        ),
        "processing_error_count": len(error_rows),
        "processing_errors_by_stage": _sorted_counts(
            row["processing_stage"]
            for row in error_rows
        ),
        "processing_errors_by_type": _sorted_counts(
            row["error_type"]
            for row in error_rows
        ),
        "chain_accounting_valid": (
            chain_accounting_valid
        ),
    }


def write_bri_task_summary_atomic(
    summary: dict[str, Any],
    output_root: str | Path,
) -> Path:
    """Write the Stage-3 task completion marker atomically."""

    if (
        summary.get("summary_schema_name")
        != BRI_TASK_SUMMARY_SCHEMA_NAME
    ):
        raise BRIRunnerError(
            "Unexpected BRI task-summary schema name"
        )

    if (
        summary.get("summary_schema_version")
        != BRI_TASK_SUMMARY_SCHEMA_VERSION
    ):
        raise BRIRunnerError(
            "Unexpected BRI task-summary schema version"
        )

    if summary.get("chain_accounting_valid") is not True:
        raise BRIRunnerError(
            "Cannot publish BRI task summary: "
            "chain accounting failed"
        )

    task_id = _safe_bri_task_id(
        summary.get("task_id", "")
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


def publish_bri_batch(
    batch: BRIBatchResult,
    *,
    output_root: str | Path,
    context: BRITaskContext,
    started_perf_counter: float,
    shard_writer: Callable[
        ...,
        dict[str, Path],
    ] = write_bri_shards,
    summary_writer: Callable[
        ...,
        Path,
    ] = write_bri_task_summary_atomic,
    utc_now: Callable[[], str] = _utc_now_text,
    perf_counter: Callable[[], float] = time.perf_counter,
    peak_memory_reader: Callable[
        [],
        int | None,
    ] = _linux_process_peak_memory_bytes,
) -> BRITaskPublication:
    """Publish terminal shards first and task summary strictly last."""

    if not batch.chain_accounting_valid:
        raise BRIRunnerError(
            "Cannot publish Stage-3 BRI task: "
            "chain accounting failed"
        )

    if (
        not isinstance(
            started_perf_counter,
            (int, float),
        )
        or isinstance(started_perf_counter, bool)
    ):
        raise BRIRunnerError(
            "started_perf_counter must be numeric"
        )

    task_id = _safe_bri_task_id(
        context.task_id
    )

    # A stale completion marker must not survive an attempted
    # re-publication.
    existing_summary = (
        Path(output_root)
        / "summaries"
        / f"task_{task_id}.json"
    )

    if existing_summary.exists():
        existing_summary.unlink()

    # Terminal Parquet artifacts are published first.
    shard_paths = shard_writer(
        batch.tables,
        output_root,
        task_id,
    )

    completed_at_utc = utc_now()
    runtime_seconds = (
        perf_counter()
        - started_perf_counter
    )

    if runtime_seconds < 0:
        raise BRIRunnerError(
            "Monotonic Stage-3 BRI runtime cannot be negative"
        )

    completed_context = replace(
        context,
        completed_at_utc=completed_at_utc,
        runtime_seconds=runtime_seconds,
        peak_memory_bytes=peak_memory_reader(),
    )

    summary = build_bri_task_summary(
        batch,
        completed_context,
    )

    # This is deliberately the final publication operation and
    # therefore acts as the task completion marker.
    summary_path = summary_writer(
        summary,
        output_root,
    )

    return BRITaskPublication(
        shard_paths=shard_paths,
        summary_path=summary_path,
        summary=summary,
    )


def execute_bri_task(
    manifest_rows: Iterable[Mapping[str, Any]],
    eligible_rows: Iterable[Mapping[str, Any]],
    *,
    upstream: UpstreamGeometricValidation,
    output_root: str | Path,
    task_id: str | int,
    bucket_url: str,
    bri_pipeline_git_commit: str,
    timeout_seconds: int = 60,
    max_retries: int = 0,
    environ: Mapping[str, str] | None = None,
    batch_processor: Callable[
        ...,
        BRIBatchResult,
    ] = process_bri_batch,
    publisher: Callable[
        ...,
        BRITaskPublication,
    ] = publish_bri_batch,
    utc_now: Callable[[], str] = _utc_now_text,
    perf_counter: Callable[[], float] = time.perf_counter,
) -> BRITaskPublication:
    """Execute one complete logical Stage-3 BRI task."""

    task_id_text = _safe_bri_task_id(
        task_id
    )

    bri_pipeline_git_commit = _validated_git_commit(
        bri_pipeline_git_commit,
        field="bri_pipeline_git_commit",
    )

    slurm_job_id, slurm_array_task_id = (
        _slurm_environment(environ)
    )

    started_at_utc = utc_now()
    started_perf_counter = perf_counter()

    batch = batch_processor(
        manifest_rows,
        eligible_rows,
        upstream=upstream,
        bucket_url=bucket_url,
        bri_pipeline_git_commit=(
            bri_pipeline_git_commit
        ),
        timeout_seconds=timeout_seconds,
        max_retries=max_retries,
    )

    context = BRITaskContext(
        task_id=task_id_text,
        snapshot=upstream.snapshot,
        cleaning_protocol=(
            upstream.cleaning_protocol
        ),
        quality_pipeline_git_commit=(
            upstream.quality_pipeline_git_commit
        ),
        geometric_validation_pipeline_git_commit=(
            upstream.geometric_validation_pipeline_git_commit
        ),
        geometric_validation_finalizer_git_commit=(
            upstream.geometric_validation_finalizer_git_commit
        ),
        bri_pipeline_git_commit=(
            bri_pipeline_git_commit
        ),
        started_at_utc=started_at_utc,
        slurm_job_id=slurm_job_id,
        slurm_array_task_id=(
            slurm_array_task_id
        ),
    )

    return publisher(
        batch,
        output_root=output_root,
        context=context,
        started_perf_counter=(
            started_perf_counter
        ),
        utc_now=utc_now,
        perf_counter=perf_counter,
    )
