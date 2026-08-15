"""Production orchestration for Stage-3 Definition 3.4 BRI.

This module consumes only the canonical population published by the
completed post-cleaning geometric-validation stage.

The scientific BRI calculation itself lives in ``pdbclean.bri``.
"""

from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
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
from pdbclean.schemas import STAGE3_ELIGIBLE_CHAIN_SCHEMA
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
