"""Batch orchestration for Protocol 3.2 quality cleaning."""

from __future__ import annotations

from collections.abc import Callable, Iterable, Mapping
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from datetime import datetime, timezone
import os
from pathlib import Path
import resource
import time
from typing import Any

from pdbclean.cleaning import (
    clean_protocol32_chain,
    is_protocol32_candidate,
)
from pdbclean.gold import (
    GoldChainRecords,
    GoldProvenance,
    GoldTables,
    QualityTaskContext,
    build_quality_task_summary,
    gold_records_to_tables,
    materialize_gold_chain,
    write_gold_quality_shards,
    write_quality_task_summary_atomic,
)
from pdbclean.mmcif_parser import (
    ChainObservation,
    MMCIFParseError,
    parse_coordinate_mmcif_bytes,
)
from pdbclean.snapshot import (
    SnapshotError,
    SnapshotTransportError,
    download_verified_s3_object_bytes,
)


@dataclass(frozen=True)
class SourceQualityResult:
    """Quality-processing outcome for one verified source mmCIF."""

    pdb_id: str
    parsed_silver_chain_count: int
    selected_silver_chain_count: int
    candidate_entry_count: int
    candidate_chain_count: int

    gold_records: tuple[GoldChainRecords, ...] = ()
    processing_errors: tuple[dict[str, Any], ...] = ()

    # True only when the source could not reach chain-level processing,
    # for example because the mmCIF could not be parsed.
    source_failed: bool = False


@dataclass(frozen=True)
class QualityBatchResult:
    """Aggregated in-memory outcome for one manifest partition."""

    input_source_object_count: int
    successful_source_object_count: int
    failed_source_object_count: int

    parsed_silver_chain_count: int
    selected_silver_chain_count: int
    candidate_entry_count: int
    candidate_chain_count: int

    tables: GoldTables


@dataclass(frozen=True)
class QualityTaskPublication:
    """Published outputs and summary for one completed quality task."""

    shard_paths: dict[str, Path]
    summary_path: Path
    summary: dict[str, Any]


@dataclass(frozen=True)
class QualityExecutionMetadata:
    """Runtime metadata captured for one task execution."""

    started_at_utc: str
    completed_at_utc: str
    runtime_seconds: float
    slurm_job_id: str | None
    slurm_array_task_id: str | None
    peak_memory_bytes: int | None


def _utc_now_text() -> str:
    """Return a timezone-explicit UTC timestamp for task provenance."""

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


def quality_stage_output_root(
    storage_output_root: str | Path,
    *,
    snapshot: str,
    protocol_version: str,
) -> Path:
    """Return the isolated output root for one quality stage."""

    if (
        not snapshot
        or "/" in snapshot
        or "\\" in snapshot
        or snapshot in {".", ".."}
    ):
        raise QualityRunnerError(
            f"Unsafe quality-stage snapshot: {snapshot!r}"
        )

    if (
        not protocol_version
        or "/" in protocol_version
        or "\\" in protocol_version
        or protocol_version in {".", ".."}
    ):
        raise QualityRunnerError(
            "Unsafe quality-stage protocol version: "
            f"{protocol_version!r}"
        )

    return (
        Path(storage_output_root)
        / snapshot
        / protocol_version
        / "quality"
    )


class QualityRunnerError(RuntimeError):
    """Raised when quality-stage orchestration cannot proceed safely."""


def select_configured_model_chains(
    chains: Iterable[ChainObservation],
    selection_config: dict[str, Any],
) -> list[ChainObservation]:
    """Select chains according to the configured structural-model policy."""

    models = selection_config.get("models")

    if not isinstance(models, dict):
        raise QualityRunnerError(
            "selection.models must be a mapping"
        )

    policy = models.get("policy")

    if policy == "all_models":
        return list(chains)

    if policy != "first_model":
        raise QualityRunnerError(
            f"Unsupported model-selection policy: {policy!r}"
        )

    model_id = models.get("model_id")

    if not isinstance(model_id, int) or isinstance(model_id, bool):
        raise QualityRunnerError(
            "selection.models.model_id must be an integer"
        )

    if model_id <= 0:
        raise QualityRunnerError(
            "selection.models.model_id must be positive"
        )

    return [
        chain
        for chain in chains
        if chain.model_id == model_id
    ]



def candidate_accounting(
    chains: Iterable[ChainObservation],
) -> tuple[int, int]:
    """Return candidate-entry and candidate-chain counts."""

    candidate_chains = [
        chain
        for chain in chains
        if is_protocol32_candidate(chain)
    ]

    candidate_entry_count = len(
        {
            chain.pdb_id
            for chain in candidate_chains
        }
    )

    return (
        candidate_entry_count,
        len(candidate_chains),
    )



def process_verified_mmcif_bytes(
    compressed_bytes: bytes,
    *,
    pdb_id: str,
    selection_config: dict[str, Any],
    provenance: GoldProvenance,
) -> SourceQualityResult:
    """Parse and quality-clean one already verified source object."""

    normalized_pdb_id = pdb_id.lower()

    try:
        parsed_chains = parse_coordinate_mmcif_bytes(
            compressed_bytes,
            pdb_id=normalized_pdb_id,
        )
    except MMCIFParseError as exc:
        return SourceQualityResult(
            pdb_id=normalized_pdb_id,
            parsed_silver_chain_count=0,
            selected_silver_chain_count=0,
            candidate_entry_count=0,
            candidate_chain_count=0,
            processing_errors=(
                {
                    "snapshot": provenance.snapshot,
                    "pdb_id": normalized_pdb_id,
                    "model_id": None,
                    "label_chain_id": None,
                    "processing_stage": "mmcif_parse",
                    "error_type": type(exc).__name__,
                    "error_message": str(exc),
                    "source_mmcif_key": provenance.source_mmcif_key,
                    "source_etag": provenance.source_etag,
                    "pipeline_git_commit": provenance.pipeline_git_commit,
                },
            ),
            source_failed=True,
        )

    selected_chains = select_configured_model_chains(
        parsed_chains,
        selection_config,
    )

    candidate_entry_count, candidate_chain_count = candidate_accounting(
        selected_chains
    )

    gold_records: list[GoldChainRecords] = []
    processing_errors: list[dict[str, Any]] = []

    for chain in selected_chains:
        try:
            cleaning_result = clean_protocol32_chain(chain)
            gold_records.append(
                materialize_gold_chain(
                    chain,
                    cleaning_result,
                    provenance,
                )
            )
        except Exception as exc:
            processing_errors.append(
                {
                    "snapshot": provenance.snapshot,
                    "pdb_id": chain.pdb_id,
                    "model_id": chain.model_id,
                    "label_chain_id": chain.label_chain_id,
                    "processing_stage": "quality_cleaning",
                    "error_type": type(exc).__name__,
                    "error_message": str(exc),
                    "source_mmcif_key": provenance.source_mmcif_key,
                    "source_etag": provenance.source_etag,
                    "pipeline_git_commit": provenance.pipeline_git_commit,
                }
            )

    return SourceQualityResult(
        pdb_id=normalized_pdb_id,
        parsed_silver_chain_count=len(parsed_chains),
        selected_silver_chain_count=len(selected_chains),
        candidate_entry_count=candidate_entry_count,
        candidate_chain_count=candidate_chain_count,
        gold_records=tuple(gold_records),
        processing_errors=tuple(processing_errors),
        source_failed=False,
    )



def process_manifest_source(
    manifest_row: Mapping[str, Any],
    *,
    bucket_url: str,
    selection_config: dict[str, Any],
    cleaning_protocol: str,
    pipeline_git_commit: str,
    timeout_seconds: int = 60,
    max_retries: int = 0,
    downloader: Callable[..., bytes] = download_verified_s3_object_bytes,
) -> SourceQualityResult:
    """Download, verify, parse, and quality-clean one manifest source row."""

    if (
        not isinstance(max_retries, int)
        or isinstance(max_retries, bool)
        or max_retries < 0
    ):
        raise QualityRunnerError(
            "max_retries must be a non-negative integer"
        )

    required_fields = (
        "snapshot",
        "pdb_id",
        "s3_key",
        "size_bytes",
        "etag",
    )

    missing = [
        field
        for field in required_fields
        if field not in manifest_row
    ]

    if missing:
        raise QualityRunnerError(
            "Manifest row missing required field(s): "
            + ", ".join(missing)
        )

    snapshot = manifest_row["snapshot"]
    pdb_id = manifest_row["pdb_id"]
    s3_key = manifest_row["s3_key"]
    size_bytes = manifest_row["size_bytes"]
    etag = manifest_row["etag"]

    if not isinstance(snapshot, str) or not snapshot:
        raise QualityRunnerError(
            "Manifest row snapshot must be a non-empty string"
        )

    if not isinstance(pdb_id, str) or not pdb_id:
        raise QualityRunnerError(
            "Manifest row pdb_id must be a non-empty string"
        )

    if not isinstance(s3_key, str) or not s3_key:
        raise QualityRunnerError(
            "Manifest row s3_key must be a non-empty string"
        )

    if (
        not isinstance(size_bytes, int)
        or isinstance(size_bytes, bool)
        or size_bytes <= 0
    ):
        raise QualityRunnerError(
            "Manifest row size_bytes must be a positive integer"
        )

    if not isinstance(etag, str) or not etag.strip().strip('"'):
        raise QualityRunnerError(
            "Manifest row etag must be a non-empty string"
        )

    normalized_pdb_id = pdb_id.lower()

    provenance = GoldProvenance(
        snapshot=snapshot,
        source_mmcif_key=s3_key,
        source_etag=etag,
        cleaning_protocol=cleaning_protocol,
        pipeline_git_commit=pipeline_git_commit,
    )

    compressed_bytes: bytes | None = None
    source_error: SnapshotError | None = None

    for attempt in range(max_retries + 1):
        try:
            compressed_bytes = downloader(
                bucket_url=bucket_url,
                s3_key=s3_key,
                expected_size_bytes=size_bytes,
                expected_etag=etag,
                timeout_seconds=timeout_seconds,
            )
            source_error = None
            break
        except SnapshotTransportError as exc:
            source_error = exc

            if attempt == max_retries:
                break
        except SnapshotError as exc:
            # Deterministic verification failures and other non-transport
            # snapshot failures must not be retried.
            source_error = exc
            break

    if source_error is not None:
        return SourceQualityResult(
            pdb_id=normalized_pdb_id,
            parsed_silver_chain_count=0,
            selected_silver_chain_count=0,
            candidate_entry_count=0,
            candidate_chain_count=0,
            processing_errors=(
                {
                    "snapshot": snapshot,
                    "pdb_id": normalized_pdb_id,
                    "model_id": None,
                    "label_chain_id": None,
                    "processing_stage": "source_download_verify",
                    "error_type": type(source_error).__name__,
                    "error_message": str(source_error),
                    "source_mmcif_key": s3_key,
                    "source_etag": etag,
                    "pipeline_git_commit": pipeline_git_commit,
                },
            ),
            source_failed=True,
        )

    if compressed_bytes is None:
        raise QualityRunnerError(
            "Downloader completed without bytes or a source error"
        )

    return process_verified_mmcif_bytes(
        compressed_bytes,
        pdb_id=normalized_pdb_id,
        selection_config=selection_config,
        provenance=provenance,
    )



def process_manifest_batch(
    manifest_rows: Iterable[Mapping[str, Any]],
    *,
    bucket_url: str,
    selection_config: dict[str, Any],
    cleaning_protocol: str,
    pipeline_git_commit: str,
    timeout_seconds: int = 60,
    max_retries: int = 0,
    download_concurrency: int = 1,
    source_processor: Callable[
        ..., SourceQualityResult
    ] = process_manifest_source,
) -> QualityBatchResult:
    """Process and aggregate one manifest partition deterministically."""

    if (
        not isinstance(max_retries, int)
        or isinstance(max_retries, bool)
        or max_retries < 0
    ):
        raise QualityRunnerError(
            "max_retries must be a non-negative integer"
        )

    if (
        not isinstance(download_concurrency, int)
        or isinstance(download_concurrency, bool)
        or download_concurrency <= 0
    ):
        raise QualityRunnerError(
            "download_concurrency must be a positive integer"
        )

    # Materialise the partition exactly once. Its existing validated
    # manifest order defines deterministic aggregation order.
    rows = list(manifest_rows)

    records: list[GoldChainRecords] = []
    processing_errors: list[dict[str, Any]] = []

    input_source_object_count = len(rows)
    successful_source_object_count = 0
    failed_source_object_count = 0

    parsed_silver_chain_count = 0
    selected_silver_chain_count = 0
    candidate_entry_count = 0
    candidate_chain_count = 0

    def process_row(
        manifest_row: Mapping[str, Any],
    ) -> SourceQualityResult:
        return source_processor(
            manifest_row,
            bucket_url=bucket_url,
            selection_config=selection_config,
            cleaning_protocol=cleaning_protocol,
            pipeline_git_commit=pipeline_git_commit,
            timeout_seconds=timeout_seconds,
            max_retries=max_retries,
        )

    # executor.map yields results in input order even when workers finish
    # out of order. Therefore concurrency cannot reorder Gold records or
    # processing-error lineage relative to the validated manifest.
    with ThreadPoolExecutor(
        max_workers=download_concurrency
    ) as executor:
        results = executor.map(process_row, rows)

        for result in results:
            if result.source_failed:
                failed_source_object_count += 1
            else:
                successful_source_object_count += 1

            parsed_silver_chain_count += (
                result.parsed_silver_chain_count
            )
            selected_silver_chain_count += (
                result.selected_silver_chain_count
            )
            candidate_entry_count += (
                result.candidate_entry_count
            )
            candidate_chain_count += (
                result.candidate_chain_count
            )

            records.extend(result.gold_records)
            processing_errors.extend(
                result.processing_errors
            )

    tables = gold_records_to_tables(
        records,
        processing_errors=processing_errors,
    )

    return QualityBatchResult(
        input_source_object_count=input_source_object_count,
        successful_source_object_count=(
            successful_source_object_count
        ),
        failed_source_object_count=failed_source_object_count,
        parsed_silver_chain_count=parsed_silver_chain_count,
        selected_silver_chain_count=selected_silver_chain_count,
        candidate_entry_count=candidate_entry_count,
        candidate_chain_count=candidate_chain_count,
        tables=tables,
    )

def publish_quality_batch(
    batch: QualityBatchResult,
    *,
    output_root: str | Path,
    task_id: str | int,
    snapshot: str,
    cleaning_protocol: str,
    pipeline_git_commit: str,
    started_at_utc: str,
    started_perf_counter: float,
    slurm_job_id: str | None = None,
    slurm_array_task_id: str | None = None,
    shard_writer: Callable[..., dict[str, Path]] = write_gold_quality_shards,
    summary_writer: Callable[..., Path] = write_quality_task_summary_atomic,
    utc_now: Callable[[], str] = _utc_now_text,
    perf_counter: Callable[[], float] = time.perf_counter,
    peak_memory_reader: Callable[
        [], int | None
    ] = _linux_process_peak_memory_bytes,
) -> QualityTaskPublication:
    """Validate and publish one quality task, writing its summary last."""

    task_id_text = str(task_id)

    # Validate path safety before any Parquet shard can be created.
    if (
        not task_id_text
        or "/" in task_id_text
        or "\\" in task_id_text
        or task_id_text in {".", ".."}
    ):
        raise QualityRunnerError(
            f"Unsafe quality-task task_id: {task_id_text!r}"
        )

    if (
        not isinstance(started_perf_counter, (int, float))
        or isinstance(started_perf_counter, bool)
    ):
        raise QualityRunnerError(
            "started_perf_counter must be numeric"
        )

    # Build a provisional context solely to validate accounting before
    # publishing any task output. Timing fields do not affect accounting.
    provisional_context = QualityTaskContext(
        task_id=task_id_text,
        snapshot=snapshot,
        cleaning_protocol=cleaning_protocol,
        pipeline_git_commit=pipeline_git_commit,
        started_at_utc=started_at_utc,
        completed_at_utc=started_at_utc,
        runtime_seconds=0.0,
        input_source_object_count=batch.input_source_object_count,
        successful_source_object_count=(
            batch.successful_source_object_count
        ),
        failed_source_object_count=batch.failed_source_object_count,
        parsed_silver_chain_count=batch.parsed_silver_chain_count,
        selected_silver_chain_count=batch.selected_silver_chain_count,
        candidate_entry_count=batch.candidate_entry_count,
        candidate_chain_count=batch.candidate_chain_count,
        slurm_job_id=slurm_job_id,
        slurm_array_task_id=slurm_array_task_id,
        peak_memory_bytes=None,
    )

    provisional_summary = build_quality_task_summary(
        batch.tables,
        provisional_context,
    )

    # Accounting must be valid before any task output is published.
    if (
        provisional_summary.get("source_object_accounting_valid")
        is not True
    ):
        raise QualityRunnerError(
            "Cannot publish quality task: "
            "source-object accounting failed"
        )

    if (
        provisional_summary.get("selected_chain_accounting_valid")
        is not True
    ):
        raise QualityRunnerError(
            "Cannot publish quality task: "
            "selected-chain accounting failed"
        )

    # All five deterministic Parquet shards are published before the
    # task-level completion marker.
    shard_paths = shard_writer(
        batch.tables,
        output_root,
        task_id_text,
    )

    # Completion metadata deliberately includes Parquet publication.
    completed_at_utc = utc_now()
    runtime_seconds = perf_counter() - started_perf_counter

    if runtime_seconds < 0:
        raise QualityRunnerError(
            "Monotonic task runtime cannot be negative"
        )

    execution_metadata = QualityExecutionMetadata(
        started_at_utc=started_at_utc,
        completed_at_utc=completed_at_utc,
        runtime_seconds=runtime_seconds,
        slurm_job_id=slurm_job_id,
        slurm_array_task_id=slurm_array_task_id,
        peak_memory_bytes=peak_memory_reader(),
    )

    context = QualityTaskContext(
        task_id=task_id_text,
        snapshot=snapshot,
        cleaning_protocol=cleaning_protocol,
        pipeline_git_commit=pipeline_git_commit,
        started_at_utc=execution_metadata.started_at_utc,
        completed_at_utc=execution_metadata.completed_at_utc,
        runtime_seconds=execution_metadata.runtime_seconds,
        input_source_object_count=batch.input_source_object_count,
        successful_source_object_count=(
            batch.successful_source_object_count
        ),
        failed_source_object_count=batch.failed_source_object_count,
        parsed_silver_chain_count=batch.parsed_silver_chain_count,
        selected_silver_chain_count=batch.selected_silver_chain_count,
        candidate_entry_count=batch.candidate_entry_count,
        candidate_chain_count=batch.candidate_chain_count,
        slurm_job_id=execution_metadata.slurm_job_id,
        slurm_array_task_id=execution_metadata.slurm_array_task_id,
        peak_memory_bytes=execution_metadata.peak_memory_bytes,
    )

    summary = build_quality_task_summary(
        batch.tables,
        context,
    )

    # The summary is intentionally published last and acts as the
    # task-level completion marker.
    summary_path = summary_writer(
        summary,
        output_root,
    )

    return QualityTaskPublication(
        shard_paths=shard_paths,
        summary_path=summary_path,
        summary=summary,
    )


def execute_quality_task(
    manifest_rows: Iterable[Mapping[str, Any]],
    *,
    output_root: str | Path,
    task_id: str | int,
    snapshot: str,
    bucket_url: str,
    selection_config: dict[str, Any],
    cleaning_protocol: str,
    pipeline_git_commit: str,
    timeout_seconds: int = 60,
    max_retries: int = 0,
    download_concurrency: int = 1,
    environ: Mapping[str, str] | None = None,
    batch_processor: Callable[..., QualityBatchResult] = process_manifest_batch,
    publisher: Callable[..., QualityTaskPublication] = publish_quality_batch,
    utc_now: Callable[[], str] = _utc_now_text,
    perf_counter: Callable[[], float] = time.perf_counter,
) -> QualityTaskPublication:
    """Execute one complete quality-cleaning task."""

    slurm_job_id, slurm_array_task_id = _slurm_environment(environ)

    # Start metadata is captured immediately before source processing.
    started_at_utc = utc_now()
    started_perf_counter = perf_counter()

    batch = batch_processor(
        manifest_rows,
        bucket_url=bucket_url,
        selection_config=selection_config,
        cleaning_protocol=cleaning_protocol,
        pipeline_git_commit=pipeline_git_commit,
        timeout_seconds=timeout_seconds,
        max_retries=max_retries,
        download_concurrency=download_concurrency,
    )

    # publish_quality_batch writes Parquet shards first, captures final
    # execution metadata, and publishes the summary completion marker last.
    return publisher(
        batch,
        output_root=output_root,
        task_id=task_id,
        snapshot=snapshot,
        cleaning_protocol=cleaning_protocol,
        pipeline_git_commit=pipeline_git_commit,
        started_at_utc=started_at_utc,
        started_perf_counter=started_perf_counter,
        slurm_job_id=slurm_job_id,
        slurm_array_task_id=slurm_array_task_id,
        utc_now=utc_now,
        perf_counter=perf_counter,
    )
