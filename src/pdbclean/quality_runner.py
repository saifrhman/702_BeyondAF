"""Batch orchestration for Protocol 3.2 quality cleaning."""

from __future__ import annotations

from collections.abc import Callable, Iterable, Mapping
from dataclasses import dataclass
from typing import Any

from pdbclean.cleaning import (
    clean_protocol32_chain,
    is_protocol32_candidate,
)
from pdbclean.gold import (
    GoldChainRecords,
    GoldProvenance,
    GoldTables,
    gold_records_to_tables,
    materialize_gold_chain,
)
from pdbclean.mmcif_parser import (
    ChainObservation,
    MMCIFParseError,
    parse_coordinate_mmcif_bytes,
)
from pdbclean.snapshot import (
    SnapshotError,
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


class QualityRunnerError(RuntimeError):
    """Raised when quality-stage orchestration cannot proceed safely."""


def select_configured_model_chains(
    chains: Iterable[ChainObservation],
    selection_config: dict[str, Any],
) -> list[ChainObservation]:
    """Select chains belonging to the configured structural model."""

    models = selection_config.get("models")

    if not isinstance(models, dict):
        raise QualityRunnerError(
            "selection.models must be a mapping"
        )

    policy = models.get("policy")

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
    downloader: Callable[..., bytes] = download_verified_s3_object_bytes,
) -> SourceQualityResult:
    """Download, verify, parse, and quality-clean one manifest source row."""

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

    try:
        compressed_bytes = downloader(
            bucket_url=bucket_url,
            s3_key=s3_key,
            expected_size_bytes=size_bytes,
            expected_etag=etag,
            timeout_seconds=timeout_seconds,
        )
    except SnapshotError as exc:
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
                    "error_type": type(exc).__name__,
                    "error_message": str(exc),
                    "source_mmcif_key": s3_key,
                    "source_etag": etag,
                    "pipeline_git_commit": pipeline_git_commit,
                },
            ),
            source_failed=True,
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
    source_processor: Callable[..., SourceQualityResult] = process_manifest_source,
) -> QualityBatchResult:
    """Process and aggregate one manifest partition in memory."""

    records: list[GoldChainRecords] = []
    processing_errors: list[dict[str, Any]] = []

    input_source_object_count = 0
    successful_source_object_count = 0
    failed_source_object_count = 0

    parsed_silver_chain_count = 0
    selected_silver_chain_count = 0
    candidate_entry_count = 0
    candidate_chain_count = 0

    for manifest_row in manifest_rows:
        input_source_object_count += 1

        result = source_processor(
            manifest_row,
            bucket_url=bucket_url,
            selection_config=selection_config,
            cleaning_protocol=cleaning_protocol,
            pipeline_git_commit=pipeline_git_commit,
            timeout_seconds=timeout_seconds,
        )

        if result.source_failed:
            failed_source_object_count += 1
        else:
            successful_source_object_count += 1

        parsed_silver_chain_count += result.parsed_silver_chain_count
        selected_silver_chain_count += result.selected_silver_chain_count
        candidate_entry_count += result.candidate_entry_count
        candidate_chain_count += result.candidate_chain_count

        records.extend(result.gold_records)
        processing_errors.extend(result.processing_errors)

    tables = gold_records_to_tables(
        records,
        processing_errors=processing_errors,
    )

    return QualityBatchResult(
        input_source_object_count=input_source_object_count,
        successful_source_object_count=successful_source_object_count,
        failed_source_object_count=failed_source_object_count,
        parsed_silver_chain_count=parsed_silver_chain_count,
        selected_silver_chain_count=selected_silver_chain_count,
        candidate_entry_count=candidate_entry_count,
        candidate_chain_count=candidate_chain_count,
        tables=tables,
    )
