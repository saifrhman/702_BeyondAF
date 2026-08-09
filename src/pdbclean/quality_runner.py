"""Batch orchestration for Protocol 3.2 quality cleaning."""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from typing import Any

from pdbclean.cleaning import (
    clean_protocol32_chain,
    is_protocol32_candidate,
)
from pdbclean.gold import (
    GoldChainRecords,
    GoldProvenance,
    materialize_gold_chain,
)
from pdbclean.mmcif_parser import (
    ChainObservation,
    MMCIFParseError,
    parse_coordinate_mmcif_bytes,
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
