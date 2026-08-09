"""Batch orchestration for Protocol 3.2 quality cleaning."""

from __future__ import annotations

from collections.abc import Iterable
from typing import Any

from pdbclean.mmcif_parser import ChainObservation


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
