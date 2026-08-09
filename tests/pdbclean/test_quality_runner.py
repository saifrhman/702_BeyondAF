"""Tests for quality-stage batch orchestration."""

import pytest

from pdbclean.mmcif_parser import ChainObservation
from pdbclean.quality_runner import (
    QualityRunnerError,
    select_configured_model_chains,
)


def _chain(
    model_id: int,
    label_chain_id: str,
) -> ChainObservation:
    return ChainObservation(
        pdb_id="test",
        model_id=model_id,
        label_chain_id=label_chain_id,
        entry_has_polypeptide=True,
    )


def test_select_configured_model_chains_keeps_only_model_one() -> None:
    chains = [
        _chain(1, "A"),
        _chain(2, "A"),
        _chain(1, "B"),
        _chain(3, "C"),
    ]

    selected = select_configured_model_chains(
        chains,
        {
            "models": {
                "policy": "first_model",
                "model_id": 1,
            }
        },
    )

    assert [
        (chain.model_id, chain.label_chain_id)
        for chain in selected
    ] == [
        (1, "A"),
        (1, "B"),
    ]


def test_select_configured_model_chains_respects_configured_model_id() -> None:
    chains = [
        _chain(1, "A"),
        _chain(2, "A"),
    ]

    selected = select_configured_model_chains(
        chains,
        {
            "models": {
                "policy": "first_model",
                "model_id": 2,
            }
        },
    )

    assert len(selected) == 1
    assert selected[0].model_id == 2


def test_select_configured_model_chains_rejects_unsupported_policy() -> None:
    with pytest.raises(
        QualityRunnerError,
        match="Unsupported model-selection policy",
    ):
        select_configured_model_chains(
            [_chain(1, "A")],
            {
                "models": {
                    "policy": "all_models",
                    "model_id": 1,
                }
            },
        )
