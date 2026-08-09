"""Tests for quality-stage batch orchestration."""

import pytest

from pdbclean.mmcif_parser import AtomObservation, ChainObservation
from pdbclean.quality_runner import (
    QualityRunnerError,
    candidate_accounting,
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



def _candidate_chain(
    pdb_id: str,
    label_chain_id: str,
) -> ChainObservation:
    atom = AtomObservation(
        model_id=1,
        label_chain_id=label_chain_id,
        auth_chain_id=label_chain_id,
        entity_id="1",
        label_seq_id=1,
        auth_seq_id="1",
        residue_name="ALA",
        atom_name="CA",
        alt_id=None,
        occupancy=1.0,
        x=0.0,
        y=0.0,
        z=0.0,
        group_pdb="ATOM",
        occupancy_raw="1.00",
    )

    return ChainObservation(
        pdb_id=pdb_id,
        model_id=1,
        label_chain_id=label_chain_id,
        entry_has_polypeptide=True,
        atoms=[atom],
    )


def _non_candidate_chain(
    pdb_id: str,
    label_chain_id: str,
) -> ChainObservation:
    atom = AtomObservation(
        model_id=1,
        label_chain_id=label_chain_id,
        auth_chain_id=label_chain_id,
        entity_id="2",
        label_seq_id=None,
        auth_seq_id=None,
        residue_name="SO4",
        atom_name="S",
        alt_id=None,
        occupancy=1.0,
        x=0.0,
        y=0.0,
        z=0.0,
        group_pdb="HETATM",
        occupancy_raw="1.00",
    )

    return ChainObservation(
        pdb_id=pdb_id,
        model_id=1,
        label_chain_id=label_chain_id,
        entry_has_polypeptide=True,
        atoms=[atom],
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


def test_candidate_accounting_counts_unique_entries_and_chains() -> None:
    chains = [
        _candidate_chain("1abc", "A"),
        _candidate_chain("1abc", "B"),
        _candidate_chain("2def", "A"),
    ]

    entry_count, chain_count = candidate_accounting(chains)

    assert entry_count == 2
    assert chain_count == 3


def test_candidate_accounting_excludes_non_candidates() -> None:
    chains = [
        _candidate_chain("1abc", "A"),
        _non_candidate_chain("1abc", "B"),
        _non_candidate_chain("2def", "A"),
    ]

    entry_count, chain_count = candidate_accounting(chains)

    assert entry_count == 1
    assert chain_count == 1


def test_candidate_accounting_empty_input() -> None:
    assert candidate_accounting([]) == (0, 0)
