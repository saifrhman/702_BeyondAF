"""Tests for mmCIF parser observation types."""

from pdbclean.mmcif_parser import (
    AtomObservation,
    ChainObservation,
)


def test_atom_observation_fields() -> None:
    atom = AtomObservation(
        model_id=1,
        label_chain_id="A",
        auth_chain_id="A",
        entity_id="1",
        label_seq_id=10,
        auth_seq_id="10",
        residue_name="GLY",
        atom_name="CA",
        alt_id=None,
        occupancy=1.0,
        x=1.0,
        y=2.0,
        z=3.0,
    )

    assert atom.model_id == 1
    assert atom.label_chain_id == "A"
    assert atom.label_seq_id == 10
    assert atom.atom_name == "CA"
    assert atom.occupancy == 1.0


def test_chain_canonical_key() -> None:
    chain = ChainObservation(
        pdb_id="100d",
        model_id=1,
        label_chain_id="A",
    )

    assert chain.canonical_key == ("100d", 1, "A")


def test_chain_atom_count() -> None:
    chain = ChainObservation(
        pdb_id="100d",
        model_id=1,
        label_chain_id="A",
    )

    chain.atoms.append(
        AtomObservation(
            model_id=1,
            label_chain_id="A",
            auth_chain_id="A",
            entity_id="1",
            label_seq_id=1,
            auth_seq_id="1",
            residue_name="ALA",
            atom_name="N",
            alt_id=None,
            occupancy=1.0,
            x=0.0,
            y=0.0,
            z=0.0,
        )
    )

    assert chain.atom_count == 1
