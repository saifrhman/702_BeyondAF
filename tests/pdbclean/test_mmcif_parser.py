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
    assert atom.occupancy_raw is None


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

import gzip

import pytest

from pdbclean.mmcif_parser import (
    MMCIFParseError,
    parse_coordinate_mmcif_bytes,
)


def _gzip_cif(text: str) -> bytes:
    return gzip.compress(text.encode("utf-8"))


def test_parse_coordinate_mmcif_into_chain() -> None:
    cif = """data_test
#
loop_
_entity_poly.entity_id
_entity_poly.type
1 'polypeptide(L)'
#
loop_
_atom_site.pdbx_PDB_model_num
_atom_site.label_asym_id
_atom_site.auth_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.auth_seq_id
_atom_site.label_comp_id
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.occupancy
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
1 A X 1 1 10 ALA N  . .    0.0 0.0 0.0
1 A X 1 1 10 ALA CA . ?    1.0 0.0 0.0
1 A X 1 1 10 ALA C  . 1.00 2.0 0.0 0.0
1 A X 1 2 11 GLY N  . 1.00 3.0 0.0 0.0
1 A X 1 2 11 GLY CA . 0.50 4.0 0.0 0.0
1 A X 1 2 11 GLY C  A 1.00 5.0 0.0 0.0
#
"""

    chains = parse_coordinate_mmcif_bytes(
        _gzip_cif(cif),
        pdb_id="TEST",
    )

    assert len(chains) == 1

    chain = chains[0]

    assert chain.pdb_id == "test"
    assert chain.model_id == 1
    assert chain.label_chain_id == "A"
    assert chain.auth_chain_id == "X"
    assert chain.entity_id == "1"
    assert chain.polymer_type == "polypeptide(L)"
    assert chain.atom_count == 6

    assert chain.atoms[0].label_seq_id == 1
    assert chain.atoms[0].residue_name == "ALA"
    assert chain.atoms[0].atom_name == "N"

    # Numeric normalization alone cannot reproduce BRI disorder_check:
    # both become None numerically, while the raw tokens remain distinct.
    assert chain.atoms[0].occupancy is None
    assert chain.atoms[0].occupancy_raw == "."
    assert chain.atoms[1].occupancy is None
    assert chain.atoms[1].occupancy_raw == "?"

    assert chain.atoms[4].occupancy == 0.5
    assert chain.atoms[4].occupancy_raw == "0.50"
    assert chain.atoms[5].alt_id == "A"


def test_invalid_gzip_is_rejected() -> None:
    with pytest.raises(
        MMCIFParseError,
        match="invalid or truncated gzip stream",
    ):
        parse_coordinate_mmcif_bytes(
            b"not-a-gzip-file",
            pdb_id="bad1",
        )
