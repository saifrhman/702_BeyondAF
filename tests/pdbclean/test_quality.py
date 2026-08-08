"""Tests for PDBClean scientific quality rules."""

from pdbclean.mmcif_parser import ChainObservation
from pdbclean.quality import evaluate_q001_protein_polymer


ALLOWED_PROTEIN_TYPES = {
    "polypeptide(L)",
    "polypeptide(D)",
}


def test_q001_accepts_allowed_protein_polymer() -> None:
    chain = ChainObservation(
        pdb_id="1aam",
        model_id=1,
        label_chain_id="A",
        polymer_type="polypeptide(L)",
    )

    result = evaluate_q001_protein_polymer(
        chain,
        allowed_polymer_types=ALLOWED_PROTEIN_TYPES,
    )

    assert result.rule_id == "Q001"
    assert result.passed is True
    assert result.reason == "allowed_protein_polymer_type"


def test_q001_rejects_missing_polymer_type() -> None:
    chain = ChainObservation(
        pdb_id="1aam",
        model_id=1,
        label_chain_id="B",
        polymer_type=None,
    )

    result = evaluate_q001_protein_polymer(
        chain,
        allowed_polymer_types=ALLOWED_PROTEIN_TYPES,
    )

    assert result.passed is False
    assert result.reason == "missing_or_nonpolymer_entity_poly_type"


def test_q001_rejects_disallowed_polymer_type() -> None:
    chain = ChainObservation(
        pdb_id="test",
        model_id=1,
        label_chain_id="A",
        polymer_type="polydeoxyribonucleotide",
    )

    result = evaluate_q001_protein_polymer(
        chain,
        allowed_polymer_types=ALLOWED_PROTEIN_TYPES,
    )

    assert result.passed is False
    assert (
        result.reason
        == "disallowed_polymer_type:polydeoxyribonucleotide"
    )

from pdbclean.mmcif_parser import AtomObservation
from pdbclean.quality import evaluate_q002_occupancy


def _atom(
    *,
    occupancy: float | None = 1.0,
    alt_id: str | None = None,
) -> AtomObservation:
    return AtomObservation(
        model_id=1,
        label_chain_id="A",
        auth_chain_id="A",
        entity_id="1",
        label_seq_id=1,
        auth_seq_id="1",
        residue_name="ALA",
        atom_name="CA",
        alt_id=alt_id,
        occupancy=occupancy,
        x=0.0,
        y=0.0,
        z=0.0,
    )


def test_q002_accepts_full_occupancy_without_altloc() -> None:
    chain = ChainObservation(
        pdb_id="test",
        model_id=1,
        label_chain_id="A",
        polymer_type="polypeptide(L)",
        atoms=[
            _atom(occupancy=1.0),
            _atom(occupancy=1.0),
        ],
    )

    result = evaluate_q002_occupancy(
        chain,
        minimum_occupancy=1.0,
        reject_alternate_locations=True,
    )

    assert result.passed is True
    assert result.reason == "occupancy_and_altloc_requirements_met"


def test_q002_rejects_missing_occupancy() -> None:
    chain = ChainObservation(
        pdb_id="test",
        model_id=1,
        label_chain_id="A",
        atoms=[
            _atom(occupancy=1.0),
            _atom(occupancy=None),
        ],
    )

    result = evaluate_q002_occupancy(
        chain,
        minimum_occupancy=1.0,
        reject_alternate_locations=True,
    )

    assert result.passed is False
    assert result.reason == "missing_occupancy:1"


def test_q002_rejects_occupancy_below_minimum() -> None:
    chain = ChainObservation(
        pdb_id="test",
        model_id=1,
        label_chain_id="A",
        atoms=[
            _atom(occupancy=1.0),
            _atom(occupancy=0.5),
        ],
    )

    result = evaluate_q002_occupancy(
        chain,
        minimum_occupancy=1.0,
        reject_alternate_locations=True,
    )

    assert result.passed is False
    assert result.reason == "occupancy_below_minimum:1"


def test_q002_rejects_alternate_locations() -> None:
    chain = ChainObservation(
        pdb_id="test",
        model_id=1,
        label_chain_id="A",
        atoms=[
            _atom(occupancy=1.0, alt_id=None),
            _atom(occupancy=1.0, alt_id="A"),
        ],
    )

    result = evaluate_q002_occupancy(
        chain,
        minimum_occupancy=1.0,
        reject_alternate_locations=True,
    )

    assert result.passed is False
    assert result.reason == "alternate_locations_present:1"
