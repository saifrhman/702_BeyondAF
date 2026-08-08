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
