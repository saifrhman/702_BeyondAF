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

from pdbclean.quality import (
    evaluate_q003_residue_continuity,
    missing_internal_label_seq_ids,
)


def _atom_at_residue(
    label_seq_id: int | None,
) -> AtomObservation:
    return AtomObservation(
        model_id=1,
        label_chain_id="A",
        auth_chain_id="A",
        entity_id="1",
        label_seq_id=label_seq_id,
        auth_seq_id=(
            str(label_seq_id)
            if label_seq_id is not None
            else None
        ),
        residue_name="ALA",
        atom_name="CA",
        alt_id=None,
        occupancy=1.0,
        x=0.0,
        y=0.0,
        z=0.0,
    )


def test_q003_accepts_consecutive_observed_residues() -> None:
    chain = ChainObservation(
        pdb_id="test",
        model_id=1,
        label_chain_id="A",
        atoms=[
            _atom_at_residue(4),
            _atom_at_residue(5),
            _atom_at_residue(6),
        ],
    )

    result = evaluate_q003_residue_continuity(
        chain,
        allow_internal_gaps=False,
    )

    assert result.passed is True
    assert result.reason == "observed_residues_consecutive"
    assert missing_internal_label_seq_ids(chain) == []


def test_q003_rejects_internal_gap() -> None:
    chain = ChainObservation(
        pdb_id="test",
        model_id=1,
        label_chain_id="A",
        atoms=[
            _atom_at_residue(4),
            _atom_at_residue(5),
            _atom_at_residue(7),
        ],
    )

    result = evaluate_q003_residue_continuity(
        chain,
        allow_internal_gaps=False,
    )

    assert result.passed is False
    assert result.reason == "internal_label_seq_id_gaps:6"
    assert missing_internal_label_seq_ids(chain) == [6]


def test_q003_accepts_internal_gap_when_configured() -> None:
    chain = ChainObservation(
        pdb_id="test",
        model_id=1,
        label_chain_id="A",
        atoms=[
            _atom_at_residue(10),
            _atom_at_residue(12),
        ],
    )

    result = evaluate_q003_residue_continuity(
        chain,
        allow_internal_gaps=True,
    )

    assert result.passed is True
    assert result.reason == "internal_gaps_allowed"


def test_q003_rejects_missing_label_seq_id() -> None:
    chain = ChainObservation(
        pdb_id="test",
        model_id=1,
        label_chain_id="A",
        atoms=[
            _atom_at_residue(1),
            _atom_at_residue(None),
        ],
    )

    result = evaluate_q003_residue_continuity(
        chain,
        allow_internal_gaps=False,
    )

    assert result.passed is False
    assert result.reason == "missing_label_seq_id:1"

from pdbclean.quality import (
    evaluate_q004_backbone_atoms,
    summarize_backbone_atom_issues,
)


def _backbone_atom(
    *,
    label_seq_id: int,
    atom_name: str,
) -> AtomObservation:
    return AtomObservation(
        model_id=1,
        label_chain_id="A",
        auth_chain_id="A",
        entity_id="1",
        label_seq_id=label_seq_id,
        auth_seq_id=str(label_seq_id),
        residue_name="ALA",
        atom_name=atom_name,
        alt_id=None,
        occupancy=1.0,
        x=0.0,
        y=0.0,
        z=0.0,
    )


def test_q004_accepts_exactly_one_backbone_atom_per_residue() -> None:
    chain = ChainObservation(
        pdb_id="test",
        model_id=1,
        label_chain_id="A",
        atoms=[
            _backbone_atom(label_seq_id=1, atom_name="N"),
            _backbone_atom(label_seq_id=1, atom_name="CA"),
            _backbone_atom(label_seq_id=1, atom_name="C"),
            _backbone_atom(label_seq_id=2, atom_name="N"),
            _backbone_atom(label_seq_id=2, atom_name="CA"),
            _backbone_atom(label_seq_id=2, atom_name="C"),
        ],
    )

    result = evaluate_q004_backbone_atoms(
        chain,
        required_atoms=("N", "CA", "C"),
        require_exactly_one=True,
    )

    assert result.passed is True
    assert result.reason == "required_backbone_atoms_present_once"


def test_q004_rejects_missing_backbone_atom() -> None:
    chain = ChainObservation(
        pdb_id="test",
        model_id=1,
        label_chain_id="A",
        atoms=[
            _backbone_atom(label_seq_id=1, atom_name="N"),
            _backbone_atom(label_seq_id=1, atom_name="CA"),
        ],
    )

    issues = summarize_backbone_atom_issues(
        chain,
        required_atoms=("N", "CA", "C"),
    )

    result = evaluate_q004_backbone_atoms(
        chain,
        required_atoms=("N", "CA", "C"),
        require_exactly_one=True,
    )

    assert issues.missing_atom_count == 1
    assert issues.missing_atoms == ("1:C",)
    assert result.passed is False
    assert result.reason == "missing_backbone_atoms:1:1:C"


def test_q004_rejects_duplicate_backbone_atom() -> None:
    chain = ChainObservation(
        pdb_id="test",
        model_id=1,
        label_chain_id="A",
        atoms=[
            _backbone_atom(label_seq_id=1, atom_name="N"),
            _backbone_atom(label_seq_id=1, atom_name="CA"),
            _backbone_atom(label_seq_id=1, atom_name="CA"),
            _backbone_atom(label_seq_id=1, atom_name="C"),
        ],
    )

    issues = summarize_backbone_atom_issues(
        chain,
        required_atoms=("N", "CA", "C"),
    )

    result = evaluate_q004_backbone_atoms(
        chain,
        required_atoms=("N", "CA", "C"),
        require_exactly_one=True,
    )

    assert issues.duplicate_atom_count == 1
    assert issues.duplicate_atoms == ("1:CA:2",)
    assert result.passed is False
    assert result.reason == "duplicate_backbone_atoms:1:1:CA:2"


def test_q004_allows_duplicates_when_not_required_exactly_once() -> None:
    chain = ChainObservation(
        pdb_id="test",
        model_id=1,
        label_chain_id="A",
        atoms=[
            _backbone_atom(label_seq_id=1, atom_name="N"),
            _backbone_atom(label_seq_id=1, atom_name="CA"),
            _backbone_atom(label_seq_id=1, atom_name="CA"),
            _backbone_atom(label_seq_id=1, atom_name="C"),
        ],
    )

    result = evaluate_q004_backbone_atoms(
        chain,
        required_atoms=("N", "CA", "C"),
        require_exactly_one=False,
    )

    assert result.passed is True

from pdbclean.quality import (
    evaluate_q005_backbone_distance,
    minimum_consecutive_backbone_distance,
)


def _backbone_atom_xyz(
    *,
    label_seq_id: int,
    atom_name: str,
    x: float,
    y: float = 0.0,
    z: float = 0.0,
) -> AtomObservation:
    return AtomObservation(
        model_id=1,
        label_chain_id="A",
        auth_chain_id="A",
        entity_id="1",
        label_seq_id=label_seq_id,
        auth_seq_id=str(label_seq_id),
        residue_name="ALA",
        atom_name=atom_name,
        alt_id=None,
        occupancy=1.0,
        x=x,
        y=y,
        z=z,
    )


def test_q005_accepts_normal_backbone_distances() -> None:
    chain = ChainObservation(
        pdb_id="test",
        model_id=1,
        label_chain_id="A",
        atoms=[
            _backbone_atom_xyz(label_seq_id=1, atom_name="N", x=0.0),
            _backbone_atom_xyz(label_seq_id=1, atom_name="CA", x=1.4),
            _backbone_atom_xyz(label_seq_id=1, atom_name="C", x=2.9),
            _backbone_atom_xyz(label_seq_id=2, atom_name="N", x=4.2),
            _backbone_atom_xyz(label_seq_id=2, atom_name="CA", x=5.6),
            _backbone_atom_xyz(label_seq_id=2, atom_name="C", x=7.1),
        ],
    )

    result = evaluate_q005_backbone_distance(
        chain,
        required_atoms=("N", "CA", "C"),
        minimum_distance_angstrom=0.01,
    )

    assert result.passed is True
    assert result.reason == "backbone_distances_meet_minimum"
    assert minimum_consecutive_backbone_distance(
        chain,
        required_atoms=("N", "CA", "C"),
    ) > 1.0


def test_q005_rejects_distance_below_threshold() -> None:
    chain = ChainObservation(
        pdb_id="test",
        model_id=1,
        label_chain_id="A",
        atoms=[
            _backbone_atom_xyz(label_seq_id=1, atom_name="N", x=0.0),
            _backbone_atom_xyz(label_seq_id=1, atom_name="CA", x=0.005),
            _backbone_atom_xyz(label_seq_id=1, atom_name="C", x=1.5),
        ],
    )

    result = evaluate_q005_backbone_distance(
        chain,
        required_atoms=("N", "CA", "C"),
        minimum_distance_angstrom=0.01,
    )

    assert result.passed is False
    assert result.reason.startswith(
        "consecutive_backbone_distance_below_minimum:"
    )


def test_q005_accepts_distance_exactly_at_threshold() -> None:
    chain = ChainObservation(
        pdb_id="test",
        model_id=1,
        label_chain_id="A",
        atoms=[
            _backbone_atom_xyz(label_seq_id=1, atom_name="N", x=0.0),
            _backbone_atom_xyz(label_seq_id=1, atom_name="CA", x=0.01),
            _backbone_atom_xyz(label_seq_id=1, atom_name="C", x=1.5),
        ],
    )

    result = evaluate_q005_backbone_distance(
        chain,
        required_atoms=("N", "CA", "C"),
        minimum_distance_angstrom=0.01,
    )

    assert result.passed is True


def test_q005_rejects_nonunique_backbone() -> None:
    chain = ChainObservation(
        pdb_id="test",
        model_id=1,
        label_chain_id="A",
        atoms=[
            _backbone_atom_xyz(label_seq_id=1, atom_name="N", x=0.0),
            _backbone_atom_xyz(label_seq_id=1, atom_name="CA", x=1.0),
            _backbone_atom_xyz(label_seq_id=1, atom_name="CA", x=1.1),
            _backbone_atom_xyz(label_seq_id=1, atom_name="C", x=2.0),
        ],
    )

    result = evaluate_q005_backbone_distance(
        chain,
        required_atoms=("N", "CA", "C"),
        minimum_distance_angstrom=0.01,
    )

    assert result.passed is False
    assert result.reason == "backbone_not_exactly_one_per_residue"
