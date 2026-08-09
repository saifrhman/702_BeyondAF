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
from pdbclean.quality import (
    evaluate_q002_disorder,
    q002_disordered_residue_ids,
)


def _atom(
    *,
    label_seq_id: int = 1,
    atom_name: str = "CA",
    occupancy: float | None = 1.0,
    occupancy_raw: str | None = "1.00",
    alt_id: str | None = None,
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
        alt_id=alt_id,
        occupancy=occupancy,
        x=0.0,
        y=0.0,
        z=0.0,
        occupancy_raw=occupancy_raw,
    )


def test_q002_accepts_full_backbone_occupancy() -> None:
    chain = ChainObservation(
        pdb_id="test",
        model_id=1,
        label_chain_id="A",
        atoms=[
            _atom(atom_name="N"),
            _atom(atom_name="CA"),
            _atom(atom_name="C"),
        ],
    )

    result = evaluate_q002_disorder(chain)

    assert result.passed is True
    assert result.reason == "no_disordered_backbone_residues"


def test_q002_accepts_dot_occupancy_like_bri() -> None:
    chain = ChainObservation(
        pdb_id="test",
        model_id=1,
        label_chain_id="A",
        atoms=[
            _atom(
                atom_name="N",
                occupancy=None,
                occupancy_raw=".",
            ),
        ],
    )

    result = evaluate_q002_disorder(chain)

    assert result.passed is True


def test_q002_rejects_question_mark_occupancy() -> None:
    chain = ChainObservation(
        pdb_id="test",
        model_id=1,
        label_chain_id="A",
        atoms=[
            _atom(
                atom_name="CA",
                occupancy=None,
                occupancy_raw="?",
            ),
        ],
    )

    result = evaluate_q002_disorder(chain)

    assert result.passed is False
    assert result.reason == "disordered_backbone_residues:1"


def test_q002_rejects_partial_backbone_occupancy() -> None:
    chain = ChainObservation(
        pdb_id="test",
        model_id=1,
        label_chain_id="A",
        atoms=[
            _atom(
                atom_name="CA",
                occupancy=0.5,
                occupancy_raw="0.50",
            ),
        ],
    )

    result = evaluate_q002_disorder(chain)

    assert result.passed is False
    assert q002_disordered_residue_ids(chain) == {1}


def test_q002_does_not_reject_single_altloc_row() -> None:
    """BRI has no independent label_alt_id rejection."""

    chain = ChainObservation(
        pdb_id="test",
        model_id=1,
        label_chain_id="A",
        atoms=[
            _atom(
                atom_name="CA",
                alt_id="A",
                occupancy=1.0,
                occupancy_raw="1.00",
            ),
        ],
    )

    result = evaluate_q002_disorder(chain)

    assert result.passed is True


def test_q002_rejects_duplicate_backbone_atom() -> None:
    """Duplicate residue/model/chain/atom rows are disorder in BRI."""

    chain = ChainObservation(
        pdb_id="test",
        model_id=1,
        label_chain_id="A",
        atoms=[
            _atom(
                label_seq_id=7,
                atom_name="CA",
                alt_id="A",
            ),
            _atom(
                label_seq_id=7,
                atom_name="CA",
                alt_id="B",
            ),
        ],
    )

    result = evaluate_q002_disorder(chain)

    assert result.passed is False
    assert q002_disordered_residue_ids(chain) == {7}


def test_q002_ignores_sidechain_partial_occupancy() -> None:
    """Protocol 3.2 disorder_check receives N/CA/C features only."""

    chain = ChainObservation(
        pdb_id="test",
        model_id=1,
        label_chain_id="A",
        atoms=[
            _atom(atom_name="N"),
            _atom(atom_name="CA"),
            _atom(atom_name="C"),
            _atom(
                atom_name="CB",
                occupancy=0.5,
                occupancy_raw="0.50",
            ),
        ],
    )

    result = evaluate_q002_disorder(chain)

    assert result.passed is True


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
    minimum_protocol32_backbone_distance,
    protocol32_backbone_distances,
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


def test_q005_two_residues_have_nine_protocol32_distances() -> None:
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

    distances = protocol32_backbone_distances(
        chain,
        required_atoms=("N", "CA", "C"),
    )

    assert len(distances) == 9

    assert {
        (item.atom1, item.atom2)
        for item in distances
    } == {
        ("N", "CA"),
        ("N", "C"),
        ("CA", "C"),
        ("N", "N+1"),
        ("CA", "N+1"),
        ("C", "N+1"),
    }


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
    assert result.reason == "backbone_clash_threshold_met"

    minimum = minimum_protocol32_backbone_distance(
        chain,
        required_atoms=("N", "CA", "C"),
    )
    assert minimum is not None
    assert minimum > 1.0


def test_q005_detects_same_residue_n_c_clash() -> None:
    """N-C is one of the BRI comparisons missed by the old implementation."""

    chain = ChainObservation(
        pdb_id="test",
        model_id=1,
        label_chain_id="A",
        atoms=[
            _backbone_atom_xyz(label_seq_id=1, atom_name="N", x=0.0),
            _backbone_atom_xyz(label_seq_id=1, atom_name="CA", x=1.4),
            _backbone_atom_xyz(label_seq_id=1, atom_name="C", x=0.005),
        ],
    )

    result = evaluate_q005_backbone_distance(
        chain,
        required_atoms=("N", "CA", "C"),
        minimum_distance_angstrom=0.01,
    )

    assert result.passed is False
    assert result.reason.startswith("backbone_clashes:1:")


def test_q005_detects_inter_residue_n_n_clash() -> None:
    """N(i)-N(i+1) is another BRI comparison absent from the old code."""

    chain = ChainObservation(
        pdb_id="test",
        model_id=1,
        label_chain_id="A",
        atoms=[
            _backbone_atom_xyz(label_seq_id=1, atom_name="N", x=0.0),
            _backbone_atom_xyz(label_seq_id=1, atom_name="CA", x=1.4),
            _backbone_atom_xyz(label_seq_id=1, atom_name="C", x=2.9),
            _backbone_atom_xyz(label_seq_id=2, atom_name="N", x=0.005),
            _backbone_atom_xyz(label_seq_id=2, atom_name="CA", x=4.2),
            _backbone_atom_xyz(label_seq_id=2, atom_name="C", x=5.7),
        ],
    )

    result = evaluate_q005_backbone_distance(
        chain,
        required_atoms=("N", "CA", "C"),
        minimum_distance_angstrom=0.01,
    )

    assert result.passed is False
    assert result.reason.startswith("backbone_clashes:1:")


def test_q005_detects_inter_residue_ca_n_clash() -> None:
    """CA(i)-N(i+1) is another comparison missed by the old code."""

    chain = ChainObservation(
        pdb_id="test",
        model_id=1,
        label_chain_id="A",
        atoms=[
            _backbone_atom_xyz(label_seq_id=1, atom_name="N", x=0.0),
            _backbone_atom_xyz(label_seq_id=1, atom_name="CA", x=1.4),
            _backbone_atom_xyz(label_seq_id=1, atom_name="C", x=2.9),
            _backbone_atom_xyz(label_seq_id=2, atom_name="N", x=1.405),
            _backbone_atom_xyz(label_seq_id=2, atom_name="CA", x=4.2),
            _backbone_atom_xyz(label_seq_id=2, atom_name="C", x=5.7),
        ],
    )

    result = evaluate_q005_backbone_distance(
        chain,
        required_atoms=("N", "CA", "C"),
        minimum_distance_angstrom=0.01,
    )

    assert result.passed is False
    assert result.reason.startswith("backbone_clashes:1:")


def test_q005_threshold_is_strict() -> None:
    """BRI uses distance < lower_bound, so exactly 0.01 Å is accepted."""

    chain = ChainObservation(
        pdb_id="test",
        model_id=1,
        label_chain_id="A",
        atoms=[
            _backbone_atom_xyz(label_seq_id=1, atom_name="N", x=0.0),
            _backbone_atom_xyz(label_seq_id=1, atom_name="CA", x=1.4),
            _backbone_atom_xyz(label_seq_id=1, atom_name="C", x=0.01),
        ],
    )

    result = evaluate_q005_backbone_distance(
        chain,
        required_atoms=("N", "CA", "C"),
        minimum_distance_angstrom=0.01,
    )

    assert result.passed is True


def test_q005_rejects_nonunique_backbone_precondition() -> None:
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
    assert result.reason.startswith("q005_precondition_failed:")

