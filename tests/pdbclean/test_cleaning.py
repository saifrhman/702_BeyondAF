"""Tests for stateful Protocol 3.2 cleaning orchestration."""

from pdbclean.cleaning import clean_protocol32_chain, is_protocol32_candidate
from pdbclean.mmcif_parser import AtomObservation, ChainObservation


def _atom(
    *,
    atom_name: str,
    group_pdb: str = "ATOM",
    residue_id: int = 1,
    x: float | None = None,
) -> AtomObservation:
    if x is None:
        atom_offset = {
            "N": 0.0,
            "CA": 1.0,
            "C": 2.0,
        }.get(atom_name, 3.0)
        x = (residue_id - 1) * 10.0 + atom_offset

    return AtomObservation(
        model_id=1,
        label_chain_id="A",
        auth_chain_id="A",
        entity_id="1",
        label_seq_id=residue_id,
        auth_seq_id=str(residue_id),
        residue_name="ALA",
        atom_name=atom_name,
        alt_id=None,
        occupancy=1.0,
        x=x,
        y=0.0,
        z=0.0,
        group_pdb=group_pdb,
        occupancy_raw="1.00",
    )


def _chain(atoms: list[AtomObservation]) -> ChainObservation:
    return ChainObservation(
        pdb_id="test",
        model_id=1,
        label_chain_id="A",
        entry_has_polypeptide=True,
        atoms=atoms,
    )


def test_protocol32_candidate_requires_nonempty_atom_backbone_projection() -> None:
    assert is_protocol32_candidate(
        _chain([_atom(atom_name="CA")])
    ) is True


def test_hetatm_only_chain_is_not_protocol32_candidate() -> None:
    assert is_protocol32_candidate(
        _chain([_atom(atom_name="CA", group_pdb="HETATM")])
    ) is False


def test_atom_chain_without_n_ca_c_is_not_protocol32_candidate() -> None:
    assert is_protocol32_candidate(
        _chain([_atom(atom_name="CB")])
    ) is False


def test_cleaner_records_empty_projection_as_non_candidate() -> None:
    result = clean_protocol32_chain(
        _chain([_atom(atom_name="CB")])
    )

    assert result.status == "non_candidate"
    assert result.reason == "empty_protocol32_backbone_projection"
    assert result.terminal_stage == "candidate_selection"
    assert result.retained_chain is None


def test_cleaner_applies_q001_after_candidate_selection() -> None:
    chain = _chain([_atom(atom_name="CA")])
    chain.entry_has_polypeptide = False

    result = clean_protocol32_chain(chain)

    assert result.status == "rejected"
    assert result.reason == "entry_contains_no_polypeptide_entity"
    assert result.terminal_stage == "Q001"
    assert result.retained_chain is None


def test_cleaner_preserves_projected_working_chain_through_q006_q004() -> None:
    chain = _chain(
        [
            _atom(atom_name="N"),
            _atom(atom_name="CA"),
            _atom(atom_name="C"),
            _atom(atom_name="CB"),
            _atom(atom_name="CA", group_pdb="HETATM"),
        ]
    )

    result = clean_protocol32_chain(chain)

    assert result.status == "accepted"
    assert result.reason == "protocol32_clean"
    assert result.terminal_stage == "completed"
    assert result.retained_chain is not None
    assert [atom.atom_name for atom in result.retained_chain.atoms] == [
        "N",
        "CA",
        "C",
    ]
    assert all(
        atom.group_pdb == "ATOM"
        for atom in result.retained_chain.atoms
    )


def test_q002_terminal_dirty_residue_is_trimmed_and_chain_survives() -> None:
    atoms = [
        _atom(atom_name=name, residue_id=residue_id)
        for residue_id in (1, 2, 3)
        for name in ("N", "CA", "C")
    ]

    # Make terminal residue 3 disordered.
    atoms[-2] = AtomObservation(
        **{
            **atoms[-2].__dict__,
            "occupancy": 0.5,
            "occupancy_raw": "0.50",
        }
    )

    result = clean_protocol32_chain(_chain(atoms))

    assert result.status == "accepted"
    assert result.reason == "protocol32_clean"
    assert result.terminal_stage == "completed"
    assert result.retained_chain is not None
    assert {
        atom.label_seq_id
        for atom in result.retained_chain.atoms
    } == {1, 2}

    assert len(result.dirty_residues) == 1
    assert result.dirty_residues[0].residue_id == 3
    assert result.dirty_residues[0].rule_id == "Q002"
    assert result.dirty_residues[0].dirty_type == "disordered"


def test_q002_internal_dirty_residue_creates_q003_chain_break() -> None:
    atoms = [
        _atom(atom_name=name, residue_id=residue_id)
        for residue_id in (1, 2, 3)
        for name in ("N", "CA", "C")
    ]

    # Make internal residue 2 disordered.
    index = next(
        i
        for i, atom in enumerate(atoms)
        if atom.label_seq_id == 2 and atom.atom_name == "CA"
    )
    atoms[index] = AtomObservation(
        **{
            **atoms[index].__dict__,
            "occupancy": 0.5,
            "occupancy_raw": "0.50",
        }
    )

    result = clean_protocol32_chain(_chain(atoms))

    assert result.status == "rejected"
    assert result.terminal_stage == "Q003_after_Q002"
    assert result.reason == "internal_label_seq_id_gaps:2"
    assert result.missing_label_seq_ids == (2,)
    assert [item.residue_id for item in result.dirty_residues] == [2]


def test_q002_can_remove_entire_working_chain() -> None:
    atoms = [
        _atom(atom_name=name, residue_id=1)
        for name in ("N", "CA", "C")
    ]

    atoms[0] = AtomObservation(
        **{
            **atoms[0].__dict__,
            "occupancy": 0.5,
            "occupancy_raw": "0.50",
        }
    )

    result = clean_protocol32_chain(_chain(atoms))

    assert result.status == "rejected"
    assert result.reason == "all_working_residues_removed"
    assert result.terminal_stage == "Q002"
    assert result.retained_chain is None
    assert [item.residue_id for item in result.dirty_residues] == [1]


def test_q006_nonstandard_wins_over_q004_for_same_residue() -> None:
    atoms = [
        _atom(atom_name=name, residue_id=residue_id)
        for residue_id in (1, 2)
        for name in ("N", "CA", "C")
    ]

    # Residue 3 is SEC -> U and is also incomplete because C is absent.
    for name in ("N", "CA"):
        atom = _atom(atom_name=name, residue_id=3)
        atoms.append(
            AtomObservation(
                **{
                    **atom.__dict__,
                    "residue_name": "SEC",
                }
            )
        )

    result = clean_protocol32_chain(_chain(atoms))

    assert result.status == "accepted"
    assert result.reason == "protocol32_clean"
    assert result.terminal_stage == "completed"
    assert result.retained_chain is not None

    records = [
        item
        for item in result.dirty_residues
        if item.residue_id == 3
    ]

    assert len(records) == 1
    assert records[0].rule_id == "Q006"
    assert records[0].dirty_type == "non-standard"
    assert records[0].deposited_residue_name == "SEC"
    assert records[0].mapped_residue_code == "U"


def test_q004_terminal_incomplete_residue_is_trimmed() -> None:
    atoms = [
        _atom(atom_name=name, residue_id=residue_id)
        for residue_id in (1, 2)
        for name in ("N", "CA", "C")
    ]

    # Canonical ALA residue 3 is incomplete: C is absent.
    atoms.extend(
        [
            _atom(atom_name="N", residue_id=3),
            _atom(atom_name="CA", residue_id=3),
        ]
    )

    result = clean_protocol32_chain(_chain(atoms))

    assert result.status == "accepted"
    assert result.reason == "protocol32_clean"
    assert result.terminal_stage == "completed"
    assert result.retained_chain is not None
    assert {
        atom.label_seq_id
        for atom in result.retained_chain.atoms
    } == {1, 2}

    records = [
        item
        for item in result.dirty_residues
        if item.residue_id == 3
    ]

    assert len(records) == 1
    assert records[0].rule_id == "Q004"
    assert records[0].dirty_type == "incomplete"
    assert records[0].mapped_residue_code == "A"


def test_q006_internal_nonstandard_residue_creates_chain_break() -> None:
    atoms = []

    for residue_id in (1, 2, 3):
        for name in ("N", "CA", "C"):
            atom = _atom(
                atom_name=name,
                residue_id=residue_id,
            )

            if residue_id == 2:
                atom = AtomObservation(
                    **{
                        **atom.__dict__,
                        "residue_name": "SEC",
                    }
                )

            atoms.append(atom)

    result = clean_protocol32_chain(_chain(atoms))

    assert result.status == "rejected"
    assert result.terminal_stage == "Q003_after_Q006_Q004"
    assert result.reason == "internal_label_seq_id_gaps:2"
    assert result.missing_label_seq_ids == (2,)

    records = [
        item
        for item in result.dirty_residues
        if item.residue_id == 2
    ]

    assert len(records) == 1
    assert records[0].rule_id == "Q006"


def test_q004_can_remove_entire_remaining_working_chain() -> None:
    atoms = [
        _atom(atom_name="N", residue_id=1),
        _atom(atom_name="CA", residue_id=1),
    ]

    result = clean_protocol32_chain(_chain(atoms))

    assert result.status == "rejected"
    assert result.reason == "all_working_residues_removed"
    assert result.terminal_stage == "Q006_Q004"
    assert result.retained_chain is None

    assert len(result.dirty_residues) == 1
    assert result.dirty_residues[0].rule_id == "Q004"
    assert result.dirty_residues[0].residue_id == 1


def test_q005_terminal_clash_residue_is_trimmed_and_chain_survives() -> None:
    atoms = [
        _atom(atom_name=name, residue_id=residue_id)
        for residue_id in (1, 2, 3)
        for name in ("N", "CA", "C")
    ]

    # Create a strict N-CA clash in terminal residue 3.
    for index, atom in enumerate(atoms):
        if atom.label_seq_id == 3 and atom.atom_name == "N":
            atoms[index] = AtomObservation(
                **{
                    **atom.__dict__,
                    "x": 20.0,
                }
            )
        elif atom.label_seq_id == 3 and atom.atom_name == "CA":
            atoms[index] = AtomObservation(
                **{
                    **atom.__dict__,
                    "x": 20.001,
                }
            )

    result = clean_protocol32_chain(_chain(atoms))

    assert result.status == "accepted"
    assert result.reason == "protocol32_clean"
    assert result.terminal_stage == "completed"
    assert result.retained_chain is not None
    assert {
        atom.label_seq_id
        for atom in result.retained_chain.atoms
    } == {1, 2}

    records = [
        item
        for item in result.dirty_residues
        if item.residue_id == 3
    ]

    assert len(records) == 1
    assert records[0].rule_id == "Q005"
    assert records[0].dirty_type == "clash"


def test_q005_internal_clash_residue_creates_final_chain_break() -> None:
    atoms = [
        _atom(atom_name=name, residue_id=residue_id)
        for residue_id in (1, 2, 3)
        for name in ("N", "CA", "C")
    ]

    # Create a strict N-CA clash in internal residue 2.
    for index, atom in enumerate(atoms):
        if atom.label_seq_id == 2 and atom.atom_name == "N":
            atoms[index] = AtomObservation(
                **{
                    **atom.__dict__,
                    "x": 10.0,
                }
            )
        elif atom.label_seq_id == 2 and atom.atom_name == "CA":
            atoms[index] = AtomObservation(
                **{
                    **atom.__dict__,
                    "x": 10.001,
                }
            )

    result = clean_protocol32_chain(_chain(atoms))

    assert result.status == "rejected"
    assert result.terminal_stage == "Q003_final"
    assert result.reason == "internal_label_seq_id_gaps:2"
    assert result.missing_label_seq_ids == (2,)

    records = [
        item
        for item in result.dirty_residues
        if item.residue_id == 2
    ]

    assert len(records) == 1
    assert records[0].rule_id == "Q005"


def test_q005_exact_threshold_is_not_a_clash() -> None:
    atoms = [
        _atom(atom_name=name, residue_id=1)
        for name in ("N", "CA", "C")
    ]

    for index, atom in enumerate(atoms):
        if atom.atom_name == "N":
            atoms[index] = AtomObservation(
                **{
                    **atom.__dict__,
                    "x": 0.0,
                }
            )
        elif atom.atom_name == "CA":
            atoms[index] = AtomObservation(
                **{
                    **atom.__dict__,
                    "x": 0.01,
                }
            )

    result = clean_protocol32_chain(_chain(atoms))

    assert result.status == "accepted"
    assert result.reason == "protocol32_clean"
    assert result.terminal_stage == "completed"
    assert not [
        item
        for item in result.dirty_residues
        if item.rule_id == "Q005"
    ]

