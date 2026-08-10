"""Differential tests against pinned BRI v1.2.2 cleaning behaviour."""

import gzip
import hashlib
import tempfile
from pathlib import Path

import pandas as pd

from bri.filter import integrated_chainwise_filter
from bri.pdbx2df import Entry

from pdbclean.cleaning import clean_protocol32_chain
from pdbclean.mmcif_parser import (
    AtomObservation,
    ChainObservation,
    parse_coordinate_mmcif_bytes,
)
from pdbclean.quality import protocol32_backbone_projection


def _atom(
    residue_id: int,
    atom_name: str,
    *,
    residue_name: str = "ALA",
    occupancy_raw: str = "1.00",
    x: float | None = None,
) -> AtomObservation:
    if x is None:
        offset = {
            "N": 0.0,
            "CA": 1.0,
            "C": 2.0,
        }[atom_name]
        x = (residue_id - 1) * 10.0 + offset

    occupancy = None
    try:
        occupancy = float(occupancy_raw)
    except ValueError:
        pass

    return AtomObservation(
        model_id=1,
        label_chain_id="A",
        auth_chain_id="A",
        entity_id="1",
        label_seq_id=residue_id,
        auth_seq_id=str(residue_id),
        residue_name=residue_name,
        atom_name=atom_name,
        alt_id=None,
        occupancy=occupancy,
        x=x,
        y=0.0,
        z=0.0,
        group_pdb="ATOM",
        occupancy_raw=occupancy_raw,
    )


def _chain(
    atoms: list[AtomObservation],
) -> ChainObservation:
    return ChainObservation(
        pdb_id="test",
        model_id=1,
        label_chain_id="A",
        auth_chain_id="A",
        entity_id="1",
        polymer_type="polypeptide(L)",
        entry_has_polypeptide=True,
        atoms=atoms,
    )


def _complete_chain(
    residue_ids: tuple[int, ...],
) -> ChainObservation:
    return _chain(
        [
            _atom(residue_id, atom_name)
            for residue_id in residue_ids
            for atom_name in ("N", "CA", "C")
        ]
    )


def _to_bri_dataframe(
    chain: ChainObservation,
) -> pd.DataFrame:
    """Construct the dataframe supplied to BRI integrated_chainwise_filter."""

    return pd.DataFrame(
        [
            {
                "model_id": atom.model_id,
                "chain_id": atom.label_chain_id,
                "residue_id": atom.label_seq_id,
                "auth_residue_id": int(atom.auth_seq_id),
                "residue_label": atom.residue_name,
                "atom": atom.atom_name,
                "x": atom.x,
                "y": atom.y,
                "z": atom.z,
                "occupancy": atom.occupancy_raw,
            }
            for atom in chain.atoms
            if atom.group_pdb == "ATOM"
            and atom.atom_name in {"N", "CA", "C"}
        ]
    )


def _run_bri(
    chain: ChainObservation,
) -> pd.DataFrame:
    result = integrated_chainwise_filter(
        _to_bri_dataframe(chain).copy(deep=True)
    )

    assert result is not None
    return result


def _bri_clean_residue_ids(
    result: pd.DataFrame,
) -> tuple[int, ...]:
    if "type" not in result.columns:
        return ()

    clean = result[result["type"] == "clean"]

    return tuple(
        sorted(
            {
                int(value)
                for value in clean["residue_id"]
            }
        )
    )


def _our_clean_residue_ids(
    chain: ChainObservation,
) -> tuple[int, ...]:
    result = clean_protocol32_chain(chain)

    if result.retained_chain is None:
        return ()

    return tuple(
        sorted(
            {
                int(atom.label_seq_id)
                for atom in result.retained_chain.atoms
                if atom.label_seq_id is not None
            }
        )
    )


def _bri_dirty_pairs(
    result: pd.DataFrame,
) -> set[tuple[int, str]]:
    if "type" not in result.columns:
        return set()

    dirty = result[result["type"] != "clean"]

    return {
        (int(row.residue_id), str(row.type))
        for row in dirty.itertuples()
    }


def _our_dirty_pairs(
    chain: ChainObservation,
) -> set[tuple[int, str]]:
    result = clean_protocol32_chain(chain)

    return {
        (record.residue_id, record.dirty_type)
        for record in result.dirty_residues
    }


def test_clean_chain_matches_bri() -> None:
    chain = _complete_chain((1, 2, 3))

    bri_result = _run_bri(chain)

    assert _our_clean_residue_ids(chain) == _bri_clean_residue_ids(
        bri_result
    )
    assert _our_clean_residue_ids(chain) == (1, 2, 3)

    assert _our_dirty_pairs(chain) == _bri_dirty_pairs(bri_result)
    assert _our_dirty_pairs(chain) == set()


def test_terminal_q002_disorder_matches_bri() -> None:
    chain = _complete_chain((1, 2, 3))

    atoms = list(chain.atoms)

    target = next(
        index
        for index, atom in enumerate(atoms)
        if atom.label_seq_id == 3
        and atom.atom_name == "CA"
    )

    original = atoms[target]

    atoms[target] = AtomObservation(
        **{
            **original.__dict__,
            "occupancy": 0.5,
            "occupancy_raw": "0.50",
        }
    )

    chain = _chain(atoms)
    bri_result = _run_bri(chain)

    assert _our_clean_residue_ids(chain) == _bri_clean_residue_ids(
        bri_result
    )
    assert _our_clean_residue_ids(chain) == (1, 2)

    assert _our_dirty_pairs(chain) == _bri_dirty_pairs(bri_result)
    assert _our_dirty_pairs(chain) == {(3, "disordered")}


def _bri_terminal_status(
    result: pd.DataFrame,
) -> str:
    types = set(result["type"].astype(str))

    if "clean" in types:
        return "accepted"

    return "rejected"


def _bri_has_chain_break(
    result: pd.DataFrame,
) -> bool:
    return bool(
        (result["type"].astype(str) == "chain-break").any()
    )


def _bri_protocol_dirty_pairs(
    result: pd.DataFrame,
) -> set[tuple[int, str]]:
    """Compare residue-level dirty types, excluding Q003 chain-break rows."""

    dirty_types = {
        "disordered",
        "non-standard",
        "incomplete",
        "clash",
    }

    return {
        (int(row.residue_id), str(row.type))
        for row in result.itertuples()
        if str(row.type) in dirty_types
    }


def test_internal_q002_disorder_chain_break_matches_bri() -> None:
    chain = _complete_chain((1, 2, 3))
    atoms = list(chain.atoms)

    target = next(
        index
        for index, atom in enumerate(atoms)
        if atom.label_seq_id == 2
        and atom.atom_name == "CA"
    )

    original = atoms[target]
    atoms[target] = AtomObservation(
        **{
            **original.__dict__,
            "occupancy": 0.5,
            "occupancy_raw": "0.50",
        }
    )

    chain = _chain(atoms)

    ours = clean_protocol32_chain(chain)
    bri_result = _run_bri(chain)

    assert ours.status == _bri_terminal_status(bri_result)
    assert ours.status == "rejected"

    assert _bri_has_chain_break(bri_result) is True
    assert ours.terminal_stage == "Q003_after_Q002"
    assert ours.missing_label_seq_ids == (2,)

    assert _our_clean_residue_ids(chain) == ()
    assert _bri_clean_residue_ids(bri_result) == ()

    assert _our_dirty_pairs(chain) == _bri_protocol_dirty_pairs(
        bri_result
    )
    assert _our_dirty_pairs(chain) == {(2, "disordered")}


def test_terminal_q004_incomplete_matches_bri() -> None:
    atoms = [
        _atom(residue_id, atom_name)
        for residue_id in (1, 2)
        for atom_name in ("N", "CA", "C")
    ]

    atoms.extend(
        [
            _atom(3, "N"),
            _atom(3, "CA"),
        ]
    )

    chain = _chain(atoms)

    ours = clean_protocol32_chain(chain)
    bri_result = _run_bri(chain)

    assert ours.status == _bri_terminal_status(bri_result)
    assert ours.status == "accepted"

    assert _our_clean_residue_ids(chain) == _bri_clean_residue_ids(
        bri_result
    )
    assert _our_clean_residue_ids(chain) == (1, 2)

    assert _our_dirty_pairs(chain) == _bri_protocol_dirty_pairs(
        bri_result
    )
    assert _our_dirty_pairs(chain) == {(3, "incomplete")}


def test_q006_precedence_over_q004_matches_bri() -> None:
    atoms = [
        _atom(residue_id, atom_name)
        for residue_id in (1, 2)
        for atom_name in ("N", "CA", "C")
    ]

    atoms.extend(
        [
            _atom(3, "N", residue_name="SEC"),
            _atom(3, "CA", residue_name="SEC"),
        ]
    )

    chain = _chain(atoms)

    ours = clean_protocol32_chain(chain)
    bri_result = _run_bri(chain)

    assert ours.status == _bri_terminal_status(bri_result)
    assert ours.status == "accepted"

    assert _our_clean_residue_ids(chain) == _bri_clean_residue_ids(
        bri_result
    )
    assert _our_clean_residue_ids(chain) == (1, 2)

    assert _our_dirty_pairs(chain) == _bri_protocol_dirty_pairs(
        bri_result
    )
    assert _our_dirty_pairs(chain) == {(3, "non-standard")}


def test_llp_mapping_matches_pinned_bri() -> None:
    atoms = []

    for residue_id, residue_name in (
        (1, "ALA"),
        (2, "LLP"),
        (3, "GLY"),
    ):
        for atom_name in ("N", "CA", "C"):
            atoms.append(
                _atom(
                    residue_id,
                    atom_name,
                    residue_name=residue_name,
                )
            )

    chain = _chain(atoms)

    ours = clean_protocol32_chain(chain)
    bri_result = _run_bri(chain)

    assert ours.status == _bri_terminal_status(bri_result)
    assert ours.status == "accepted"

    assert _our_clean_residue_ids(chain) == _bri_clean_residue_ids(
        bri_result
    )
    assert _our_clean_residue_ids(chain) == (1, 2, 3)

    assert _our_dirty_pairs(chain) == _bri_protocol_dirty_pairs(
        bri_result
    )
    assert _our_dirty_pairs(chain) == set()


def test_terminal_q005_clash_matches_bri() -> None:
    chain = _complete_chain((1, 2, 3))
    atoms = list(chain.atoms)

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

    chain = _chain(atoms)

    ours = clean_protocol32_chain(chain)
    bri_result = _run_bri(chain)

    assert ours.status == _bri_terminal_status(bri_result)
    assert ours.status == "accepted"

    assert _our_clean_residue_ids(chain) == _bri_clean_residue_ids(
        bri_result
    )
    assert _our_clean_residue_ids(chain) == (1, 2)

    assert _our_dirty_pairs(chain) == _bri_protocol_dirty_pairs(
        bri_result
    )
    assert _our_dirty_pairs(chain) == {(3, "clash")}


def test_internal_q005_clash_chain_break_matches_bri() -> None:
    chain = _complete_chain((1, 2, 3))
    atoms = list(chain.atoms)

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

    chain = _chain(atoms)

    ours = clean_protocol32_chain(chain)
    bri_result = _run_bri(chain)

    assert ours.status == _bri_terminal_status(bri_result)
    assert ours.status == "rejected"

    assert _bri_has_chain_break(bri_result) is True
    assert ours.terminal_stage == "Q003_final"
    assert ours.missing_label_seq_ids == (2,)

    assert _our_dirty_pairs(chain) == _bri_protocol_dirty_pairs(
        bri_result
    )
    assert _our_dirty_pairs(chain) == {(2, "clash")}


def test_q005_exact_threshold_matches_bri() -> None:
    chain = _complete_chain((1,))
    atoms = list(chain.atoms)

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

    chain = _chain(atoms)

    ours = clean_protocol32_chain(chain)
    bri_result = _run_bri(chain)

    assert ours.status == _bri_terminal_status(bri_result)
    assert ours.status == "accepted"

    assert _our_dirty_pairs(chain) == _bri_protocol_dirty_pairs(
        bri_result
    )
    assert _our_dirty_pairs(chain) == set()


def test_preexisting_internal_q003_gap_matches_bri() -> None:
    chain = _complete_chain((1, 3))

    ours = clean_protocol32_chain(chain)
    bri_result = _run_bri(chain)

    assert ours.status == _bri_terminal_status(bri_result)
    assert ours.status == "rejected"

    assert _bri_has_chain_break(bri_result) is True
    assert ours.terminal_stage == "Q003_after_Q002"
    assert ours.missing_label_seq_ids == (2,)

    assert _our_clean_residue_ids(chain) == ()
    assert _bri_clean_residue_ids(bri_result) == ()

def test_1aam_revision_projection_and_outcome_match_pinned_bri() -> None:
    """PDB revision, not special-case logic, changes the 1AAM-A outcome."""

    fixture_root = (
        Path(__file__).parent
        / "fixtures"
        / "1aam_revisions"
    )

    cases = (
        (
            "v1.3",
            "pdb_00001aam_xyz_v1-3.cif.gz",
            "bfe37d729a1feac8d4b2d8a57b7d1d8d088bd8f34d0817c5788699042ec1b012",
            396,
            True,
            "accepted",
        ),
        (
            "v2.0",
            "pdb_00001aam_xyz_v2-0.cif.gz",
            "57dd5b1c7bd03b38b1fcc3348805b553d762c3dfb0d2bd30e3aa13678a42d06f",
            395,
            False,
            "rejected",
        ),
    )

    for (
        version,
        filename,
        expected_sha256,
        expected_residue_count,
        residue_246_present,
        expected_status,
    ) in cases:
        compressed = (fixture_root / filename).read_bytes()

        assert hashlib.sha256(compressed).hexdigest() == expected_sha256

        # Pinned BRI parses and projects the deposited revision itself.
        with tempfile.TemporaryDirectory() as tmp:
            cif_path = Path(tmp) / "1aam.cif"
            cif_path.write_bytes(gzip.decompress(compressed))

            bri_entry = Entry(str(cif_path))
            assert len(bri_entry.chains) == 1

            bri_features = bri_entry.chains[0].get_feature(
                "features",
                HETATM=False,
            ).copy()

            bri_ids = tuple(
                sorted(
                    {
                        int(value)
                        for value in bri_features["residue_id"]
                    }
                )
            )

            bri_result = integrated_chainwise_filter(
                bri_features.copy(deep=True)
            )

        assert bri_result is not None

        # PDBClean parses the exact same compressed mmCIF revision.
        parsed = parse_coordinate_mmcif_bytes(
            compressed,
            pdb_id="1aam",
        )

        our_chain = next(
            chain
            for chain in parsed
            if chain.model_id == 1
            and chain.label_chain_id == "A"
        )

        projected = protocol32_backbone_projection(our_chain)

        our_ids = tuple(
            sorted(
                {
                    int(atom.label_seq_id)
                    for atom in projected.atoms
                    if atom.label_seq_id is not None
                }
            )
        )

        ours = clean_protocol32_chain(our_chain)

        # Both implementations must see exactly the same backbone input.
        assert our_ids == bri_ids
        assert len(our_ids) == expected_residue_count
        assert (246 in our_ids) is residue_246_present

        # Both implementations must reach the same final decision.
        assert ours.status == expected_status
        assert ours.status == _bri_terminal_status(bri_result)

        if version == "v1.3":
            assert our_ids == tuple(range(1, 397))
            assert _our_clean_residue_ids(our_chain) == tuple(range(1, 397))
            assert _bri_clean_residue_ids(bri_result) == tuple(range(1, 397))

        if version == "v2.0":
            assert _bri_has_chain_break(bri_result) is True
            assert ours.terminal_stage == "Q003_after_Q002"
            assert ours.missing_label_seq_ids == (246,)
            assert _our_clean_residue_ids(our_chain) == ()
            assert _bri_clean_residue_ids(bri_result) == ()
