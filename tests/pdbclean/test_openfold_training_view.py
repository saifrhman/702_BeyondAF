"""The training-view projection decides which residues OpenFold ever sees.

Every failure mode here is silent at run time.  Feed OpenFold the deposited
polymer and it builds structure features one length and MSA features another;
shift the numbering origin by one and every residue is paired with the wrong
coordinates while all the lengths stay right.  Neither raises.

So these tests pin the arithmetic that places Gold's ``label_seq_id`` values on
OpenFold's seqres index, and assert the projection refuses rather than guesses
whenever the lineage does not fit the deposited sequence.

``project_mmcif_object`` imports OpenFold lazily; OpenFold is not installed in
the pdbclean environment, so a stub stands in for the one dataclass it needs.
The stub is a faithful stand-in: the projection only ever constructs an
MmcifObject, never interprets one.
"""

from __future__ import annotations

import sys
import types
from dataclasses import dataclass
from typing import Any, Mapping

import pytest

from pdbclean.openfold_training_view import (
    ProjectionError,
    contiguous_ranges,
    expand_ranges,
    project_mmcif_object,
    projected_indices,
    projection_entry,
    seqres_start_number,
)


# --------------------------------------------------------------------------
# a stand-in for OpenFold's MmcifObject
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class StubMmcifObject:
    file_id: str
    header: Any
    structure: Any
    chain_to_seqres: Mapping[str, str]
    seqres_to_structure: Mapping[str, Mapping[int, Any]]
    raw_string: Any


@pytest.fixture(autouse=True)
def stub_openfold(monkeypatch):
    module = types.ModuleType("openfold")
    data = types.ModuleType("openfold.data")
    parsing = types.ModuleType("openfold.data.mmcif_parsing")

    parsing.MmcifObject = StubMmcifObject
    data.mmcif_parsing = parsing
    module.data = data

    monkeypatch.setitem(sys.modules, "openfold", module)
    monkeypatch.setitem(sys.modules, "openfold.data", data)
    monkeypatch.setitem(sys.modules, "openfold.data.mmcif_parsing", parsing)


SEQRES = "MVLSEGEWQLVLHVWAKVEAD"          # 21 residues, label_seq_id 1..21


def make_object(seqres: str = SEQRES, start_num: int = 1, chain: str = "A"):
    """A deposited chain whose residue i carries a recognisable marker."""

    return StubMmcifObject(
        file_id="1abc",
        header={"resolution": 1.5, "release_date": "2020-01-01"},
        structure=object(),
        chain_to_seqres={chain: seqres, "B": "AAAA"},
        seqres_to_structure={
            chain: {i: f"res{i}" for i in range(len(seqres))},
            "B": {i: f"other{i}" for i in range(4)},
        },
        raw_string={
            "_entity_poly_seq.num": [
                str(start_num + i) for i in range(len(seqres))
            ],
            "_entity_poly_seq.entity_id": ["1"] * len(seqres),
        },
    )


# --------------------------------------------------------------------------
# run encoding
# --------------------------------------------------------------------------

def test_a_contiguous_lineage_is_one_run():
    assert contiguous_ranges([1, 2, 3, 4]) == [(1, 4)]


def test_an_internal_gap_splits_the_run():
    assert contiguous_ranges([1, 2, 5, 6, 9]) == [(1, 2), (5, 6), (9, 9)]


def test_a_single_residue_is_a_run():
    assert contiguous_ranges([7]) == [(7, 7)]


def test_no_residues_is_no_runs():
    assert contiguous_ranges([]) == []


@pytest.mark.parametrize(
    "residues",
    [[1, 2, 3], [5], [1, 2, 9, 10, 11], list(range(1, 200))],
)
def test_runs_round_trip(residues):
    assert expand_ranges(projection_entry(residues)) == residues


# --------------------------------------------------------------------------
# placing Gold ids on OpenFold's index
# --------------------------------------------------------------------------

def test_numbering_origin_is_read_from_the_deposited_file():
    assert seqres_start_number(make_object(start_num=1), "1") == 1
    assert seqres_start_number(make_object(start_num=7), "1") == 7


def test_a_missing_entity_is_refused_not_defaulted():
    with pytest.raises(ProjectionError, match="entity"):
        seqres_start_number(make_object(), "99")


def test_indices_are_offset_by_the_numbering_origin():
    assert projected_indices([1, 2, 3], 1) == [0, 1, 2]
    assert projected_indices([7, 8, 9], 7) == [0, 1, 2]


def test_duplicate_residues_are_refused():
    with pytest.raises(ProjectionError, match="duplicates"):
        projected_indices([1, 1, 2], 1)


def test_an_empty_lineage_is_refused():
    with pytest.raises(ProjectionError, match="empty"):
        projected_indices([], 1)


# --------------------------------------------------------------------------
# the projection itself
# --------------------------------------------------------------------------

def test_projection_keeps_exactly_the_retained_residues():
    obj = make_object()

    projected = project_mmcif_object(obj, "A", [3, 4, 5], entity_id="1")

    assert projected.chain_to_seqres["A"] == SEQRES[2:5]


def test_projection_reindexes_densely_from_zero():
    """get_atom_coords walks range(num_res); sparse keys would KeyError."""

    obj = make_object()

    projected = project_mmcif_object(obj, "A", [3, 4, 5], entity_id="1")

    assert sorted(projected.seqres_to_structure["A"]) == [0, 1, 2]


def test_projection_preserves_which_residue_each_index_points_at():
    """The check that catches an off-by-one in the numbering origin."""

    obj = make_object()

    projected = project_mmcif_object(obj, "A", [3, 4, 5], entity_id="1")

    # Gold residue 3 is seqres index 2, whose marker is "res2".
    assert projected.seqres_to_structure["A"][0] == "res2"
    assert projected.seqres_to_structure["A"][1] == "res3"
    assert projected.seqres_to_structure["A"][2] == "res4"


def test_a_nonstandard_numbering_origin_shifts_correctly():
    obj = make_object(start_num=101)

    projected = project_mmcif_object(obj, "A", [103, 104], entity_id="1")

    assert projected.chain_to_seqres["A"] == SEQRES[2:4]
    assert projected.seqres_to_structure["A"][0] == "res2"


def test_projection_drops_every_other_chain():
    """A chain Gold did not retain must not be featurisable by accident."""

    obj = make_object()

    projected = project_mmcif_object(obj, "A", [1, 2], entity_id="1")

    assert set(projected.chain_to_seqres) == {"A"}
    assert set(projected.seqres_to_structure) == {"A"}


def test_coordinates_and_header_are_untouched():
    obj = make_object()

    projected = project_mmcif_object(obj, "A", [1, 2], entity_id="1")

    assert projected.structure is obj.structure
    assert projected.header is obj.header
    assert projected.file_id == obj.file_id


def test_a_full_lineage_reproduces_the_deposited_chain():
    obj = make_object()

    projected = project_mmcif_object(
        obj, "A", list(range(1, len(SEQRES) + 1)), entity_id="1"
    )

    assert projected.chain_to_seqres["A"] == SEQRES


# --------------------------------------------------------------------------
# refusal rather than silent damage
# --------------------------------------------------------------------------

def test_a_lineage_past_the_end_is_refused():
    obj = make_object()

    with pytest.raises(ProjectionError, match="outside the deposited"):
        project_mmcif_object(obj, "A", [20, 21, 22, 23], entity_id="1")


def test_a_lineage_before_the_start_is_refused():
    obj = make_object()

    with pytest.raises(ProjectionError, match="outside the deposited"):
        project_mmcif_object(obj, "A", [0], entity_id="1")


def test_an_absent_auth_chain_is_refused():
    obj = make_object()

    with pytest.raises(ProjectionError, match="absent"):
        project_mmcif_object(obj, "Z", [1, 2], entity_id="1")


def test_a_residue_with_no_structure_mapping_is_refused():
    obj = make_object()

    broken = StubMmcifObject(
        file_id=obj.file_id,
        header=obj.header,
        structure=obj.structure,
        chain_to_seqres=obj.chain_to_seqres,
        seqres_to_structure={"A": {0: "res0"}},
        raw_string=obj.raw_string,
    )

    with pytest.raises(ProjectionError, match="no structure mapping"):
        project_mmcif_object(broken, "A", [1, 2, 3], entity_id="1")


def test_projection_requires_a_numbering_origin():
    obj = make_object()

    with pytest.raises(ProjectionError, match="entity_id or seq_start_num"):
        project_mmcif_object(obj, "A", [1, 2])
