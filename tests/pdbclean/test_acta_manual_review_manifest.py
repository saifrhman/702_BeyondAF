import pytest

from pdbclean.acta_manual_review_manifest import (
    ActaManualReviewManifestError,
    MANIFEST_SCHEMA,
    _canonical_deposition_pair,
)


def test_canonical_pair_preserves_forward_orientation():
    result = _canonical_deposition_pair(
        "1abc",
        "A",
        "2xyz",
        "B",
    )

    assert result == (
        "1abc",
        "2xyz",
        "A",
        "B",
    )


def test_canonical_pair_swaps_reverse_orientation():
    result = _canonical_deposition_pair(
        "9xyz",
        "Q",
        "1abc",
        "S",
    )

    assert result == (
        "1abc",
        "9xyz",
        "S",
        "Q",
    )


def test_same_deposition_is_forbidden():
    with pytest.raises(
        ActaManualReviewManifestError
    ):
        _canonical_deposition_pair(
            "1abc",
            "A",
            "1ABC",
            "B",
        )


def test_manifest_explicitly_records_no_filtering():
    metadata = MANIFEST_SCHEMA.metadata

    assert (
        metadata[
            b"scientific_filtering"
        ]
        == b"none"
    )

    assert (
        metadata[
            b"grouping_unit"
        ]
        == b"unordered_PDB_deposition_pair"
    )



def test_stage11_v1_publication_layout_is_preserved():
    from pdbclean.acta_manual_review_manifest import (
        _publication_layout,
    )

    assert _publication_layout("1.0") == (
        "acta_downstream_investigation",
        "acta_manual_review_manifest",
    )


def test_stage11_v2_publication_layout_is_isolated():
    from pdbclean.acta_manual_review_manifest import (
        _publication_layout,
    )

    assert _publication_layout("2.0") == (
        "acta_downstream_investigation_v2",
        "acta_manual_review_manifest_v2",
    )
