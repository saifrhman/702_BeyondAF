"""Tests for immutable PDB snapshot manifest validation."""

import pyarrow as pa
import pytest

from pdbclean.manifest import ManifestError, validate_manifest_table


def make_table(
    rows: list[dict[str, object]],
) -> pa.Table:
    return pa.Table.from_pylist(rows)


VALID_ROWS = [
    {
        "snapshot": "20260101",
        "pdb_id": "100d",
        "s3_key": (
            "20260101/pub/pdb/data/structures/"
            "divided/mmCIF/00/100d.cif.gz"
        ),
        "size_bytes": 23548,
        "etag": "etag-100d",
    },
    {
        "snapshot": "20260101",
        "pdb_id": "200l",
        "s3_key": (
            "20260101/pub/pdb/data/structures/"
            "divided/mmCIF/00/200l.cif.gz"
        ),
        "size_bytes": 48190,
        "etag": "etag-200l",
    },
]


def test_valid_manifest() -> None:
    summary = validate_manifest_table(
        make_table(VALID_ROWS),
        expected_snapshot="20260101",
        expected_count=2,
        expected_total_bytes=71738,
    )

    assert summary.row_count == 2
    assert summary.total_bytes == 71738
    assert summary.unique_pdb_ids == 2
    assert summary.unique_s3_keys == 2


def test_duplicate_pdb_id_is_rejected() -> None:
    rows = [
        VALID_ROWS[0],
        {
            **VALID_ROWS[0],
            "s3_key": (
                "20260101/pub/pdb/data/structures/"
                "divided/mmCIF/00/copy.cif.gz"
            ),
        },
    ]

    with pytest.raises(
        ManifestError,
        match="duplicate PDB ID",
    ):
        validate_manifest_table(
            make_table(rows),
            expected_snapshot="20260101",
        )


def test_incorrect_snapshot_prefix_is_rejected() -> None:
    rows = [
        {
            **VALID_ROWS[0],
            "s3_key": (
                "20250101/pub/pdb/data/structures/"
                "divided/mmCIF/00/100d.cif.gz"
            ),
        }
    ]

    with pytest.raises(
        ManifestError,
        match="outside snapshot prefix",
    ):
        validate_manifest_table(
            make_table(rows),
            expected_snapshot="20260101",
        )


def test_incorrect_divided_folder_is_rejected() -> None:
    rows = [
        {
            **VALID_ROWS[0],
            "s3_key": (
                "20260101/pub/pdb/data/structures/"
                "divided/mmCIF/10/100d.cif.gz"
            ),
        }
    ]

    with pytest.raises(
        ManifestError,
        match="does not match PDB layout",
    ):
        validate_manifest_table(
            make_table(rows),
            expected_snapshot="20260101",
        )


def test_expected_count_is_enforced() -> None:
    with pytest.raises(
        ManifestError,
        match="Expected 3 rows, found 2",
    ):
        validate_manifest_table(
            make_table(VALID_ROWS),
            expected_snapshot="20260101",
            expected_count=3,
        )


def test_expected_total_bytes_is_enforced() -> None:
    with pytest.raises(
        ManifestError,
        match="Expected total size 1, found 71738",
    ):
        validate_manifest_table(
            make_table(VALID_ROWS),
            expected_snapshot="20260101",
            expected_total_bytes=1,
        )
