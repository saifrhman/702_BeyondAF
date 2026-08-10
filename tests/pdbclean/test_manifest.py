"""Tests for immutable PDB snapshot manifest validation."""

import pyarrow as pa
import pytest

from pdbclean.manifest import (
    ManifestError,
    manifest_partition_count,
    select_manifest_partition,
    validate_manifest_table,
)


def make_table(
    rows: list[dict[str, object]],
) -> pa.Table:
    return pa.Table.from_pylist(rows)


VALID_ROWS = [
    {
        "snapshot": "20260101",
        "source_layout": "canonical_divided_mmcif",
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
        "source_layout": "canonical_divided_mmcif",
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


def _partition_test_table(row_count: int) -> pa.Table:
    return pa.Table.from_pylist(
        [{"row_id": index} for index in range(row_count)]
    )


@pytest.mark.parametrize(
    ("row_count", "batch_size", "expected"),
    [
        (1, 500, 1),
        (499, 500, 1),
        (500, 500, 1),
        (501, 500, 2),
        (1000, 500, 2),
        (1001, 500, 3),
        (17, 4, 5),
        (0, 500, 0),
    ],
)
def test_manifest_partition_count_is_derived_from_runtime_size(
    row_count: int,
    batch_size: int,
    expected: int,
) -> None:
    assert manifest_partition_count(row_count, batch_size) == expected


def test_select_manifest_partition_first_middle_and_final() -> None:
    table = _partition_test_table(12)

    first = select_manifest_partition(
        table,
        task_id=0,
        batch_size=5,
    )
    middle = select_manifest_partition(
        table,
        task_id=1,
        batch_size=5,
    )
    final = select_manifest_partition(
        table,
        task_id=2,
        batch_size=5,
    )

    assert first["row_id"].to_pylist() == [0, 1, 2, 3, 4]
    assert middle["row_id"].to_pylist() == [5, 6, 7, 8, 9]
    assert final["row_id"].to_pylist() == [10, 11]


def test_manifest_partitions_cover_every_row_once() -> None:
    table = _partition_test_table(23)
    batch_size = 7

    partition_count = manifest_partition_count(
        table.num_rows,
        batch_size,
    )

    observed = []

    for task_id in range(partition_count):
        partition = select_manifest_partition(
            table,
            task_id=task_id,
            batch_size=batch_size,
        )
        observed.extend(partition["row_id"].to_pylist())

    assert observed == list(range(23))


def test_select_manifest_partition_rejects_out_of_range_task() -> None:
    table = _partition_test_table(12)

    with pytest.raises(
        ManifestError,
        match=r"task_id 3 is out of range",
    ):
        select_manifest_partition(
            table,
            task_id=3,
            batch_size=5,
        )


def test_select_manifest_partition_rejects_empty_manifest() -> None:
    table = _partition_test_table(0)

    with pytest.raises(
        ManifestError,
        match="empty manifest",
    ):
        select_manifest_partition(
            table,
            task_id=0,
            batch_size=5,
        )


@pytest.mark.parametrize("task_id", [-1, True, 1.5, "1"])
def test_select_manifest_partition_rejects_invalid_task_id(
    task_id,
) -> None:
    table = _partition_test_table(10)

    with pytest.raises(
        ManifestError,
        match="task_id must be a non-negative integer",
    ):
        select_manifest_partition(
            table,
            task_id=task_id,
            batch_size=5,
        )


@pytest.mark.parametrize("batch_size", [0, -1, True, 1.5, "5"])
def test_manifest_partition_rejects_invalid_batch_size(
    batch_size,
) -> None:
    with pytest.raises(
        ManifestError,
        match="batch_size must be a positive integer",
    ):
        manifest_partition_count(
            10,
            batch_size,
        )
