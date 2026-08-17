import pytest

from pdbclean.downstream_metadata_task import (
    DownstreamMetadataTaskError,
    _task_count,
)


@pytest.mark.parametrize(
    ("rows", "batch_size", "expected"),
    [
        (1, 500, 1),
        (500, 500, 1),
        (501, 500, 2),
        (7259, 500, 15),
    ],
)
def test_task_count(
    rows,
    batch_size,
    expected,
):
    assert (
        _task_count(
            rows,
            batch_size,
        )
        == expected
    )


def test_task_count_rejects_empty_population():
    with pytest.raises(
        DownstreamMetadataTaskError
    ):
        _task_count(
            0,
            500,
        )


def test_task_count_rejects_invalid_batch_size():
    with pytest.raises(
        DownstreamMetadataTaskError
    ):
        _task_count(
            10,
            0,
        )
