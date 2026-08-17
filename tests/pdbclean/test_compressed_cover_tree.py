import numpy as np
import pytest

from pdbclean.compressed_cover_tree import (
    CompressedCoverTree,
    linf_distance_mA,
)


def _bruteforce(points, query, k):
    ranked = sorted(
        (
            linf_distance_mA(query, point),
            index,
        )
        for index, point in enumerate(points)
    )
    return ranked[:k]


def test_tree_axioms_hold_on_integer_linf_points():
    rng = np.random.default_rng(702)

    points = rng.integers(
        -500,
        501,
        size=(40, 12),
        dtype=np.int64,
    )

    tree = CompressedCoverTree(points)
    tree.validate()


def test_exact_knn_matches_bruteforce():
    rng = np.random.default_rng(703)

    points = rng.integers(
        -1000,
        1001,
        size=(60, 9),
        dtype=np.int64,
    )

    tree = CompressedCoverTree(points)
    tree.validate()

    for _ in range(20):
        query = rng.integers(
            -1000,
            1001,
            size=9,
            dtype=np.int64,
        )

        for k in (1, 2, 5, 13, 60):
            expected = _bruteforce(
                points,
                query,
                k,
            )

            result = tree.knn(
                query,
                k,
            )

            observed = list(
                zip(
                    result.distances_mA.tolist(),
                    result.indices.tolist(),
                )
            )

            assert observed == expected


def test_radius_search_matches_bruteforce():
    rng = np.random.default_rng(704)

    points = rng.integers(
        -30,
        31,
        size=(80, 15),
        dtype=np.int64,
    )

    tree = CompressedCoverTree(points)
    tree.validate()

    for _ in range(15):
        query = rng.integers(
            -30,
            31,
            size=15,
            dtype=np.int64,
        )

        for radius in (0, 5, 10, 20, 40, 80):
            expected = sorted(
                (
                    linf_distance_mA(
                        query,
                        point,
                    ),
                    index,
                )
                for index, point
                in enumerate(points)
                if linf_distance_mA(
                    query,
                    point,
                )
                <= radius
            )

            result = tree.radius_neighbors(
                query,
                radius,
            )

            observed = list(
                zip(
                    result.distances_mA.tolist(),
                    result.indices.tolist(),
                )
            )

            assert observed == expected


def test_query_may_equal_reference_point():
    points = np.asarray(
        [
            [0, 0, 0],
            [1, 0, 0],
            [5, 5, 5],
            [10, 10, 10],
        ],
        dtype=np.int64,
    )

    tree = CompressedCoverTree(points)
    tree.validate()

    result = tree.knn(
        points[2],
        3,
    )

    expected = _bruteforce(
        points,
        points[2],
        3,
    )

    assert list(
        zip(
            result.distances_mA.tolist(),
            result.indices.tolist(),
        )
    ) == expected


def test_duplicate_reference_vectors_must_be_collapsed_first():
    points = np.asarray(
        [
            [1, 2, 3],
            [1, 2, 3],
        ],
        dtype=np.int64,
    )

    with pytest.raises(
        ValueError,
        match="duplicate metric vectors",
    ):
        CompressedCoverTree(points)
