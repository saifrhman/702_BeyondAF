"""Tests for MATCH Definition 5.1 Brain."""

from __future__ import annotations

import numpy as np
import pytest

from pdbclean.brain import (
    BRAIN_DIMENSION,
    compute_brain,
)


def test_brain_is_mean_of_rows_after_first() -> None:
    bri = np.asarray(
        [
            [99.0] * 9,
            [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0],
            [3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0],
        ],
        dtype=np.float64,
    )

    observed = compute_brain(bri)

    expected = np.asarray(
        [2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0],
        dtype=np.float64,
    )

    assert observed.shape == (BRAIN_DIMENSION,)
    assert observed.dtype == np.float64
    assert np.array_equal(observed, expected)


def test_brain_excludes_first_row() -> None:
    tail = np.asarray(
        [
            [1.0] * 9,
            [3.0] * 9,
        ],
        dtype=np.float64,
    )

    first_a = np.asarray(
        [[0.0] * 9],
        dtype=np.float64,
    )
    first_b = np.asarray(
        [[1000.0] * 9],
        dtype=np.float64,
    )

    assert np.array_equal(
        compute_brain(
            np.vstack((first_a, tail))
        ),
        compute_brain(
            np.vstack((first_b, tail))
        ),
    )


def test_brain_is_not_rounded_after_mean() -> None:
    bri = np.asarray(
        [
            [0.0] * 9,
            [0.001] * 9,
            [0.002] * 9,
            [0.002] * 9,
        ],
        dtype=np.float64,
    )

    observed = compute_brain(bri)

    expected = np.full(
        9,
        0.005 / 3.0,
        dtype=np.float64,
    )

    assert np.array_equal(
        observed,
        expected,
    )

    assert not np.array_equal(
        observed,
        np.around(observed, 3),
    )


def test_brain_for_two_residues_is_second_bri_row() -> None:
    second = np.asarray(
        [
            0.101,
            -0.202,
            0.303,
            1.404,
            -1.505,
            2.606,
            -2.707,
            3.808,
            -3.909,
        ],
        dtype=np.float64,
    )

    bri = np.vstack(
        (
            np.zeros(9, dtype=np.float64),
            second,
        )
    )

    assert np.array_equal(
        compute_brain(bri),
        second,
    )


def test_brain_rejects_m1_instead_of_returning_nan() -> None:
    with pytest.raises(
        ValueError,
        match="undefined for m < 2",
    ):
        compute_brain(
            np.zeros(
                (1, 9),
                dtype=np.float64,
            )
        )


@pytest.mark.parametrize(
    "shape",
    [
        (2, 8),
        (2, 10),
    ],
)
def test_brain_rejects_wrong_coordinate_dimension(
    shape: tuple[int, int],
) -> None:
    with pytest.raises(
        ValueError,
        match="exactly 9 columns",
    ):
        compute_brain(
            np.zeros(
                shape,
                dtype=np.float64,
            )
        )


def test_brain_rejects_noncanonical_bri_precision() -> None:
    bri = np.zeros(
        (2, 9),
        dtype=np.float64,
    )
    bri[1, 0] = 0.0005

    with pytest.raises(
        ValueError,
        match="canonical 3-decimal",
    ):
        compute_brain(bri)


def test_brain_rejects_nonfinite_bri() -> None:
    bri = np.zeros(
        (2, 9),
        dtype=np.float64,
    )
    bri[1, 0] = np.nan

    with pytest.raises(
        ValueError,
        match="finite",
    ):
        compute_brain(bri)


def test_brain_satisfies_lemma_6_1_on_example() -> None:
    left = np.asarray(
        [
            [0.0] * 9,
            [1.000] * 9,
            [2.000] * 9,
            [3.000] * 9,
        ],
        dtype=np.float64,
    )

    right = np.asarray(
        [
            [0.0] * 9,
            [1.005] * 9,
            [1.993] * 9,
            [3.010] * 9,
        ],
        dtype=np.float64,
    )

    full_distance = float(
        np.max(
            np.abs(left - right)
        )
    )

    brain_distance = float(
        np.max(
            np.abs(
                compute_brain(left)
                - compute_brain(right)
            )
        )
    )

    assert brain_distance <= full_distance
