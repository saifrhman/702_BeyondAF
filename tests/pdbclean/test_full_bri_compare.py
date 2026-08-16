import numpy as np
import pytest


def test_identical_bri_has_zero_distance():
    from pdbclean.full_bri_compare import (
        full_bri_distance,
    )

    bri = np.array(
        [
            [1.001, 0.0, 0.0, 0.0, 0.0, 0.0, 2.002, 3.003, 0.0],
            [1.100, 2.200, 3.300, 4.400, 5.500, 6.600, 7.700, 8.800, 9.900],
        ],
        dtype=np.float64,
    )

    result = full_bri_distance(
        bri,
        bri.copy(),
    )

    assert result.retained_residue_count == 2
    assert result.d_bri_mA == 0
    assert result.d_bri == 0.0


def test_exact_0010_angstrom_distance_is_10_mA():
    from pdbclean.full_bri_compare import (
        full_bri_distance,
    )

    query = np.zeros(
        (3, 9),
        dtype=np.float64,
    )
    subject = query.copy()

    subject[2, 8] = 0.010

    result = full_bri_distance(
        query,
        subject,
    )

    assert result.d_bri_mA == 10
    assert result.d_bri == 0.010


def test_0011_angstrom_distance_is_11_mA():
    from pdbclean.full_bri_compare import (
        full_bri_distance,
    )

    query = np.zeros(
        (2, 9),
        dtype=np.float64,
    )
    subject = query.copy()

    subject[1, 4] = -0.011

    result = full_bri_distance(
        query,
        subject,
    )

    assert result.d_bri_mA == 11
    assert result.d_bri == 0.011


def test_uses_max_over_all_rows_and_coordinates():
    from pdbclean.full_bri_compare import (
        full_bri_distance,
    )

    query = np.zeros(
        (4, 9),
        dtype=np.float64,
    )
    subject = query.copy()

    subject[0, 0] = 0.003
    subject[1, 2] = -0.007
    subject[3, 8] = 0.013

    result = full_bri_distance(
        query,
        subject,
    )

    assert result.d_bri_mA == 13


def test_m1_is_valid_for_direct_full_bri_comparison():
    from pdbclean.full_bri_compare import (
        full_bri_distance,
    )

    query = np.zeros(
        (1, 9),
        dtype=np.float64,
    )
    subject = query.copy()

    subject[0, 6] = 0.004

    result = full_bri_distance(
        query,
        subject,
    )

    assert result.retained_residue_count == 1
    assert result.d_bri_mA == 4
    assert result.d_bri == 0.004


def test_different_lengths_are_rejected():
    from pdbclean.full_bri_compare import (
        FullBRICompareError,
        full_bri_distance,
    )

    with pytest.raises(
        FullBRICompareError,
        match="identical BRI shape",
    ):
        full_bri_distance(
            np.zeros((2, 9)),
            np.zeros((3, 9)),
        )


def test_noncanonical_input_is_rejected():
    from pdbclean.full_bri_compare import (
        FullBRICompareError,
        full_bri_distance,
    )

    query = np.zeros(
        (2, 9),
        dtype=np.float64,
    )
    query[1, 0] = 0.0005

    with pytest.raises(
        FullBRICompareError,
        match="canonical 3dp",
    ):
        full_bri_distance(
            query,
            np.zeros((2, 9)),
        )


def test_integer_conversion_is_exact():
    from pdbclean.full_bri_compare import (
        bri_to_integer_mA,
    )

    bri = np.array(
        [
            [-134.481, 119.415, 0.001, -0.001, 0.0, 1.234, -5.678, 9.999, 10.001],
        ],
        dtype=np.float64,
    )

    observed = bri_to_integer_mA(
        bri
    )

    expected = np.array(
        [
            [-134481, 119415, 1, -1, 0, 1234, -5678, 9999, 10001],
        ],
        dtype=np.int64,
    )

    assert np.array_equal(
        observed,
        expected,
    )
