import numpy as np
import pytest


def _brain_from_sums(
    sums,
    *,
    m,
):
    return (
        np.asarray(
            sums,
            dtype=np.float64,
        )
        / (1000.0 * (m - 1))
    )


def test_exact_zero_distance_pair_is_retained():
    from pdbclean.brain_prefilter import (
        brain_candidate_pairs,
    )

    brains = np.vstack(
        [
            _brain_from_sums([100] * 9, m=3),
            _brain_from_sums([100] * 9, m=3),
        ]
    )

    result = brain_candidate_pairs(
        brains,
        m=3,
    )

    assert result.pairs.tolist() == [[0, 1]]


def test_exact_0010_angstrom_boundary_is_retained():
    from pdbclean.brain_prefilter import (
        brain_candidate_pairs,
    )

    # m=3 => n=2 => exact radius = 20 sum-mA units.
    brains = np.vstack(
        [
            _brain_from_sums([0] * 9, m=3),
            _brain_from_sums([20] + [0] * 8, m=3),
        ]
    )

    result = brain_candidate_pairs(
        brains,
        m=3,
    )

    assert result.pairs.tolist() == [[0, 1]]


def test_just_outside_boundary_is_rejected():
    from pdbclean.brain_prefilter import (
        brain_candidate_pairs,
    )

    brains = np.vstack(
        [
            _brain_from_sums([0] * 9, m=3),
            _brain_from_sums([21] + [0] * 8, m=3),
        ]
    )

    result = brain_candidate_pairs(
        brains,
        m=3,
    )

    assert result.pair_count == 0


def test_chebyshev_uses_maximum_of_nine_coordinates():
    from pdbclean.brain_prefilter import (
        brain_candidate_pairs,
    )

    brains = np.vstack(
        [
            _brain_from_sums([0] * 9, m=6),
            _brain_from_sums(
                [50, -50, 10, 0, 1, 2, 3, 4, 5],
                m=6,
            ),
        ]
    )

    result = brain_candidate_pairs(
        brains,
        m=6,
    )

    assert result.pairs.tolist() == [[0, 1]]


def test_m1_is_not_sent_through_brain_search():
    from pdbclean.brain_prefilter import (
        BrainPrefilterError,
        brain_candidate_pairs,
    )

    with pytest.raises(
        BrainPrefilterError,
        match="m >= 2",
    ):
        brain_candidate_pairs(
            np.empty((0, 9)),
            m=1,
        )


def test_single_chain_bucket_has_no_pairs():
    from pdbclean.brain_prefilter import (
        brain_candidate_pairs,
    )

    result = brain_candidate_pairs(
        _brain_from_sums(
            [[10] * 9],
            m=2,
        ),
        m=2,
    )

    assert result.pairs.shape == (0, 2)


@pytest.mark.parametrize(
    "m,seed,n",
    [
        (2, 1, 25),
        (10, 2, 40),
        (99, 3, 60),
        (231, 4, 50),
    ],
)
def test_ckdtree_exactly_matches_brute_force(
    m,
    seed,
    n,
):
    from pdbclean.brain_prefilter import (
        brain_candidate_pairs,
        brute_force_brain_candidate_pairs,
    )

    rng = np.random.default_rng(
        seed
    )

    sums = rng.integers(
        -2000,
        2001,
        size=(n, 9),
        dtype=np.int64,
    )

    # Inject guaranteed near/boundary pairs.
    if n >= 4:
        sums[1] = sums[0]
        sums[2] = sums[0]
        sums[2, 0] += 10 * (m - 1)
        sums[3] = sums[0]
        sums[3, 0] += 10 * (m - 1) + 1

    brains = (
        sums.astype(np.float64)
        / (1000.0 * (m - 1))
    )

    optimized = brain_candidate_pairs(
        brains,
        m=m,
    )
    brute = brute_force_brain_candidate_pairs(
        brains,
        m=m,
    )

    assert np.array_equal(
        optimized.pairs,
        brute.pairs,
    )
