"""Stage 7 lossless same-length Brain prefilter."""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations

import numpy as np
from scipy.spatial import cKDTree


BRAIN_PREFILTER_TAU_ANGSTROM = 0.010
BRAIN_PREFILTER_TAU_MA = 10
BRAIN_DIMENSION = 9


class BrainPrefilterError(ValueError):
    """Raised when Stage-7 Brain candidate search is invalid."""


@dataclass(frozen=True)
class BrainCandidatePairs:
    """Candidate pair indices for one exact-m bucket."""

    m: int
    pairs: np.ndarray

    @property
    def pair_count(self) -> int:
        return int(self.pairs.shape[0])


def _validated_brain_matrix(
    brains: np.ndarray,
) -> np.ndarray:
    value = np.asarray(
        brains,
        dtype=np.float64,
    )

    if value.ndim != 2 or value.shape[1] != BRAIN_DIMENSION:
        raise BrainPrefilterError(
            "Brain matrix must have shape (n, 9)"
        )

    if not np.isfinite(value).all():
        raise BrainPrefilterError(
            "Brain matrix must be finite"
        )

    return value


def brain_integer_numerators(
    brains: np.ndarray,
    *,
    m: int,
) -> np.ndarray:
    """Recover exact sums of canonical BRI milliangstrom coordinates."""

    if (
        not isinstance(m, int)
        or isinstance(m, bool)
        or m < 2
    ):
        raise BrainPrefilterError(
            "Brain-defined search requires integer m >= 2"
        )

    value = _validated_brain_matrix(
        brains
    )

    n = m - 1

    scaled = (
        value
        * 1000.0
        * float(n)
    )

    numerators = np.rint(
        scaled
    ).astype(np.int64)

    # Stage-5 Brain values came from means of canonical 3dp
    # BRI coordinates, so these must reconstruct their integer
    # milliangstrom column sums to floating-point precision.
    reconstructed = (
        numerators.astype(np.float64)
        / (1000.0 * float(n))
    )

    if not np.allclose(
        value,
        reconstructed,
        rtol=0.0,
        atol=1e-12,
    ):
        raise BrainPrefilterError(
            "Brain vector is inconsistent with canonical "
            "3dp BRI integer-sum representation"
        )

    return numerators


def brain_candidate_pairs(
    brains: np.ndarray,
    *,
    m: int,
) -> BrainCandidatePairs:
    """Return all and only same-m pairs satisfying d_Brain <= 0.010 A."""

    numerators = brain_integer_numerators(
        brains,
        m=m,
    )

    chain_count = numerators.shape[0]

    if chain_count < 2:
        return BrainCandidatePairs(
            m=m,
            pairs=np.empty(
                (0, 2),
                dtype=np.int64,
            ),
        )

    # In one exact-m bucket the denominator is common.
    # Thus the scientific threshold tau=0.010 A is exactly
    # 10*(m-1) milliangstrom-sum units.
    exact_radius = (
        BRAIN_PREFILTER_TAU_MA
        * (m - 1)
    )

    tree = cKDTree(
        numerators.astype(
            np.float64,
            copy=False,
        )
    )

    pairs = tree.query_pairs(
        r=float(exact_radius),
        p=np.inf,
        eps=0.0,
        output_type="ndarray",
    )

    pairs = np.asarray(
        pairs,
        dtype=np.int64,
    ).reshape(-1, 2)

    if pairs.shape[0]:
        pairs = pairs[
            np.lexsort(
                (
                    pairs[:, 1],
                    pairs[:, 0],
                )
            )
        ]

        # Exact integer post-filter. This is the scientific
        # d_Brain <= 0.010 decision, not an expanded radius.
        differences = np.max(
            np.abs(
                numerators[pairs[:, 0]]
                - numerators[pairs[:, 1]]
            ),
            axis=1,
        )

        pairs = pairs[
            differences <= exact_radius
        ]

    return BrainCandidatePairs(
        m=m,
        pairs=pairs,
    )


def brute_force_brain_candidate_pairs(
    brains: np.ndarray,
    *,
    m: int,
) -> BrainCandidatePairs:
    """Independent O(n^2) correctness oracle for small buckets."""

    numerators = brain_integer_numerators(
        brains,
        m=m,
    )

    radius = (
        BRAIN_PREFILTER_TAU_MA
        * (m - 1)
    )

    accepted = []

    for i, j in combinations(
        range(numerators.shape[0]),
        2,
    ):
        distance = int(
            np.max(
                np.abs(
                    numerators[i]
                    - numerators[j]
                )
            )
        )

        if distance <= radius:
            accepted.append(
                (i, j)
            )

    return BrainCandidatePairs(
        m=m,
        pairs=np.asarray(
            accepted,
            dtype=np.int64,
        ).reshape(-1, 2),
    )
