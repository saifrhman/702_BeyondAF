"""Stage 8 exact full-BRI L-infinity comparison."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


BRI_COORDINATE_DIMENSION = 9
BRI_MA_PER_ANGSTROM = 1000


class FullBRICompareError(ValueError):
    """Raised when an exact full-BRI comparison is invalid."""


@dataclass(frozen=True)
class FullBRIDistance:
    """Exact L-infinity distance between two same-length BRIs."""

    retained_residue_count: int
    d_bri_mA: int

    @property
    def d_bri(self) -> float:
        return self.d_bri_mA / BRI_MA_PER_ANGSTROM


def _validated_canonical_bri(
    bri: np.ndarray,
    *,
    name: str,
) -> np.ndarray:
    value = np.asarray(
        bri,
        dtype=np.float64,
    )

    if (
        value.ndim != 2
        or value.shape[1] != BRI_COORDINATE_DIMENSION
        or value.shape[0] < 1
    ):
        raise FullBRICompareError(
            f"{name} BRI must have shape (m, 9) with m >= 1"
        )

    if not np.isfinite(value).all():
        raise FullBRICompareError(
            f"{name} BRI must be finite"
        )

    if not np.array_equal(
        value,
        np.around(value, 3),
    ):
        raise FullBRICompareError(
            f"{name} BRI must be canonical 3dp"
        )

    return value


def bri_to_integer_mA(
    bri: np.ndarray,
    *,
    name: str = "input",
) -> np.ndarray:
    """Convert canonical 3dp BRI coordinates losslessly to integer mA."""

    value = _validated_canonical_bri(
        bri,
        name=name,
    )

    scaled = (
        value
        * BRI_MA_PER_ANGSTROM
    )

    integer = np.rint(
        scaled
    ).astype(
        np.int64
    )

    reconstructed = (
        integer.astype(np.float64)
        / BRI_MA_PER_ANGSTROM
    )

    if not np.array_equal(
        reconstructed,
        value,
    ):
        raise FullBRICompareError(
            f"{name} BRI is not losslessly representable in integer mA"
        )

    return integer


def full_bri_distance(
    query_bri: np.ndarray,
    subject_bri: np.ndarray,
) -> FullBRIDistance:
    """Compute exact same-length full-BRI L-infinity distance."""

    query = bri_to_integer_mA(
        query_bri,
        name="query",
    )

    subject = bri_to_integer_mA(
        subject_bri,
        name="subject",
    )

    if query.shape != subject.shape:
        raise FullBRICompareError(
            "Full-BRI comparison requires identical BRI shape / m"
        )

    difference = np.abs(
        query - subject
    )

    d_bri_mA = int(
        np.max(difference)
    )

    return FullBRIDistance(
        retained_residue_count=query.shape[0],
        d_bri_mA=d_bri_mA,
    )
