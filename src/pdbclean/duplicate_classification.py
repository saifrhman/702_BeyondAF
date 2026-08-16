"""Stage 10 classification from exact full-BRI distance."""

from __future__ import annotations

from dataclasses import dataclass


PAPER_NEAR_DUPLICATE_THRESHOLD_MA = 10


class DuplicateClassificationError(ValueError):
    """Raised when an invalid exact BRI distance is classified."""


@dataclass(frozen=True)
class DuplicateClassification:
    d_bri_mA: int
    is_zero_duplicate: bool
    is_paper_near_duplicate: bool
    is_nonzero_near_duplicate: bool


def classify_bri_distance(
    d_bri_mA: int,
) -> DuplicateClassification:
    """Classify one exact integer-mA full-BRI distance."""

    if (
        not isinstance(d_bri_mA, int)
        or isinstance(d_bri_mA, bool)
        or d_bri_mA < 0
    ):
        raise DuplicateClassificationError(
            "d_bri_mA must be a nonnegative integer"
        )

    zero = d_bri_mA == 0

    paper_near = (
        d_bri_mA
        <= PAPER_NEAR_DUPLICATE_THRESHOLD_MA
    )

    nonzero_near = (
        0
        < d_bri_mA
        <= PAPER_NEAR_DUPLICATE_THRESHOLD_MA
    )

    return DuplicateClassification(
        d_bri_mA=d_bri_mA,
        is_zero_duplicate=zero,
        is_paper_near_duplicate=paper_near,
        is_nonzero_near_duplicate=nonzero_near,
    )
