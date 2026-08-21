"""Stage 10 classification from exact full-BRI distance."""

from __future__ import annotations

from dataclasses import dataclass


#: Validated default near-duplicate threshold in exact integer milliangstroms.
#: This is the value that produced the frozen COMP702 20260101 release and it
#: remains the default for every call site.  A caller may pass ``threshold_mA``
#: to use a configured value instead; with no argument the behaviour is
#: identical to the frozen pipeline.  The comparison is inclusive (``<=``).
PAPER_NEAR_DUPLICATE_THRESHOLD_MA = 10


class DuplicateClassificationError(ValueError):
    """Raised when an invalid exact BRI distance is classified."""


def _validated_threshold_mA(threshold_mA: int | None) -> int:
    if threshold_mA is None:
        return PAPER_NEAR_DUPLICATE_THRESHOLD_MA

    if isinstance(threshold_mA, bool) or not isinstance(threshold_mA, int):
        raise DuplicateClassificationError(
            "threshold_mA must be an integer number of milliangstroms"
        )

    if threshold_mA < 0:
        raise DuplicateClassificationError(
            "threshold_mA must be non-negative"
        )

    return threshold_mA


@dataclass(frozen=True)
class DuplicateClassification:
    d_bri_mA: int
    is_zero_duplicate: bool
    is_paper_near_duplicate: bool
    is_nonzero_near_duplicate: bool


def classify_bri_distance(
    d_bri_mA: int,
    *,
    threshold_mA: int | None = None,
) -> DuplicateClassification:
    """Classify one exact integer-mA complete-BRI distance.

    ``threshold_mA`` defaults to
    :data:`PAPER_NEAR_DUPLICATE_THRESHOLD_MA` (10 mA = 0.010 A).
    """

    threshold = _validated_threshold_mA(threshold_mA)

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
        <= threshold
    )

    nonzero_near = (
        0
        < d_bri_mA
        <= threshold
    )

    return DuplicateClassification(
        d_bri_mA=d_bri_mA,
        is_zero_duplicate=zero,
        is_paper_near_duplicate=paper_near,
        is_nonzero_near_duplicate=nonzero_near,
    )
