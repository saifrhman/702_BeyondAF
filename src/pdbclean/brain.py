"""Backbone Rigid Average Invariant (Brain).

This module implements the nine-coordinate average invariant from
MATCH Definition 5.1 for an already-canonical Stage-3 BRI matrix.

No BRI computation, cleaning, duplicate classification, or candidate
search occurs here.
"""

from __future__ import annotations

import numpy as np


BRAIN_DIMENSION = 9


def compute_brain(
    bri: np.ndarray,
) -> np.ndarray:
    """Compute Definition 5.1 Brain from one canonical BRI matrix.

    Brain is the vector of nine column means of BRI excluding its
    first row. Therefore it is defined only for backbones containing
    at least two residues.

    The canonical Stage-3 BRI values are used directly. The resulting
    mathematical means are deliberately not rounded.
    """

    value = np.asarray(
        bri,
        dtype=np.float64,
    )

    if value.ndim != 2:
        raise ValueError(
            "Brain requires a two-dimensional BRI matrix"
        )

    if value.shape[1] != BRAIN_DIMENSION:
        raise ValueError(
            "Brain requires a BRI matrix with exactly 9 columns"
        )

    if value.shape[0] < 2:
        raise ValueError(
            "Definition 5.1 Brain is undefined for m < 2"
        )

    if not np.isfinite(value).all():
        raise ValueError(
            "Brain requires finite BRI coordinates"
        )

    if not np.array_equal(
        value,
        np.around(value, 3),
    ):
        raise ValueError(
            "Brain requires canonical 3-decimal Stage-3 BRI"
        )

    brain = np.mean(
        value[1:, :],
        axis=0,
        dtype=np.float64,
    )

    if (
        brain.shape != (BRAIN_DIMENSION,)
        or brain.dtype != np.float64
        or not np.isfinite(brain).all()
    ):
        raise ValueError(
            "Definition 5.1 Brain computation produced an "
            "invalid result"
        )

    return brain
