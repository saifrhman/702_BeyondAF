import pytest


@pytest.mark.parametrize(
    "distance,zero,paper_near,nonzero_near",
    [
        (0, True, True, False),
        (1, False, True, True),
        (9, False, True, True),
        (10, False, True, True),
        (11, False, False, False),
        (1000, False, False, False),
    ],
)
def test_exact_stage10_classification(
    distance,
    zero,
    paper_near,
    nonzero_near,
):
    from pdbclean.duplicate_classification import (
        classify_bri_distance,
    )

    result = classify_bri_distance(
        distance
    )

    assert result.is_zero_duplicate is zero
    assert result.is_paper_near_duplicate is paper_near
    assert result.is_nonzero_near_duplicate is nonzero_near


@pytest.mark.parametrize(
    "value",
    [-1, 0.0, True, "10"],
)
def test_invalid_distance_is_rejected(
    value,
):
    from pdbclean.duplicate_classification import (
        DuplicateClassificationError,
        classify_bri_distance,
    )

    with pytest.raises(
        DuplicateClassificationError,
    ):
        classify_bri_distance(
            value
        )
