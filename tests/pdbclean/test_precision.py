"""BRI representation precision as an explicit experimental control.

Precision p is configuration, not a source edit. These tests pin that:

* the validated default is 0.001 A and reproduces the frozen semantics;
* changing p produces a *distinct scientific configuration*;
* thresholds are validated against the configured grid, never silently rounded;
* an experimental precision can never reuse validated output;
* p and tau are separate axes and are never conflated.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from pdbclean.defaults import (
    IMPLEMENTED_REPRESENTATION_PRECISION_ANGSTROM,
    VALIDATED_DEFAULTS,
    PrecisionNotImplementedError,
    brain_filter_threshold_units,
    near_duplicate_threshold_units,
    precision_is_implemented,
    quantise_angstrom_to_units,
    representation_precision_angstrom,
    representation_unit_label,
    require_implemented_precision,
    validated_defaults,
)
from pdbclean.runconfig import RunConfigError, resolve_run_config


REPO = Path(__file__).resolve().parents[2]
PROFILE = REPO / "config" / "pdbclean" / "profiles" / "comp702_frozen_20260101.yaml"


def _experimental(precision, *, tau=0.010, brain=0.010):
    return resolve_run_config(
        overrides={
            "bri.representation_precision_angstrom": precision,
            "duplicate_search.near_duplicate_threshold_angstrom": tau,
            "brain_filter.threshold_angstrom": brain,
        }
    )


# --------------------------------------------------------------------------
# A. Default precision
# --------------------------------------------------------------------------


def test_validated_default_precision_is_one_milliangstrom():
    assert VALIDATED_DEFAULTS["bri"]["representation_precision_angstrom"] == 0.001
    assert IMPLEMENTED_REPRESENTATION_PRECISION_ANGSTROM == 0.001

    resolved = resolve_run_config()

    assert resolved.get("bri.representation_precision_angstrom") == 0.001
    assert representation_precision_angstrom(resolved.data) == 0.001


def test_default_precision_reproduces_the_frozen_unit_semantics():
    """At p = 0.001, one unit is 1 mA and tau = 0.010 A is 10 units."""

    defaults = validated_defaults()

    assert near_duplicate_threshold_units(defaults) == 10
    assert brain_filter_threshold_units(defaults) == 10
    assert representation_unit_label(0.001) == "mA"


def test_default_precision_is_the_implemented_one():
    assert precision_is_implemented(validated_defaults()) is True
    assert require_implemented_precision(
        validated_defaults(), stage="test"
    ) == 0.001


def test_frozen_profile_keeps_the_validated_precision():
    resolved = resolve_run_config(config_path=PROFILE)

    assert resolved.get("bri.representation_precision_angstrom") == 0.001
    assert resolved.near_duplicate_threshold_mA == 10
    assert resolved.brain_threshold_mA == 10


def test_frozen_scientific_identity_is_unchanged_by_exposing_precision():
    """Making p configurable must not move the frozen run's identity."""

    resolved = resolve_run_config(config_path=PROFILE)

    assert resolved.scientific_sha256 == (
        "25b8e62a87cb90797af41cd4149dfd4280e3a7aed99428e70fa97117c5bababa"
    )


# --------------------------------------------------------------------------
# B. Precision override is a distinct scientific configuration
# --------------------------------------------------------------------------


def test_changing_precision_changes_the_scientific_identity():
    baseline = resolve_run_config()
    experimental = _experimental(0.002)

    assert experimental.scientific_sha256 != baseline.scientific_sha256


def test_precision_is_written_into_the_resolved_configuration():
    experimental = _experimental(0.002)

    assert experimental.get("bri.representation_precision_angstrom") == 0.002
    assert experimental.to_dict()["bri"][
        "representation_precision_angstrom"
    ] == 0.002


def test_precision_override_records_its_origin():
    experimental = _experimental(0.002)

    assert experimental.source_of(
        "bri.representation_precision_angstrom"
    ).startswith("override:")


def test_precision_is_recorded_in_run_provenance(tmp_path):
    from pdbclean.run_provenance import RunProvenance

    experimental = resolve_run_config(
        overrides={
            "snapshot.mode": "fixed",
            "snapshot.snapshot_id": "20260101",
            "bri.representation_precision_angstrom": 0.002,
            "duplicate_search.near_duplicate_threshold_angstrom": 0.010,
            "brain_filter.threshold_angstrom": 0.010,
        }
    )

    run = RunProvenance.create(
        resolved=experimental,
        run_root=tmp_path / "runs",
        repo_root=tmp_path,
    )

    record = run.record

    assert record["resolved_config"]["bri"][
        "representation_precision_angstrom"
    ] == 0.002
    assert record["scientific_config_sha256"] == (
        experimental.scientific_sha256
    )
    assert record["config_value_sources"][
        "bri.representation_precision_angstrom"
    ].startswith("override:")


def test_precision_is_not_an_execution_setting():
    """p must sit in the scientific projection, not with the executor."""

    from pdbclean.runconfig import scientific_projection

    projection = scientific_projection(resolve_run_config().data)

    assert "bri" in projection
    assert projection["bri"]["representation_precision_angstrom"] == 0.001


# --------------------------------------------------------------------------
# C. Grid validation -- never silently round
# --------------------------------------------------------------------------


@pytest.mark.parametrize(
    "precision,threshold,units",
    [
        (0.001, 0.010, 10),
        (0.002, 0.010, 5),
        (0.005, 0.010, 2),
        (0.010, 0.010, 1),
        (0.001, 0.005, 5),
    ],
)
def test_compatible_threshold_and_precision_are_accepted(
    precision, threshold, units
):
    assert quantise_angstrom_to_units(threshold, precision) == units

    resolved = _experimental(precision, tau=threshold, brain=threshold)

    assert resolved.near_duplicate_threshold_mA == units


@pytest.mark.parametrize(
    "precision,threshold",
    [
        (0.003, 0.010),
        (0.004, 0.010),
        (0.007, 0.010),
        (0.002, 0.0055),
    ],
)
def test_incompatible_threshold_and_precision_are_rejected(
    precision, threshold
):
    with pytest.raises(ValueError, match="not exactly representable"):
        quantise_angstrom_to_units(threshold, precision)

    with pytest.raises(RunConfigError):
        _experimental(precision, tau=threshold, brain=threshold)


def test_rejection_is_not_a_silent_rounding():
    """0.003 / 0.010 must fail, never quietly become another threshold."""

    with pytest.raises(RunConfigError) as excinfo:
        _experimental(0.003, tau=0.010, brain=0.010)

    message = str(excinfo.value)

    assert "not exactly representable" in message
    assert "0.003" in message


def test_threshold_below_one_unit_is_rejected():
    with pytest.raises(ValueError):
        quantise_angstrom_to_units(0.0005, 0.001)


def test_non_default_precision_uses_generic_unit_terminology():
    """Only the validated grid may be described in milliangstroms."""

    assert representation_unit_label(0.001) == "mA"
    assert representation_unit_label(0.002) == "unit"
    assert representation_unit_label(0.01) == "unit"


# --------------------------------------------------------------------------
# Execution gate: an unimplemented grid never silently produces output
# --------------------------------------------------------------------------


def test_unimplemented_precision_refuses_to_execute():
    config = validated_defaults()
    config["bri"]["representation_precision_angstrom"] = 0.002

    assert precision_is_implemented(config) is False

    with pytest.raises(PrecisionNotImplementedError, match="0.002"):
        require_implemented_precision(config, stage="stage3_complete_bri")


def test_execution_gate_names_the_stage_and_the_reason():
    config = validated_defaults()
    config["bri"]["representation_precision_angstrom"] = 0.005

    with pytest.raises(PrecisionNotImplementedError) as excinfo:
        require_implemented_precision(config, stage="stage7_brain_prefilter")

    message = str(excinfo.value)

    assert "stage7_brain_prefilter" in message
    assert "1.2.2" in message
    assert "not implemented" in message


# --------------------------------------------------------------------------
# Reuse safety: an experimental precision never reuses validated output
# --------------------------------------------------------------------------


def test_experimental_precision_cannot_reuse_validated_output():
    from pdbclean.pipeline import ACTION_REUSE, plan_pipeline
    from pdbclean.stage_registry import PRECISION_DEPENDENT_STAGES

    experimental = resolve_run_config(
        config_path=PROFILE,
        overrides={
            "bri.representation_precision_angstrom": 0.002,
            "duplicate_search.near_duplicate_threshold_angstrom": 0.010,
            "brain_filter.threshold_angstrom": 0.010,
        },
    )

    plan = plan_pipeline(experimental, repo_root=REPO)

    by_id = plan.by_id

    for stage_id in PRECISION_DEPENDENT_STAGES:
        observation = by_id.get(stage_id)

        if observation is None:
            continue

        assert observation.action != ACTION_REUSE, stage_id


def test_validated_precision_still_reuses_validated_output():
    from pdbclean.pipeline import ACTION_NOT_APPLICABLE, ACTION_REUSE, plan_pipeline

    plan = plan_pipeline(
        resolve_run_config(config_path=PROFILE), repo_root=REPO
    )

    if not any(observation.exists for observation in plan.observations):
        pytest.skip("frozen outputs are not present in this checkout")

    assert all(
        observation.action in {ACTION_REUSE, ACTION_NOT_APPLICABLE}
        for observation in plan.observations
    )


# --------------------------------------------------------------------------
# p and tau are separate axes (section 47)
# --------------------------------------------------------------------------


def test_precision_and_threshold_are_independent_axes():
    baseline = resolve_run_config()

    precision_only = _experimental(0.002)
    threshold_only = resolve_run_config(
        overrides={
            "duplicate_search.near_duplicate_threshold_angstrom": 0.005,
        }
    )

    identities = {
        baseline.scientific_sha256,
        precision_only.scientific_sha256,
        threshold_only.scientific_sha256,
    }

    assert len(identities) == 3


def test_precision_is_not_the_duplicate_threshold():
    resolved = resolve_run_config()

    assert resolved.get("bri.representation_precision_angstrom") == 0.001
    assert resolved.get(
        "duplicate_search.near_duplicate_threshold_angstrom"
    ) == 0.010
    assert resolved.get(
        "bri.representation_precision_angstrom"
    ) != resolved.get("duplicate_search.near_duplicate_threshold_angstrom")


# --------------------------------------------------------------------------
# Brain filter safety is preserved under configurable precision
# --------------------------------------------------------------------------


def test_brain_filter_may_not_be_tighter_than_the_classifier_at_any_p():
    """The lossless-prefilter guarantee holds on every configured grid."""

    for precision, brain, tau in ((0.001, 0.005, 0.010), (0.002, 0.004, 0.010)):
        with pytest.raises(RunConfigError):
            _experimental(precision, tau=tau, brain=brain)


def test_brain_filter_equal_to_the_classifier_is_accepted_at_any_p():
    for precision in (0.001, 0.002, 0.005):
        resolved = _experimental(precision, tau=0.010, brain=0.010)

        assert resolved.brain_threshold_mA == resolved.near_duplicate_threshold_mA


def test_brain_prefilter_remains_lossless_against_the_oracle():
    """Candidate safety, verified against the brute-force oracle."""

    import numpy as np

    from pdbclean.brain_prefilter import (
        brain_candidate_pairs,
        brute_force_brain_candidate_pairs,
    )

    rng = np.random.default_rng(101)

    for m, tau_units in ((3, 10), (8, 10), (8, 5), (25, 20)):
        brains = np.around(rng.random((40, 9)) * 0.03, 3)
        brains[5] = brains[4]
        brains[9] = np.around(brains[8] + 0.0004, 3)

        fast = brain_candidate_pairs(brains, m=m, tau_mA=tau_units)
        oracle = brute_force_brain_candidate_pairs(
            brains, m=m, tau_mA=tau_units
        )

        assert np.array_equal(fast.pairs, oracle.pairs), (m, tau_units)


def test_brain_prefilter_boundary_pairs_are_retained():
    """A pair exactly at the threshold must survive the prefilter."""

    import numpy as np

    from pdbclean.brain_prefilter import brain_candidate_pairs

    for tau_units in (5, 10, 20):
        step = tau_units / 1000.0

        brains = np.array(
            [
                [0.0] * 9,
                [step] * 9,             # exactly at tau
                [step + 0.001] * 9,     # one unit outside
            ]
        )

        pairs = {
            tuple(pair)
            for pair in brain_candidate_pairs(
                brains, m=4, tau_mA=tau_units
            ).pairs
        }

        assert (0, 1) in pairs, tau_units
        assert (0, 2) not in pairs, tau_units


# --------------------------------------------------------------------------
# D. UI / CLI equivalence for precision
# --------------------------------------------------------------------------


def test_same_precision_resolves_identically_from_cli_and_ui():
    cli = resolve_run_config(
        config_path=PROFILE,
        overrides={"bri.representation_precision_angstrom": 0.002,
                   "duplicate_search.near_duplicate_threshold_angstrom": 0.010,
                   "brain_filter.threshold_angstrom": 0.010},
        override_origin="cli",
    )
    ui = resolve_run_config(
        config_path=PROFILE,
        overrides={"bri.representation_precision_angstrom": 0.002,
                   "duplicate_search.near_duplicate_threshold_angstrom": 0.010,
                   "brain_filter.threshold_angstrom": 0.010},
        override_origin="ui",
    )

    assert cli.scientific_sha256 == ui.scientific_sha256
    assert cli.to_dict() == ui.to_dict()
