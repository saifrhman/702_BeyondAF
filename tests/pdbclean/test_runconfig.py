"""Configuration-resolution regression tests.

These tests pin the *scientific* behaviour of the configuration layer:

* the built-in defaults reproduce the validated COMP702 methodology;
* precedence is deterministic (defaults -> file -> explicit override);
* the resolved configuration, not the input file, is what gets hashed;
* the scientific hash is stable across infrastructure differences and moves
  whenever a scientific parameter moves;
* configurations that would silently change the method are rejected.
"""

from __future__ import annotations

import pytest
import yaml

from pdbclean import defaults as defaults_module
from pdbclean.runconfig import (
    LAYER_BUILTIN,
    RunConfigError,
    config_sha256,
    resolve_run_config,
    scientific_projection,
    scientific_sha256,
    validate_resolved_config,
    runtime_environment,
    with_resolved_snapshot,
)


# --------------------------------------------------------------------------
# Validated defaults
# --------------------------------------------------------------------------


def test_defaults_reproduce_validated_methodology():
    resolved = resolve_run_config()

    assert resolved.get("selection.models.policy") == "first_model"
    assert resolved.get("selection.models.model_id") == 1

    assert resolved.get(
        "quality_rules.backbone_distance.minimum_distance_angstrom"
    ) == 0.01
    assert resolved.get(
        "post_cleaning_geometric_validation.minimum_triangle_angle_degrees"
    ) == 3.0

    assert resolved.get("bri.implementation_version") == "1.2.2"
    assert resolved.get("bri.representation_precision_angstrom") == 0.001

    assert resolved.get("brain.dimension") == 9
    assert resolved.get("brain.definition") == "average_complete_BRI"
    assert resolved.get("brain.role") == "filtering_and_indexing_only"

    assert resolved.get("brain_filter.threshold_angstrom") == 0.010
    assert resolved.get("brain_filter.grouping") == "exact_chain_length"
    assert resolved.get("brain_filter.operator") == "less_than_or_equal"

    assert resolved.get(
        "duplicate_search.near_duplicate_threshold_angstrom"
    ) == 0.010
    assert resolved.get("duplicate_search.metric") == "L_infinity"
    assert resolved.get(
        "duplicate_search.final_classification_basis"
    ) == "complete_BRI"
    assert resolved.get("duplicate_search.operator") == "less_than_or_equal"
    assert resolved.get(
        "duplicate_search.exact_duplicate_criterion"
    ) == "d_bri_mA == 0"

    assert resolved.get("graph.require_direct_edge_for_removal") is True
    assert resolved.get(
        "graph.connected_component_is_duplicate_equivalence"
    ) is False
    assert resolved.get("graph.automatic_transitive_removal") is False

    assert resolved.get("representative_selection.m1_policy") == "retain_all"
    assert resolved.get(
        "representative_selection.minimum_deduplicated_chain_length"
    ) == 2


def test_default_thresholds_are_exact_integer_milliangstroms():
    resolved = resolve_run_config()

    assert resolved.near_duplicate_threshold_mA == 10
    assert resolved.brain_threshold_mA == 10


def test_default_snapshot_mode_is_latest_complete():
    resolved = resolve_run_config()

    assert resolved.get("snapshot.mode") == "latest_complete"
    assert resolved.get("snapshot.snapshot_id") is None


def test_defaults_declare_no_expectation_counts():
    """Population counts belong to a snapshot, never to the built-ins."""

    resolved = resolve_run_config()

    expectations = resolved.get("expectations") or {}

    assert expectations
    assert all(value is None for value in expectations.values())


def test_defaults_are_not_mutated_by_resolution():
    first = resolve_run_config(
        overrides=["duplicate_search.near_duplicate_threshold_angstrom=0.005"]
    )

    assert first.near_duplicate_threshold_mA == 5

    second = resolve_run_config()

    assert second.near_duplicate_threshold_mA == 10
    assert (
        defaults_module.VALIDATED_DEFAULTS["duplicate_search"][
            "near_duplicate_threshold_angstrom"
        ]
        == 0.010
    )


# --------------------------------------------------------------------------
# Precedence
# --------------------------------------------------------------------------


def _write(tmp_path, name, payload):
    path = tmp_path / name
    path.write_text(yaml.safe_dump(payload), encoding="utf-8")
    return path


def test_precedence_default_then_file_then_override(tmp_path):
    config = _write(
        tmp_path,
        "profile.yaml",
        {
            "snapshot": {"mode": "fixed", "snapshot_id": "20250101"},
            "selection": {"models": {"model_id": 1}},
        },
    )

    resolved = resolve_run_config(
        config_path=config,
        overrides=["snapshot.snapshot_id=20260101"],
    )

    assert resolved.get("snapshot.snapshot_id") == "20260101"
    assert resolved.source_of("snapshot.snapshot_id").startswith("override:")

    assert resolved.get("snapshot.mode") == "fixed"
    assert resolved.source_of("snapshot.mode").startswith("config_file:")

    assert resolved.get("brain.dimension") == 9
    assert resolved.source_of("brain.dimension") == LAYER_BUILTIN


def test_untouched_values_keep_their_validated_defaults(tmp_path):
    config = _write(
        tmp_path,
        "profile.yaml",
        {"execution": {"executor": "slurm"}},
    )

    resolved = resolve_run_config(config_path=config)

    assert resolved.near_duplicate_threshold_mA == 10
    assert resolved.brain_threshold_mA == 10
    assert resolved.get("selection.models.model_id") == 1


def test_mappings_merge_recursively(tmp_path):
    config = _write(
        tmp_path,
        "profile.yaml",
        {"quality_rules": {"backbone_distance": {
            "minimum_distance_angstrom": 0.02
        }}},
    )

    resolved = resolve_run_config(config_path=config)

    assert resolved.get(
        "quality_rules.backbone_distance.minimum_distance_angstrom"
    ) == 0.02
    # Sibling rules survive the partial layer, and so do sibling leaves.
    assert resolved.get("quality_rules.entry_protein.rule_id") == "Q001"
    assert resolved.get("quality_rules.backbone_distance.rule_id") == "Q005"


def test_mapping_overrides_accepted_as_dict():
    resolved = resolve_run_config(
        overrides={"snapshot.snapshot_id": "20260101"},
        override_origin="ui",
    )

    assert resolved.get("snapshot.snapshot_id") == "20260101"
    assert resolved.source_of("snapshot.snapshot_id") == "override:ui"


def test_unquoted_snapshot_id_selects_the_same_snapshot(tmp_path):
    """YAML reads an unquoted date as an integer; it is the same snapshot."""

    config = _write(
        tmp_path,
        "profile.yaml",
        {"snapshot": {"mode": "fixed", "snapshot_id": 20260101}},
    )

    from_file = resolve_run_config(config_path=config)
    from_cli = resolve_run_config(
        overrides=["snapshot.mode=fixed", "snapshot.snapshot_id=20260101"]
    )
    quoted = resolve_run_config(
        overrides={"snapshot.mode": "fixed", "snapshot.snapshot_id": "20260101"}
    )

    assert from_file.get("snapshot.snapshot_id") == "20260101"
    assert from_cli.get("snapshot.snapshot_id") == "20260101"
    assert (
        from_file.scientific_sha256
        == from_cli.scientific_sha256
        == quoted.scientific_sha256
    )


def test_missing_config_file_is_an_error(tmp_path):
    with pytest.raises(RunConfigError):
        resolve_run_config(config_path=tmp_path / "absent.yaml")


# --------------------------------------------------------------------------
# Hashing and run identity
# --------------------------------------------------------------------------


def test_config_hash_is_order_independent():
    assert config_sha256({"a": 1, "b": 2}) == config_sha256({"b": 2, "a": 1})


def test_canonical_config_keeps_the_runtime_template(monkeypatch):
    """``${TMPDIR}/pdbclean`` is a policy, identical on every host."""

    monkeypatch.setenv("TMPDIR", "/tmp/node-a")

    resolved = resolve_run_config()

    assert resolved.get("storage.temporary_root") == "${TMPDIR}/pdbclean"


def test_canonical_hash_is_stable_across_machines(monkeypatch):
    """A compute node's $TMPDIR must not change the run's identity."""

    monkeypatch.setenv("TMPDIR", "/tmp/node-a")
    first = resolve_run_config()

    monkeypatch.setenv("TMPDIR", "/scratch/node-b/999")
    second = resolve_run_config()

    assert first.get("storage.temporary_root") == second.get(
        "storage.temporary_root"
    )
    assert first.sha256 == second.sha256
    assert first.scientific_sha256 == second.scientific_sha256


def test_non_runtime_variables_still_expand(monkeypatch):
    """Only host-runtime variables are deferred; the rest resolve as before."""

    monkeypatch.setenv("PDBCLEAN_TEST_ROOT", "/data/example")

    resolved = resolve_run_config(
        overrides={"storage.output_root": "${PDBCLEAN_TEST_ROOT}/out"}
    )

    assert resolved.get("storage.output_root") == "/data/example/out"


def test_runtime_environment_resolves_the_template_for_execution(monkeypatch):
    monkeypatch.setenv("TMPDIR", "/scratch/node-b/999")

    resolved = resolve_run_config()
    runtime = runtime_environment(resolved)

    assert runtime["storage"]["temporary_root"] == "/scratch/node-b/999/pdbclean"
    assert runtime["environment"]["TMPDIR"] == "/scratch/node-b/999"
    assert runtime["hostname"]

    # Resolving it for execution must not mutate the canonical configuration.
    assert resolved.get("storage.temporary_root") == "${TMPDIR}/pdbclean"


def test_protocol_projection_resolves_runtime_templates(monkeypatch):
    """The execution handoff -- and only it -- gets concrete host paths."""

    monkeypatch.setenv("TMPDIR", "/scratch/node-b/999")

    resolved = resolve_run_config(
        overrides=["snapshot.mode=fixed", "snapshot.snapshot_id=20260101"]
    )

    projected = resolved.to_protocol_config()

    assert projected["storage"]["temporary_root"] == (
        "/scratch/node-b/999/pdbclean"
    )
    assert resolved.get("storage.temporary_root") == "${TMPDIR}/pdbclean"
    assert resolved.sha256 == resolve_run_config(
        overrides=["snapshot.mode=fixed", "snapshot.snapshot_id=20260101"]
    ).sha256


def test_scientific_hash_is_stable_across_machines(monkeypatch):
    monkeypatch.setenv("TMPDIR", "/tmp/node-a")
    first = resolve_run_config()

    monkeypatch.setenv("TMPDIR", "/scratch/node-b/999")
    second = resolve_run_config()

    assert first.scientific_sha256 == second.scientific_sha256


def test_scientific_hash_ignores_executor_choice():
    baseline = resolve_run_config()
    slurm = resolve_run_config(overrides=["execution.executor=slurm"])

    assert baseline.scientific_sha256 == slurm.scientific_sha256


@pytest.mark.parametrize(
    "override",
    [
        "duplicate_search.near_duplicate_threshold_angstrom=0.005",
        "brain_filter.threshold_angstrom=0.020",
        "quality_rules.backbone_distance.minimum_distance_angstrom=0.02",
        "post_cleaning_geometric_validation.minimum_triangle_angle_degrees=5.0",
        "selection.models.model_id=2",
        "representative_selection.minimum_deduplicated_chain_length=3",
    ],
)
def test_scientific_hash_moves_when_the_science_moves(override):
    baseline = resolve_run_config()
    changed = resolve_run_config(overrides=[override], validate=False)

    assert baseline.scientific_sha256 != changed.scientific_sha256


def test_scientific_hash_tracks_the_snapshot_identity():
    baseline = resolve_run_config(overrides=["snapshot.snapshot_id=20260101"])
    other = resolve_run_config(overrides=["snapshot.snapshot_id=20250101"])

    assert baseline.scientific_sha256 != other.scientific_sha256


def test_scientific_hash_ignores_how_the_snapshot_was_selected():
    """Latest-complete and explicit selection of the same snapshot agree."""

    explicit = resolve_run_config(
        overrides=["snapshot.mode=fixed", "snapshot.snapshot_id=20260101"]
    )

    discovered = with_resolved_snapshot(
        resolve_run_config(),
        snapshot_id="20260101",
        selection_mode="latest_complete",
    )

    assert explicit.scientific_sha256 == discovered.scientific_sha256


def test_scientific_projection_excludes_infrastructure():
    resolved = resolve_run_config()

    projection = scientific_projection(resolved.data)

    for section in ("storage", "execution", "observability", "expectations"):
        assert section not in projection

    for section in ("bri", "brain", "brain_filter", "duplicate_search",
                    "graph", "representative_selection", "quality_rules"):
        assert section in projection

    assert scientific_sha256(resolved.data) == config_sha256(projection)


def test_scientific_projection_keeps_only_snapshot_identity():
    resolved = resolve_run_config(overrides=["snapshot.snapshot_id=20260101"])

    projection = scientific_projection(resolved.data)

    assert projection["snapshot"] == {
        "snapshot_id": "20260101",
        "bucket_url": resolved.get("snapshot.bucket_url"),
    }


def test_expectations_do_not_change_the_scientific_identity():
    """Expected counts are assertions about a run, not parameters of it."""

    baseline = resolve_run_config(overrides=["snapshot.snapshot_id=20260101"])
    gated = resolve_run_config(
        overrides=[
            "snapshot.snapshot_id=20260101",
            "expectations.retained_chain_count=499770",
        ]
    )

    assert baseline.scientific_sha256 == gated.scientific_sha256
    assert baseline.sha256 != gated.sha256


# --------------------------------------------------------------------------
# Guard rails
# --------------------------------------------------------------------------


def test_non_inclusive_operator_is_rejected():
    with pytest.raises(RunConfigError, match="inclusive"):
        resolve_run_config(overrides=["duplicate_search.operator=less_than"])


def test_brain_dimension_is_fixed_by_definition():
    with pytest.raises(RunConfigError, match="9"):
        resolve_run_config(overrides=["brain.dimension=8"])


def test_metric_must_stay_l_infinity():
    with pytest.raises(RunConfigError, match="L_infinity"):
        resolve_run_config(overrides=["duplicate_search.metric=euclidean"])


def test_final_classification_basis_must_stay_complete_bri():
    with pytest.raises(RunConfigError, match="complete_BRI"):
        resolve_run_config(
            overrides=["duplicate_search.final_classification_basis=brain"]
        )


def test_brain_filter_may_not_be_lossy():
    """A Brain prefilter tighter than the classifier would drop true pairs."""

    with pytest.raises(RunConfigError):
        resolve_run_config(
            overrides=[
                "brain_filter.threshold_angstrom=0.005",
                "duplicate_search.near_duplicate_threshold_angstrom=0.010",
            ]
        )


def test_thresholds_must_land_on_the_representation_grid():
    with pytest.raises(RunConfigError):
        resolve_run_config(
            overrides=[
                "duplicate_search.near_duplicate_threshold_angstrom=0.0105"
            ]
        )


def test_validate_accepts_the_validated_defaults():
    validate_resolved_config(resolve_run_config(validate=False).data)


# --------------------------------------------------------------------------
# Projection onto the legacy protocol schema
# --------------------------------------------------------------------------


def test_protocol_projection_requires_a_pinned_snapshot():
    with pytest.raises(RunConfigError, match="snapshot"):
        resolve_run_config().to_protocol_config()


def test_protocol_projection_pins_the_snapshot_as_fixed():
    resolved = with_resolved_snapshot(
        resolve_run_config(),
        snapshot_id="20260101",
        selection_mode="latest_complete",
    )

    projected = resolved.to_protocol_config()

    assert projected["snapshot"]["mode"] == "fixed"
    assert projected["snapshot"]["snapshot_id"] == "20260101"


def test_protocol_projection_carries_the_resolved_threshold():
    resolved = with_resolved_snapshot(
        resolve_run_config(),
        snapshot_id="20260101",
    )

    projected = resolved.to_protocol_config()

    assert projected["duplicate_search"][
        "near_duplicate_threshold_angstrom"
    ] == 0.010
    assert projected["graph"]["require_direct_edge_for_removal"] is True


def test_resolved_yaml_round_trips():
    resolved = resolve_run_config(overrides=["snapshot.snapshot_id=20260101"])

    reloaded = yaml.safe_load(resolved.to_yaml())

    assert config_sha256(reloaded) == resolved.sha256
