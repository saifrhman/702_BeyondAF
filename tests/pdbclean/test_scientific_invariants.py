"""Scientific-invariant regression tests for the configuration refactor.

Two kinds of test live here.

**Invariants.** Properties of the method that no amount of configuration may
change: Brain is a filtering layer, the metric is L-infinity, the criterion is
inclusive, removal requires a direct edge, connected components are not
equivalence classes, and so on.  These run everywhere.

**Frozen oracles.** Checks against the frozen COMP702 20260101 artefacts, which
are large and gitignored.  They skip cleanly when the outputs are not present
(a fresh clone, or CI), and assert exact equality when they are.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from pdbclean.brain import BRAIN_DIMENSION, compute_brain
from pdbclean.brain_prefilter import (
    brain_candidate_pairs,
    brute_force_brain_candidate_pairs,
)
from pdbclean.defaults import (
    FROZEN_PROTOCOL_CONFIG_SHA256,
    FROZEN_REPRESENTATIVE_POLICY_SHA256,
    VALIDATED_DEFAULTS,
    brain_filter_threshold_milliangstrom,
    near_duplicate_threshold_milliangstrom,
)
from pdbclean.duplicate_classification import (
    PAPER_NEAR_DUPLICATE_THRESHOLD_MA,
    classify_bri_distance,
)
from pdbclean.run_provenance import file_sha256
from pdbclean.runconfig import resolve_run_config


REPO_ROOT = Path(__file__).resolve().parents[2]

PROTOCOL_ROOT = (
    REPO_ROOT / "outputs" / "pdbclean" / "20260101" / "protocol3.2-comp702-v1"
)

FROZEN_RELEASE = (
    REPO_ROOT
    / "outputs"
    / "releases"
    / "PDBClean-20260101-protocol3.2-comp702-v1-dedup-v1"
)


def _frozen(path: Path):
    if not path.exists():
        pytest.skip(f"frozen artefact not present in this checkout: {path}")

    return path


def _summary(path: Path) -> dict:
    return json.loads(_frozen(path).read_text(encoding="utf-8"))


# ==========================================================================
# Threshold semantics
# ==========================================================================


def test_production_criterion_is_inclusive():
    """`d_bri_mA <= 10`, not `< 10`.  10 mA is a near duplicate."""

    at_threshold = classify_bri_distance(10)

    assert at_threshold.is_paper_near_duplicate is True
    assert at_threshold.is_nonzero_near_duplicate is True
    assert at_threshold.is_zero_duplicate is False

    assert classify_bri_distance(11).is_paper_near_duplicate is False


def test_exact_duplicate_is_distance_zero_only():
    zero = classify_bri_distance(0)

    assert zero.is_zero_duplicate is True
    assert zero.is_paper_near_duplicate is True
    assert zero.is_nonzero_near_duplicate is False

    assert classify_bri_distance(1).is_zero_duplicate is False


def test_default_threshold_is_ten_milliangstroms():
    assert PAPER_NEAR_DUPLICATE_THRESHOLD_MA == 10
    assert near_duplicate_threshold_milliangstrom(VALIDATED_DEFAULTS) == 10
    assert brain_filter_threshold_milliangstrom(VALIDATED_DEFAULTS) == 10


def test_no_argument_call_still_reproduces_the_frozen_behaviour():
    """The refactor added a keyword; the default must be unchanged."""

    for distance in range(0, 25):
        assert classify_bri_distance(distance) == classify_bri_distance(
            distance, threshold_mA=PAPER_NEAR_DUPLICATE_THRESHOLD_MA
        )


def test_configured_threshold_is_honoured_and_still_inclusive():
    assert classify_bri_distance(
        5, threshold_mA=5
    ).is_paper_near_duplicate is True
    assert classify_bri_distance(
        6, threshold_mA=5
    ).is_paper_near_duplicate is False


def test_representation_precision_is_not_the_threshold():
    """0.001 A is how complete BRI is represented, not what it compares to."""

    resolved = resolve_run_config()

    assert resolved.get("bri.representation_precision_angstrom") == 0.001
    assert resolved.get(
        "duplicate_search.near_duplicate_threshold_angstrom"
    ) == 0.010
    assert resolved.near_duplicate_threshold_mA == 10


# ==========================================================================
# Brain: filtering only, 9-D, exact-length buckets
# ==========================================================================


def test_brain_is_nine_dimensional():
    assert BRAIN_DIMENSION == 9

    rng = np.random.default_rng(0)
    bri = np.around(rng.random((12, 9)), 3)

    assert compute_brain(bri).shape == (9,)


def test_brain_excludes_the_first_row():
    bri = np.around(
        np.tile(np.arange(1, 5, dtype=np.float64).reshape(-1, 1), (1, 9)), 3
    )

    expected = bri[1:].mean(axis=0)

    assert np.allclose(compute_brain(bri), expected)


def test_brain_is_undefined_for_single_residue_backbones():
    with pytest.raises(ValueError):
        compute_brain(np.zeros((1, 9)))


def test_brain_is_configured_as_a_filtering_layer_only():
    resolved = resolve_run_config()

    assert resolved.get("brain.role") == "filtering_and_indexing_only"
    assert resolved.get(
        "duplicate_search.final_classification_basis"
    ) == "complete_BRI"


def test_brain_prefilter_matches_the_brute_force_oracle():
    rng = np.random.default_rng(17)

    m = 6
    brains = np.around(rng.random((60, 9)) * 0.05, 3)
    # A deliberate near-collision so the bucket is not trivially empty.
    brains[3] = brains[2]
    brains[7] = np.around(brains[6] + 0.000_5, 3)

    fast = brain_candidate_pairs(brains, m=m)
    oracle = brute_force_brain_candidate_pairs(brains, m=m)

    assert np.array_equal(fast.pairs, oracle.pairs)


def test_brain_prefilter_default_tau_is_the_frozen_tau():
    rng = np.random.default_rng(23)

    m = 4
    brains = np.around(rng.random((40, 9)) * 0.02, 3)

    assert np.array_equal(
        brain_candidate_pairs(brains, m=m).pairs,
        brain_candidate_pairs(brains, m=m, tau_mA=10).pairs,
    )


def test_brain_prefilter_radius_scales_with_the_bucket_length():
    """tau is exactly tau_mA * (m - 1) integer milliangstrom-sum units."""

    brains = np.array(
        [
            [0.000] * 9,
            [0.010] * 9,  # exactly at tau
            [0.011] * 9,  # just outside
        ]
    )

    for m in (2, 5, 40):
        pairs = brain_candidate_pairs(brains, m=m).pairs

        assert (0, 1) in {tuple(pair) for pair in pairs}
        assert (0, 2) not in {tuple(pair) for pair in pairs}


def test_candidate_pairs_never_cross_length_buckets():
    """A bucket call is scoped to one exact m by construction."""

    brains = np.around(np.zeros((3, 9)), 3)

    for m in (2, 3):
        result = brain_candidate_pairs(brains, m=m)
        assert result.m == m


# ==========================================================================
# Configuration-level invariants (section 22 of the productisation brief)
# ==========================================================================


def test_grouping_is_exact_chain_length():
    assert resolve_run_config().get(
        "brain_filter.grouping"
    ) == "exact_chain_length"


def test_model_scope_is_model_one():
    resolved = resolve_run_config()

    assert resolved.get("selection.models.policy") == "first_model"
    assert resolved.get("selection.models.model_id") == 1


def test_components_are_not_equivalence_classes():
    resolved = resolve_run_config()

    assert resolved.get(
        "graph.connected_component_is_duplicate_equivalence"
    ) is False


def test_removal_requires_a_direct_edge():
    resolved = resolve_run_config()

    assert resolved.get("graph.require_direct_edge_for_removal") is True
    assert resolved.get("graph.automatic_transitive_removal") is False


def test_m1_chains_are_retained():
    assert resolve_run_config().get(
        "representative_selection.m1_policy"
    ) == "retain_all"


def test_deduplication_starts_at_length_two():
    assert resolve_run_config().get(
        "representative_selection.minimum_deduplicated_chain_length"
    ) == 2


def test_historical_geometric_search_block_is_not_a_default():
    """The PDB707K `L_inf < 1.0` design is superseded, not silently applied."""

    assert "geometric_search" not in VALIDATED_DEFAULTS

    resolved = resolve_run_config()

    assert resolved.get("geometric_search") is None
    assert resolved.get("duplicate_search.operator") == "less_than_or_equal"


# ==========================================================================
# Frozen configuration files remain byte-identical
# ==========================================================================


def test_frozen_protocol_config_is_unmodified():
    path = REPO_ROOT / "config" / "pdbclean" / "protocol_3_2_comp702_v1.yaml"

    if not path.is_file():
        pytest.skip("frozen protocol configuration not present")

    assert file_sha256(path) == FROZEN_PROTOCOL_CONFIG_SHA256


def test_frozen_representative_policy_is_unmodified():
    path = (
        REPO_ROOT
        / "config"
        / "pdbclean"
        / "stage14_representative_policy_v1.yaml"
    )

    if not path.is_file():
        pytest.skip("frozen representative policy not present")

    assert file_sha256(path) == FROZEN_REPRESENTATIVE_POLICY_SHA256


# ==========================================================================
# Frozen 20260101 oracles
# ==========================================================================


def test_frozen_classification_summary_matches_the_default_threshold():
    summary = _summary(
        PROTOCOL_ROOT / "duplicate_classification" / "global_summary.json"
    )

    resolved = resolve_run_config(
        overrides=["snapshot.mode=fixed", "snapshot.snapshot_id=20260101"]
    )

    assert summary["paper_near_duplicate_threshold_mA"] == 10
    assert summary["paper_near_duplicate_threshold_angstrom"] == 0.010
    assert summary["classification_basis"] == (
        "exact_full_bri_integer_milliangstrom"
    )
    assert resolved.near_duplicate_threshold_mA == 10


def test_frozen_classification_counts_are_the_documented_population():
    summary = _summary(
        PROTOCOL_ROOT / "duplicate_classification" / "global_summary.json"
    )

    assert summary["paper_near_duplicate_pair_count"] == 1_072_751
    assert summary["zero_duplicate_pair_count"] == 17_373
    assert summary["nonzero_near_duplicate_pair_count"] == 1_055_378
    assert (
        summary["zero_duplicate_pair_count"]
        + summary["nonzero_near_duplicate_pair_count"]
        == summary["paper_near_duplicate_pair_count"]
    )


def test_frozen_brain_prefilter_threshold_is_ten_milliangstroms():
    summary = _summary(
        PROTOCOL_ROOT / "brain_prefilter" / "global_summary.json"
    )

    assert summary["brain_threshold_angstrom"] == 0.010
    assert summary["metric"] == "L_infinity"
    assert summary["search_implementation"] == (
        "scipy_cKDTree_p_inf_eps_0_exact_integer_postfilter"
    )
    assert resolve_run_config().brain_threshold_mA == 10


def test_frozen_graph_threshold_string_is_inclusive():
    summary = _summary(
        PROTOCOL_ROOT / "stage14_geometric_graph" / "global_summary.json"
    )

    assert summary["threshold"] == "d_bri_mA <= 10"


def test_frozen_release_population_closes():
    summary = _summary(
        PROTOCOL_ROOT
        / "stage14_representative_selection_v1"
        / "global_summary.json"
    )

    retained = summary["final_retained_chain_count"]
    removed = summary["final_removed_chain_count"]

    assert retained == 499_770
    assert removed == 78_754
    assert retained + removed == 578_524


def test_frozen_selection_used_the_inclusive_criterion():
    summary = _summary(
        PROTOCOL_ROOT
        / "stage14_representative_selection_v1"
        / "global_summary.json"
    )

    assert summary["duplicate_threshold"] == "d_bri_mA <= 10"
    assert summary["policy_name"] == resolve_run_config().get(
        "representative_selection.policy_name"
    )
    assert summary["policy_config_sha256"] == (
        FROZEN_REPRESENTATIVE_POLICY_SHA256
    )


def test_frozen_selection_removed_only_chains_with_a_direct_edge():
    summary = _summary(
        PROTOCOL_ROOT
        / "stage14_representative_selection_v1"
        / "global_summary.json"
    )

    assert summary["every_removed_chain_has_direct_representative_edge"] is True


def test_frozen_selection_did_not_treat_components_as_equivalence_classes():
    summary = _summary(
        PROTOCOL_ROOT
        / "stage14_representative_selection_v1"
        / "global_summary.json"
    )

    assert summary["connectedness_treated_as_equivalence"] is False

    # Non-clique components with more than one representative are the direct
    # evidence that connectedness was not collapsed to one survivor.
    assert summary["nonclique_components_with_multiple_representatives"] > 0
    assert summary["nonclique_representatives"] > summary["nonclique_components"]


def test_frozen_selection_retained_every_m1_chain():
    summary = _summary(
        PROTOCOL_ROOT
        / "stage14_representative_selection_v1"
        / "global_summary.json"
    )

    assert summary["m1_deduplication_performed"] is False
    assert summary["retained_m1_total"] == 764

    prefilter = _summary(
        PROTOCOL_ROOT / "brain_prefilter" / "global_summary.json"
    )

    assert prefilter["m1_bypass_chain_count"] == summary["retained_m1_total"]


def test_frozen_selection_did_not_use_stage13_as_the_global_edge_set():
    summary = _summary(
        PROTOCOL_ROOT
        / "stage14_representative_selection_v1"
        / "global_summary.json"
    )

    assert summary["stage13_review_subset_used_as_global_edge_set"] is False


def test_frozen_selection_did_not_compare_against_an_old_snapshot():
    summary = _summary(
        PROTOCOL_ROOT
        / "stage14_representative_selection_v1"
        / "global_summary.json"
    )

    assert summary["old_snapshot_comparison_used"] is False


def test_frozen_selection_accounting_closes_on_every_axis():
    summary = _summary(
        PROTOCOL_ROOT
        / "stage14_representative_selection_v1"
        / "global_summary.json"
    )

    assert summary["action_histogram"]["remove"] == 78_754
    assert summary["action_histogram"]["retain_representative"] == 21_100

    assert (
        summary["retained_mge2_representatives"]
        + summary["retained_mge2_no_edge_chains"]
        == summary["retained_mge2_total"]
    )
    assert (
        summary["retained_mge2_total"] + summary["retained_m1_total"]
        == summary["final_retained_chain_count"]
    )
    assert summary["removed_mge2_chains"] == summary["final_removed_chain_count"]


def test_frozen_pair_accounting_closes():
    summary = _summary(
        PROTOCOL_ROOT / "duplicate_classification" / "global_summary.json"
    )

    assert summary["pair_accounting_valid"] is True
    assert summary["input_pair_count"] == 3_531_895
    assert summary["not_near_duplicate_pair_count"] == 2_459_144
    assert (
        summary["paper_near_duplicate_pair_count"]
        + summary["not_near_duplicate_pair_count"]
        == summary["input_pair_count"]
    )


def test_frozen_brain_prefilter_population_closes():
    summary = _summary(
        PROTOCOL_ROOT / "brain_prefilter" / "global_summary.json"
    )

    assert summary["population_accounting_valid"] is True
    assert summary["input_chain_count"] == 578_524
    assert summary["brain_defined_chain_count"] == 577_760
    assert (
        summary["brain_defined_chain_count"] + summary["m1_bypass_chain_count"]
        == summary["input_chain_count"]
    )


def test_frozen_release_manifest_agrees_with_the_selection_summary():
    manifest = _summary(FROZEN_RELEASE / "release_manifest.json")

    assert manifest["retained_chain_count"] == 499_770
    assert manifest["removed_chain_count"] == 78_754
    assert manifest["canonical_input_chain_count"] == 578_524
    assert manifest["m1_retained_chain_count"] == 764
    assert manifest["m1_removed_chain_count"] == 0


def test_frozen_release_manifest_records_the_method_it_used():
    manifest = _summary(FROZEN_RELEASE / "release_manifest.json")

    assert manifest["snapshot"] == "20260101"
    assert manifest["model_scope"] == "model_1"
    assert manifest["near_duplicate_threshold"] == "d_bri_mA <= 10"
    assert manifest["automatic_transitive_removal"] is False
    assert manifest["connectedness_treated_as_duplicate_equivalence"] is False
    assert manifest["every_removed_chain_has_direct_representative_edge"] is True
    assert manifest["m1_deduplication_performed"] is False
    assert manifest["old_snapshot_comparison_used"] is False
    assert manifest["stage13_review_subset_used_as_global_edge_set"] is False


def test_release_manifest_semantics_match_the_validated_defaults():
    """What the frozen run recorded is what the defaults still resolve to."""

    manifest = _summary(FROZEN_RELEASE / "release_manifest.json")

    resolved = resolve_run_config(
        overrides=["snapshot.mode=fixed", "snapshot.snapshot_id=20260101"]
    )

    assert manifest["near_duplicate_threshold"] == (
        f"d_bri_mA <= {resolved.near_duplicate_threshold_mA}"
    )
    assert manifest["model_scope"] == (
        f"model_{resolved.get('selection.models.model_id')}"
    )
    assert manifest["automatic_transitive_removal"] == resolved.get(
        "graph.automatic_transitive_removal"
    )
    assert manifest["connectedness_treated_as_duplicate_equivalence"] == (
        resolved.get("graph.connected_component_is_duplicate_equivalence")
    )
    assert manifest["protocol"] == resolved.get("release.protocol_version")
    assert manifest["snapshot"] == resolved.get("snapshot.snapshot_id")
    assert manifest["representative_policy_version"] == resolved.get(
        "representative_selection.policy_version"
    )
    assert manifest["representative_policy_sha256"] == (
        FROZEN_REPRESENTATIVE_POLICY_SHA256
    )


def test_frozen_release_is_reachable_from_the_default_configuration():
    """The validated defaults plus the frozen snapshot name the frozen release."""

    from pdbclean.pipeline import PipelinePaths

    resolved = resolve_run_config(
        overrides=["snapshot.mode=fixed", "snapshot.snapshot_id=20260101"]
    )

    paths = PipelinePaths.from_config(resolved, repo_root=REPO_ROOT)

    assert paths.release == FROZEN_RELEASE.name


# ==========================================================================
# The superseded historical `geometric_search` block must be inert
# ==========================================================================
#
# `config/pdbclean/protocol_3_2_comp702_v1.yaml` is frozen byte-for-byte
# because its SHA256 is embedded in frozen provenance. It carries a historical
# PDB707K `geometric_search` block:
#
#     query_radius: 1.0
#     retain_when: {operator: less_than, value: 1.0}
#
# That is NOT the production classification semantics. These tests prove the
# generic runner cannot activate it, without editing the frozen file.


FROZEN_PROTOCOL_CONFIG = (
    REPO_ROOT / "config" / "pdbclean" / "protocol_3_2_comp702_v1.yaml"
)


def _frozen_protocol_yaml() -> dict:
    import yaml

    return yaml.safe_load(
        _frozen(FROZEN_PROTOCOL_CONFIG).read_text(encoding="utf-8")
    )


def test_the_stale_block_really_is_present_in_the_frozen_file():
    """Guard the guard: if this ever fails the premise has changed."""

    block = _frozen_protocol_yaml()["geometric_search"]

    assert block["query_radius"] == 1.0
    assert block["retain_when"] == {"operator": "less_than", "value": 1.0}


def test_stale_block_is_not_in_the_validated_defaults():
    assert "geometric_search" not in VALIDATED_DEFAULTS


def test_stale_block_is_not_in_a_resolved_configuration():
    resolved = resolve_run_config(
        config_path=REPO_ROOT
        / "config"
        / "pdbclean"
        / "profiles"
        / "comp702_frozen_20260101.yaml"
    )

    assert "geometric_search" not in resolved.to_dict()
    assert resolved.get("geometric_search") is None


def test_query_radius_one_angstrom_is_nowhere_in_a_resolved_configuration():
    """1.0 A must not appear as any threshold the runner would use."""

    resolved = resolve_run_config()

    assert resolved.get(
        "duplicate_search.near_duplicate_threshold_angstrom"
    ) != 1.0
    assert resolved.get("brain_filter.threshold_angstrom") != 1.0

    assert resolved.near_duplicate_threshold_mA == 10
    assert resolved.near_duplicate_threshold_mA != 1000


def test_final_classification_resolves_to_the_explicit_duplicate_search():
    """The threshold comes from duplicate_search, never from the stale block."""

    resolved = resolve_run_config()

    assert resolved.source_of(
        "duplicate_search.near_duplicate_threshold_angstrom"
    ) == "builtin_default"
    assert resolved.get(
        "duplicate_search.near_duplicate_threshold_angstrom"
    ) == 0.010
    assert resolved.get("duplicate_search.operator") == "less_than_or_equal"


def test_less_than_semantics_cannot_be_activated():
    """`retain_when: less_than` is exactly what the validator refuses."""

    from pdbclean.runconfig import RunConfigError

    with pytest.raises(RunConfigError, match="inclusive"):
        resolve_run_config(overrides=["duplicate_search.operator=less_than"])


def test_complete_bri_remains_the_final_classification_basis():
    from pdbclean.runconfig import RunConfigError

    resolved = resolve_run_config()

    assert resolved.get(
        "duplicate_search.final_classification_basis"
    ) == "complete_BRI"
    assert resolved.get("duplicate_search.metric") == "L_infinity"

    with pytest.raises(RunConfigError):
        resolve_run_config(
            overrides=["duplicate_search.final_classification_basis=brain"]
        )


def test_loading_the_frozen_file_as_a_layer_still_cannot_activate_it():
    """Even the worst case -- layering the frozen file itself -- is inert."""

    resolved = resolve_run_config(config_path=FROZEN_PROTOCOL_CONFIG)

    # The stale block rides along as inert data ...
    assert resolved.get("geometric_search.query_radius") == 1.0

    # ... but nothing the runner classifies with is touched by it.
    assert resolved.near_duplicate_threshold_mA == 10
    assert resolved.get("duplicate_search.operator") == "less_than_or_equal"
    assert resolved.get(
        "duplicate_search.final_classification_basis"
    ) == "complete_BRI"

    # And it is excluded from the scientific identity, so it can never make a
    # scientifically identical run look different.
    from pdbclean.runconfig import scientific_projection

    assert "geometric_search" not in scientific_projection(resolved.data)
    assert resolved.scientific_sha256 == resolve_run_config(
        overrides=[
            f"release.protocol_version={resolved.get('release.protocol_version')}"
        ]
    ).scientific_sha256


def test_the_stage_that_classifies_never_receives_the_stale_block():
    """Stage 10's argv carries the resolved threshold, not query_radius."""

    from pdbclean.cli import stage_command
    from pdbclean.pipeline import PipelinePaths

    resolved = resolve_run_config(
        overrides=["snapshot.mode=fixed", "snapshot.snapshot_id=20260101"]
    )
    paths = PipelinePaths.from_config(resolved, repo_root=REPO_ROOT)

    for stage_id in ("candidate_filtering", "complete_bri_nn",
                     "duplicate_classification", "redundancy_graph"):
        argv = stage_command(stage_id, resolved, paths)

        assert argv is not None
        assert "1.0" not in argv
        assert "--query-radius" not in argv

    graph = stage_command("redundancy_graph", resolved, paths)

    assert "--threshold-mA" in graph
    assert graph[graph.index("--threshold-mA") + 1] == "10"
