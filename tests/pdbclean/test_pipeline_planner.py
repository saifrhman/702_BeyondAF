"""Pipeline-orchestration regression tests.

The orchestrator is plumbing, but two of its properties are scientific
safeguards and are pinned here:

* restartability is decided from configuration and stage summaries, never from
  directory names;
* an output produced under a different scientific configuration is never
  silently reused.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from pdbclean import pipeline as pipeline_module
from pdbclean.pipeline import (
    ACTION_BLOCKED,
    ACTION_NOT_APPLICABLE,
    ACTION_REUSE,
    ACTION_RUN,
    BLOCKED,
    COMPLETE,
    EXECUTION_COMPLETE,
    NOT_APPLICABLE,
    PARTIAL,
    PENDING,
    VALIDATION_FAIL,
    VALIDATION_PASS,
    DryRunExecutor,
    PipelineError,
    PipelinePaths,
    build_executor,
    gold_release_summary,
    inspect_stage,
    plan_pipeline,
)
from pdbclean.runconfig import resolve_run_config
from pdbclean.stage_registry import (
    LAYER_BRONZE,
    LAYER_GOLD,
    LAYER_SILVER,
    LAYER_SNAPSHOT,
    STAGES,
    stage_order,
    stages_for_layer,
)


def _config(tmp_path, *extra):
    return resolve_run_config(
        overrides=[
            "snapshot.mode=fixed",
            "snapshot.snapshot_id=20260101",
            f"storage.output_root={tmp_path / 'outputs'}",
            f"storage.release_root={tmp_path / 'releases'}",
            f"storage.run_root={tmp_path / 'runs'}",
            *extra,
        ]
    )


def _spec(stage_id):
    for stage in STAGES:
        if stage.stage_id == stage_id:
            return stage

    raise AssertionError(f"no such stage: {stage_id}")


def _materialise(paths, stage_id, summary, *, success=True, output=True):
    stage = _spec(stage_id)
    directory = paths.stage_directory(stage)
    directory.mkdir(parents=True, exist_ok=True)

    if stage.summary:
        (directory / stage.summary).write_text(
            json.dumps(summary), encoding="utf-8"
        )

    if output and stage.primary_output:
        target = directory / stage.primary_output
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(b"artefact")

    if success and stage.success_marker:
        (directory / stage.success_marker).write_text("", encoding="utf-8")

    return directory


# --------------------------------------------------------------------------
# Stage registry
# --------------------------------------------------------------------------


def test_every_stage_1_to_14_is_present():
    order = stage_order()

    assert len(order) == 15  # snapshot selection + Stages 1-14
    assert order[0] == "snapshot"
    assert order[-1] == "gold_release"
    assert len(set(order)) == len(order)


def test_stage_ordinals_are_contiguous():
    ordinals = sorted(stage.ordinal for stage in STAGES)

    assert ordinals == list(range(1, len(STAGES) + 1))


def test_bronze_silver_gold_are_explicit():
    assert stages_for_layer(LAYER_SNAPSHOT)
    assert stages_for_layer(LAYER_BRONZE)
    assert stages_for_layer(LAYER_SILVER)
    assert stages_for_layer(LAYER_GOLD)

    assert {stage.layer for stage in STAGES} == {
        LAYER_SNAPSHOT,
        LAYER_BRONZE,
        LAYER_SILVER,
        LAYER_GOLD,
    }


def test_silver_is_declared_as_not_persisted():
    silver = stages_for_layer(LAYER_SILVER)

    assert silver
    assert all(stage.persisted is False for stage in silver)


def test_dependencies_only_point_backwards():
    ordinal = {stage.stage_id: stage.ordinal for stage in STAGES}

    for stage in STAGES:
        for dependency in stage.depends_on:
            assert dependency in ordinal, stage.stage_id
            assert ordinal[dependency] < stage.ordinal, stage.stage_id


def test_every_stage_declares_a_validation_gate():
    for stage in STAGES:
        assert stage.validation, stage.stage_id
        assert stage.purpose, stage.stage_id


def test_compatibility_keys_reference_real_configuration_keys():
    resolved = resolve_run_config()

    for stage in STAGES:
        for config_key in stage.compatibility.values():
            if config_key.startswith("__"):
                continue  # derived pseudo-key, resolved by the planner

            assert resolved.get(config_key, "__missing__") != "__missing__", (
                f"{stage.stage_id}: {config_key}"
            )


def test_scientific_parameters_reference_real_configuration_keys():
    resolved = resolve_run_config()

    for stage in STAGES:
        for key in stage.scientific_parameters:
            assert resolved.get(key, "__missing__") != "__missing__", (
                f"{stage.stage_id}: {key}"
            )


# --------------------------------------------------------------------------
# Paths
# --------------------------------------------------------------------------


def test_paths_require_a_pinned_snapshot(tmp_path):
    with pytest.raises(PipelineError, match="snapshot"):
        PipelinePaths.from_config(resolve_run_config(), repo_root=tmp_path)


def test_release_name_carries_snapshot_and_protocol(tmp_path):
    paths = PipelinePaths.from_config(_config(tmp_path), repo_root=tmp_path)

    assert "20260101" in paths.release
    assert paths.protocol in paths.release


# --------------------------------------------------------------------------
# Stage inspection
# --------------------------------------------------------------------------


def test_snapshot_stage_passes_once_pinned(tmp_path):
    resolved = _config(tmp_path)
    paths = PipelinePaths.from_config(resolved, repo_root=tmp_path)

    observation = inspect_stage(
        _spec("snapshot"), paths=paths, resolved=resolved
    )

    assert observation.status == COMPLETE
    assert observation.action == ACTION_REUSE


def test_snapshot_stage_is_pending_until_pinned(tmp_path):
    resolved = resolve_run_config()
    paths = PipelinePaths.from_config(_config(tmp_path), repo_root=tmp_path)

    observation = inspect_stage(
        _spec("snapshot"), paths=paths, resolved=resolved
    )

    assert observation.status == PENDING
    assert observation.action == ACTION_RUN


def test_missing_output_is_pending(tmp_path):
    resolved = _config(tmp_path)
    paths = PipelinePaths.from_config(resolved, repo_root=tmp_path)

    observation = inspect_stage(_spec("brain"), paths=paths, resolved=resolved)

    assert observation.status == PENDING
    assert observation.action == ACTION_RUN


def test_output_without_a_success_marker_is_partial(tmp_path):
    resolved = _config(tmp_path)
    paths = PipelinePaths.from_config(resolved, repo_root=tmp_path)

    _materialise(
        paths,
        "brain",
        {"snapshot": "20260101", "cleaning_protocol": resolved.get(
            "release.protocol_version"
        )},
        success=False,
    )

    observation = inspect_stage(_spec("brain"), paths=paths, resolved=resolved)

    assert observation.status == PARTIAL
    assert observation.action == ACTION_RUN


def test_complete_output_is_reused(tmp_path):
    resolved = _config(tmp_path)
    paths = PipelinePaths.from_config(resolved, repo_root=tmp_path)

    _materialise(
        paths,
        "brain",
        {
            "snapshot": "20260101",
            "cleaning_protocol": resolved.get("release.protocol_version"),
            "brain_chain_count": 577760,
        },
    )

    observation = inspect_stage(_spec("brain"), paths=paths, resolved=resolved)

    assert observation.status == COMPLETE
    assert observation.validation == VALIDATION_PASS
    assert observation.action == ACTION_REUSE
    assert observation.output_count == 577760


def test_output_from_another_snapshot_is_never_reused(tmp_path):
    resolved = _config(tmp_path)
    paths = PipelinePaths.from_config(resolved, repo_root=tmp_path)

    _materialise(
        paths,
        "brain",
        {
            "snapshot": "20250101",
            "cleaning_protocol": resolved.get("release.protocol_version"),
        },
    )

    observation = inspect_stage(_spec("brain"), paths=paths, resolved=resolved)

    assert observation.status == EXECUTION_COMPLETE
    assert observation.validation == VALIDATION_FAIL
    assert observation.action == ACTION_RUN
    assert observation.incompatibilities


def test_reuse_is_decided_from_summaries_not_directory_names(tmp_path):
    """A correctly named directory with a foreign summary is still rejected."""

    resolved = _config(tmp_path)
    paths = PipelinePaths.from_config(resolved, repo_root=tmp_path)

    directory = _materialise(
        paths,
        "brain",
        {"snapshot": "20250101", "cleaning_protocol": "protocol9.9"},
    )

    assert "20260101" in str(directory)

    observation = inspect_stage(_spec("brain"), paths=paths, resolved=resolved)

    assert observation.action == ACTION_RUN


def test_silver_stage_is_not_applicable(tmp_path):
    resolved = _config(tmp_path)
    paths = PipelinePaths.from_config(resolved, repo_root=tmp_path)

    silver = stages_for_layer(LAYER_SILVER)[0]

    observation = inspect_stage(silver, paths=paths, resolved=resolved)

    assert observation.status == NOT_APPLICABLE
    assert observation.action == ACTION_NOT_APPLICABLE


# --------------------------------------------------------------------------
# Planning
# --------------------------------------------------------------------------


def test_empty_tree_plans_every_stage(tmp_path):
    resolved = _config(tmp_path)

    plan = plan_pipeline(resolved, repo_root=tmp_path)

    assert plan.complete is False

    by_id = plan.by_id

    assert by_id["snapshot"].action == ACTION_REUSE
    assert by_id["bronze_source_manifest"].action == ACTION_RUN


def test_downstream_stages_are_blocked_until_upstream_passes(tmp_path):
    resolved = _config(tmp_path)

    plan = plan_pipeline(resolved, repo_root=tmp_path)

    by_id = plan.by_id

    assert by_id["representative_selection"].status == BLOCKED
    assert by_id["representative_selection"].action == ACTION_BLOCKED
    assert plan.blocked


def test_plan_reports_both_hashes(tmp_path):
    resolved = _config(tmp_path)

    payload = plan_pipeline(resolved, repo_root=tmp_path).to_dict()

    assert payload["resolved_config_sha256"] == resolved.sha256
    assert payload["scientific_config_sha256"] == resolved.scientific_sha256
    assert payload["snapshot"] == "20260101"
    assert len(payload["stages"]) == len(STAGES)


# --------------------------------------------------------------------------
# Gold release presentation
# --------------------------------------------------------------------------


def test_no_release_figures_are_shown_before_the_release_passes(tmp_path):
    resolved = _config(tmp_path)

    plan = plan_pipeline(resolved, repo_root=tmp_path)

    assert gold_release_summary(plan) == {}


# --------------------------------------------------------------------------
# Executors
# --------------------------------------------------------------------------


def test_dry_run_is_the_default_executor(tmp_path):
    executor = build_executor("dry-run")

    assert isinstance(executor, DryRunExecutor)

    result = executor.run("brain", ["echo", "hello"], cwd=tmp_path)

    assert result.executed is False
    assert result.returncode is None


def test_unknown_executor_is_refused():
    with pytest.raises(PipelineError):
        build_executor("magic")


def test_local_executor_runs_the_command(tmp_path):
    executor = build_executor("local")

    result = executor.run("brain", ["true"], cwd=tmp_path)

    assert result.returncode == 0


def test_status_vocabulary_is_complete():
    """The UI must be able to distinguish every documented state."""

    vocabulary = {
        pipeline_module.PENDING,
        pipeline_module.BLOCKED,
        pipeline_module.RUNNING,
        pipeline_module.EXECUTION_COMPLETE,
        pipeline_module.VALIDATING,
        pipeline_module.PARTIAL,
        pipeline_module.VALIDATION_FAIL,
        pipeline_module.VALIDATION_PASS,
        pipeline_module.COMPLETE,
        pipeline_module.NOT_APPLICABLE,
    }

    assert len(vocabulary) == 10


# --------------------------------------------------------------------------
# Canonical stage identity
# --------------------------------------------------------------------------
#
# The registry's `ordinal` is an execution order index. The project's
# scientific vocabulary is Stage 1..14, fixed by the frozen findings document
# and by the `pdbclean_stage<N>_*` schema names inside the frozen artefacts.
# These tests pin the mapping so a stage can never be silently renumbered.


CANONICAL_MAPPING = {
    "snapshot": "Prerequisite A",
    "bronze_source_manifest": "Prerequisite B",
    "silver_parse": "Prerequisite C",
    "structural_cleaning": "Stage 1",
    "geometric_validation": "Stage 2",
    "complete_bri": "Stage 3-4",
    "brain": "Stage 5",
    "length_buckets": "Stage 6",
    "candidate_filtering": "Stage 7",
    "complete_bri_nn": "Stage 8-9",
    "duplicate_classification": "Stage 10",
    "downstream_metadata": "Stage 14 input",
    "redundancy_graph": "Stage 14a",
    "representative_selection": "Stage 14b",
    "gold_release": "Stage 14c",
}


def test_canonical_stage_mapping_is_exact():
    from pdbclean.stage_registry import STAGES as ALL

    assert {s.stage_id: s.canonical_stage for s in ALL} == CANONICAL_MAPPING


def test_no_scientific_stage_is_renumbered():
    """Each canonical Stage N maps to exactly the entry that implements it."""

    by_canonical = {}

    for stage in STAGES:
        by_canonical.setdefault(stage.canonical_stage, []).append(
            stage.stage_id
        )

    # Every canonical scientific label is claimed by exactly one entry.
    for label, owners in by_canonical.items():
        if label == "prerequisite":
            continue

        assert len(owners) == 1, f"{label} claimed by {owners}"


def test_bronze_prerequisite_precedes_every_scientific_stage():
    scientific = [
        s for s in STAGES if not s.canonical_stage.startswith("Prerequisite")
    ]
    prerequisites = [
        s for s in STAGES if s.canonical_stage.startswith("Prerequisite")
    ]

    assert max(s.ordinal for s in prerequisites) < min(
        s.ordinal for s in scientific
    )

    # The Bronze source inventory exists before Stage 1.
    bronze = next(s for s in STAGES if s.stage_id == "bronze_source_manifest")
    stage_one = next(
        s for s in STAGES if s.canonical_stage == "Stage 1"
    )

    assert bronze.ordinal < stage_one.ordinal
    assert bronze.layer == LAYER_BRONZE


def test_the_graph_is_stage_14_not_stage_13():
    """Canonical Stage 13 is the manual review subset, not the deletion set."""

    graph = next(s for s in STAGES if s.stage_id == "redundancy_graph")

    assert graph.canonical_stage.startswith("Stage 14")
    assert graph.canonical_stage != "Stage 13"

    # No orchestrated entry claims Stage 11, 12 or 13 at all.
    claimed = {s.canonical_stage for s in STAGES}

    for label in ("Stage 11", "Stage 12", "Stage 13"):
        assert label not in claimed


def test_investigation_stages_are_declared_and_kept_off_the_release_path():
    from pdbclean.stage_registry import NON_ORCHESTRATED_CANONICAL_STAGES

    labels = {
        entry["canonical_stage"] for entry in NON_ORCHESTRATED_CANONICAL_STAGES
    }

    assert labels == {"Stage 11", "Stage 12", "Stage 13"}

    for entry in NON_ORCHESTRATED_CANONICAL_STAGES:
        assert entry["role"] in {"investigation", "validation"}
        assert "deletion relation" in entry["note"]

    # Nothing on the release path depends on them.
    orchestrated = {s.stage_id for s in STAGES}

    for stage in STAGES:
        for dependency in stage.depends_on:
            assert dependency in orchestrated


def test_no_extra_scientific_stage_was_invented():
    """Only Stages 1-10 and 14 are orchestrated; nothing beyond Stage 14."""

    scientific = {
        s.canonical_stage
        for s in STAGES
        if not s.canonical_stage.startswith("Prerequisite")
        and s.canonical_stage != "Stage 14 input"
    }

    assert scientific == {
        "Stage 1", "Stage 2", "Stage 3-4", "Stage 5", "Stage 6",
        "Stage 7", "Stage 8-9", "Stage 10",
        "Stage 14a", "Stage 14b", "Stage 14c",
    }


def test_plan_and_catalogue_report_the_canonical_label(tmp_path):
    from pdbclean.stage_registry import stage_catalogue

    for entry in stage_catalogue():
        assert entry["canonical_stage"] == CANONICAL_MAPPING[entry["stage_id"]]

    plan = plan_pipeline(_config(tmp_path), repo_root=tmp_path).to_dict()

    for entry in plan["stages"]:
        assert entry["canonical_stage"] == CANONICAL_MAPPING[entry["stage_id"]]


# --------------------------------------------------------------------------
# Frozen dataset counts vs generic snapshots
# --------------------------------------------------------------------------
#
# 578524 / 499770 / 78754 / 99854 and friends are the known-correct counts for
# the frozen 2026-01-01 experiment. They are acceptance gates for THAT run, not
# scientific defaults. A new snapshot must derive its own counts and must never
# be required to match 2026.


FROZEN_2026_COUNTS = {
    578524, 577760, 499770, 78754, 99854, 1068256, 4495,
    1072751, 17373, 1055378,
}

REPO = Path(__file__).resolve().parents[2]

FROZEN_PROFILE = (
    REPO / "config" / "pdbclean" / "profiles" / "comp702_frozen_20260101.yaml"
)

GATED_STAGES = (
    "redundancy_graph",
    "representative_selection",
    "gold_release",
)


def _argv(stage_id, resolved):
    from pdbclean.cli import stage_command

    paths = PipelinePaths.from_config(resolved, repo_root=REPO)

    return stage_command(stage_id, resolved, paths)


def test_validated_defaults_carry_no_frozen_counts():
    from pdbclean.defaults import VALIDATED_DEFAULTS

    rendered = json.dumps(VALIDATED_DEFAULTS)

    for count in FROZEN_2026_COUNTS:
        assert str(count) not in rendered, count


def test_generic_snapshot_inherits_no_frozen_count():
    """The whole point: a new snapshot is not measured against 2026."""

    resolved = resolve_run_config(
        overrides=["snapshot.mode=fixed", "snapshot.snapshot_id=20260701"]
    )

    assert all(
        value is None for value in (resolved.get("expectations") or {}).values()
    )

    for stage_id in GATED_STAGES:
        argv = _argv(stage_id, resolved)

        assert argv is not None

        for token in argv:
            assert token not in {str(c) for c in FROZEN_2026_COUNTS}, (
                f"{stage_id} leaked a frozen 2026 count: {token}"
            )


def test_generic_snapshot_acknowledges_the_absence_of_gates():
    """A gate must never be dropped silently -- it is stated explicitly."""

    resolved = resolve_run_config(
        overrides=["snapshot.mode=fixed", "snapshot.snapshot_id=20260701"]
    )

    assert "--no-expectation-gate" in _argv("redundancy_graph", resolved)
    assert "--no-expectation-gate" in _argv("gold_release", resolved)


def test_generic_snapshot_derives_its_population_from_upstream_output():
    """The m >= 2 population comes from this snapshot's own Stage-6 summary."""

    resolved = resolve_run_config(
        overrides=["snapshot.mode=fixed", "snapshot.snapshot_id=20260701"]
    )

    argv = _argv("redundancy_graph", resolved)

    assert "--length-buckets-summary" in argv

    summary = argv[argv.index("--length-buckets-summary") + 1]

    assert summary.endswith("length_buckets/global_summary.json")
    assert "20260701" in summary
    assert "20260101" not in summary

    assert "--expected-mge2-nodes" not in argv


def test_frozen_profile_still_gates_on_the_frozen_counts():
    """Explicit 2026 regression mode keeps every expectation, unchanged."""

    resolved = resolve_run_config(config_path=FROZEN_PROFILE)

    graph = _argv("redundancy_graph", resolved)

    assert graph[graph.index("--expected-edges") + 1] == "1068256"
    assert graph[graph.index("--expected-m1-edges") + 1] == "4495"
    assert graph[graph.index("--expected-mge2-nodes") + 1] == "577760"
    assert "--no-expectation-gate" not in graph

    reps = _argv("representative_selection", resolved)

    assert reps[reps.index("--expected-canonical-input-chains") + 1] == "578524"

    release = _argv("gold_release", resolved)

    assert release[release.index("--expected-retained-chains") + 1] == "499770"
    assert release[release.index("--expected-removed-chains") + 1] == "78754"
    assert "--no-expectation-gate" not in release


def test_expectations_never_change_the_scientific_identity():
    """Gates assert a result; they are not parameters that produce one."""

    frozen = resolve_run_config(config_path=FROZEN_PROFILE)
    ungated = resolve_run_config(
        overrides=["snapshot.mode=fixed", "snapshot.snapshot_id=20260101"]
    )

    assert frozen.scientific_sha256 == ungated.scientific_sha256
    assert frozen.sha256 != ungated.sha256


def test_frozen_and_generic_runs_share_the_same_scientific_parameters():
    """Only the dataset and its gates differ -- never the method."""

    frozen = resolve_run_config(config_path=FROZEN_PROFILE)
    generic = resolve_run_config()

    for key in (
        "duplicate_search.near_duplicate_threshold_angstrom",
        "duplicate_search.operator",
        "duplicate_search.metric",
        "duplicate_search.final_classification_basis",
        "brain_filter.threshold_angstrom",
        "brain.dimension",
        "selection.models.model_id",
        "graph.require_direct_edge_for_removal",
        "representative_selection.minimum_deduplicated_chain_length",
    ):
        assert frozen.get(key) == generic.get(key), key

    assert frozen.near_duplicate_threshold_mA == generic.near_duplicate_threshold_mA == 10
