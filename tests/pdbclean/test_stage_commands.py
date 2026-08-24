"""Generated stage commands must actually run.

The defect these tests exist for::

    $ pdbclean stage-command --stage candidate_filtering
    python -m pdbclean.brain_prefilter_production --config .../protocol_3_2...yaml

    $ <that command>
    brain_prefilter_production.py: error: the following arguments are
    required: --pipeline-git-commit

Two things were wrong: a required argument was missing, and the configuration
named was the byte-frozen base file rather than the run's own.

The required-argument set below is **re-derived from each entry point's source
code**, not restated.  Adding a new ``required=True`` option to a stage makes
these tests fail until command generation supplies it, so the registry cannot
drift away from the code.
"""

from __future__ import annotations

import argparse
import ast
from pathlib import Path

import pytest

from pdbclean.cli import (
    preview_stage_command,
    stage_command,
    stage_invocation,
)
from pdbclean.pipeline import PipelineError, PipelinePaths
from pdbclean.run_provenance import RunProvenance, FrozenRun
from pdbclean.runconfig import resolve_run_config
from pdbclean.stage_registry import (
    EXECUTION_NONE,
    STAGES,
    STAGES_BY_ID,
    stage_execution,
)


REPO_ROOT = Path(__file__).resolve().parents[2]

FROZEN_BASE_CONFIG = "config/pdbclean/protocol_3_2_comp702_v1.yaml"

#: Where each orchestrated stage's required arguments are declared.
ENTRY_POINT_SOURCES = {
    "bronze_source_manifest": "scripts/pdbclean/create_manifest.py",
    "geometric_validation": (
        "scripts/pdbclean/finalize_geometric_validation.py"
    ),
    "complete_bri": "scripts/pdbclean/finalize_bri.py",
    "brain": "src/pdbclean/brain_finalize_cli.py",
    "length_buckets": "src/pdbclean/length_buckets_cli.py",
    "candidate_filtering": "src/pdbclean/brain_prefilter_production.py",
    "complete_bri_nn": "src/pdbclean/full_bri_nn_production.py",
    "duplicate_classification": (
        "src/pdbclean/duplicate_classification_production.py"
    ),
    "downstream_metadata": "src/pdbclean/downstream_metadata_finalize.py",
    "redundancy_graph": "scripts/build_stage14_geometric_graph.py",
    "representative_selection": "scripts/select_stage14_representatives.py",
    "gold_release": "scripts/build_stage14_final_release.py",
}


def declared_arguments(source_path: str) -> list[dict]:
    """Read a stage entry point's own ``add_argument`` calls.

    Static, so importing a stage module (which can be expensive, and in one
    case executes work at import) is never required.
    """

    tree = ast.parse((REPO_ROOT / source_path).read_text(encoding="utf-8"))

    arguments: list[dict] = []

    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue

        function = node.func

        if not (
            isinstance(function, ast.Attribute)
            and function.attr == "add_argument"
        ):
            continue

        options = [
            argument.value
            for argument in node.args
            if isinstance(argument, ast.Constant)
            and isinstance(argument.value, str)
            and argument.value.startswith("-")
        ]

        if not options:
            continue

        entry: dict = {"options": options, "required": False, "action": None}

        for keyword in node.keywords:
            if keyword.arg == "required" and isinstance(
                keyword.value, ast.Constant
            ):
                entry["required"] = bool(keyword.value.value)

            if keyword.arg == "action" and isinstance(
                keyword.value, ast.Constant
            ):
                entry["action"] = keyword.value.value

        arguments.append(entry)

    return arguments


def mirror_parser(source_path: str) -> argparse.ArgumentParser:
    """Rebuild the entry point's parser, so its argv can really be parsed."""

    parser = argparse.ArgumentParser(add_help=False)

    for entry in declared_arguments(source_path):
        action = entry["action"] or "store"

        if action in {"store_true", "store_false", "count"}:
            parser.add_argument(*entry["options"], action=action)
        else:
            parser.add_argument(
                *entry["options"],
                required=entry["required"],
                action=action,
            )

    return parser


ORCHESTRATED = sorted(
    stage.stage_id
    for stage in STAGES
    if stage_execution(stage.stage_id).mode != EXECUTION_NONE
    and stage.stage_id in ENTRY_POINT_SOURCES
)


@pytest.fixture()
def resolved():
    return resolve_run_config(
        overrides=["snapshot.mode=fixed", "snapshot.snapshot_id=20260101"]
    )


@pytest.fixture()
def paths(resolved):
    return PipelinePaths.from_config(resolved, repo_root=REPO_ROOT)


def _argv(stage_id, resolved, paths, **kwargs):
    return stage_command(
        stage_id,
        resolved,
        paths,
        config_path=kwargs.pop("config_path", "/frozen/stage_config.yaml"),
        pipeline_git_commit=kwargs.pop("pipeline_git_commit", "a" * 40),
        **kwargs,
    )


# --------------------------------------------------------------------------
# Every declared requirement is supplied
# --------------------------------------------------------------------------


@pytest.mark.parametrize("stage_id", ORCHESTRATED)
def test_generated_command_supplies_every_required_argument(
    stage_id, resolved, paths
):
    argv = _argv(stage_id, resolved, paths)

    assert argv is not None, f"{stage_id} generates no command"

    required = [
        entry["options"][0]
        for entry in declared_arguments(ENTRY_POINT_SOURCES[stage_id])
        if entry["required"]
    ]

    assert required, f"{stage_id} declares no required arguments; check the map"

    missing = [option for option in required if option not in argv]

    assert not missing, (
        f"{stage_id} would abort in argparse: missing {missing}"
    )


@pytest.mark.parametrize("stage_id", ORCHESTRATED)
def test_generated_command_parses_under_the_entry_points_own_parser(
    stage_id, resolved, paths
):
    """Not just present: actually parseable, with no unknown options."""

    argv = _argv(stage_id, resolved, paths)

    parser = mirror_parser(ENTRY_POINT_SOURCES[stage_id])

    # Drop the interpreter / -m / module or script path prefix.
    if argv[1] == "-m":
        arguments = argv[3:]
    else:
        arguments = argv[2:]

    namespace, unknown = parser.parse_known_args(arguments)

    assert not unknown, f"{stage_id} passes unknown options: {unknown}"
    assert namespace is not None


@pytest.mark.parametrize("stage_id", ORCHESTRATED)
def test_registry_requirements_match_the_entry_point_source(stage_id):
    """The registry's declared requirements are the code's, not a guess."""

    declared = {
        entry["options"][0]
        for entry in declared_arguments(ENTRY_POINT_SOURCES[stage_id])
        if entry["required"]
    }

    registered = set(stage_execution(stage_id).required_arguments)

    assert declared <= registered, (
        f"{stage_id}: the entry point requires {sorted(declared - registered)} "
        "but the registry does not list it"
    )


def test_missing_required_argument_is_refused_not_emitted(monkeypatch):
    """Command generation fails loudly rather than emitting a broken argv."""

    from pdbclean import cli as cli_module
    from pdbclean.stage_registry import STAGE_EXECUTION
    import dataclasses

    original = STAGE_EXECUTION["candidate_filtering"]

    monkeypatch.setitem(
        STAGE_EXECUTION,
        "candidate_filtering",
        dataclasses.replace(
            original,
            required_arguments=original.required_arguments + ("--invented",),
        ),
    )

    resolved = resolve_run_config(
        overrides=["snapshot.mode=fixed", "snapshot.snapshot_id=20260101"]
    )
    paths = PipelinePaths.from_config(resolved, repo_root=REPO_ROOT)

    with pytest.raises(PipelineError, match="--invented"):
        _argv("candidate_filtering", resolved, paths)


# --------------------------------------------------------------------------
# The configuration a command names
# --------------------------------------------------------------------------


@pytest.mark.parametrize("stage_id", ORCHESTRATED + ["structural_cleaning"])
def test_no_command_points_at_the_frozen_base_configuration(
    stage_id, resolved, paths
):
    """The base YAML's thresholds are not this run's thresholds."""

    argv = _argv(stage_id, resolved, paths)

    assert argv is not None

    for token in argv:
        assert not token.endswith(FROZEN_BASE_CONFIG.split("/")[-1]), (
            f"{stage_id} still drives a stage from {FROZEN_BASE_CONFIG}"
        )


def test_preview_command_names_a_projection_of_the_entered_values(tmp_path):
    """`pdbclean stage-command` prints something runnable and correct."""

    from pdbclean.config import load_config
    from pdbclean.stage_config import executed_values

    resolved = resolve_run_config(
        overrides=[
            "snapshot.mode=fixed",
            "snapshot.snapshot_id=20260101",
            "brain_filter.threshold_angstrom=0.005",
            "duplicate_search.near_duplicate_threshold_angstrom=0.005",
        ]
    )

    paths = PipelinePaths.from_config(resolved, repo_root=tmp_path)

    argv = preview_stage_command("candidate_filtering", resolved, paths)

    named = argv[argv.index("--config") + 1]

    values = executed_values(load_config(named).data)

    assert values["brain_filter_threshold_angstrom"] == 0.005
    assert values["complete_bri_near_duplicate_threshold_angstrom"] == 0.005


def test_representative_policy_comes_from_the_resolved_configuration(
    resolved, paths, tmp_path
):
    """A run that selects another policy must run that policy.

    The path used to come from an environment variable, so
    ``representative_selection.policy_config`` -- a resolved scientific value
    whose SHA256 is recorded in the release manifest -- was ignored.
    """

    policy = tmp_path / "other_policy.yaml"
    policy.write_text("policy_version: '2.0'\n", encoding="utf-8")

    chosen = resolve_run_config(
        overrides=[
            "snapshot.mode=fixed",
            "snapshot.snapshot_id=20260101",
            f"representative_selection.policy_config={policy}",
        ]
    )

    argv = _argv("representative_selection", chosen, paths)

    assert argv[argv.index("--config") + 1] == str(policy)


# --------------------------------------------------------------------------
# What a frozen run executes
# --------------------------------------------------------------------------


@pytest.fixture()
def frozen(tmp_path):
    resolved = resolve_run_config(
        overrides=[
            "snapshot.mode=fixed",
            "snapshot.snapshot_id=20260101",
            "brain_filter.threshold_angstrom=0.005",
            "duplicate_search.near_duplicate_threshold_angstrom=0.005",
        ]
    )

    provenance = RunProvenance.create(
        resolved=resolved,
        run_root=tmp_path / "runs",
        repo_root=REPO_ROOT,
        snapshot={"snapshot_id": "20260101", "display": "2026-01-01"},
    )

    return FrozenRun.load(provenance.run_dir)


#: Stages the worker drives from a Protocol 3.2 configuration file. For these
#: the file itself carries the scientific values, so it must be the run's own.
PROTOCOL_CONFIG_STAGES = (
    "bronze_source_manifest",
    "structural_cleaning",
    "geometric_validation",
    "complete_bri",
    "brain",
    "length_buckets",
    "candidate_filtering",
    "complete_bri_nn",
    "duplicate_classification",
    "downstream_metadata",
)

#: The Stage-14 scripts take their scientific values as explicit flags rather
#: than from a configuration file. Their `--config`, where present, names the
#: representative policy, which is a separate scientific input.
EXPLICIT_ARGUMENT_STAGES = (
    "redundancy_graph",
    "representative_selection",
    "gold_release",
)


@pytest.mark.parametrize("stage_id", PROTOCOL_CONFIG_STAGES)
def test_frozen_run_commands_name_the_runs_own_configuration(frozen, stage_id):
    invocation = stage_invocation(
        stage_id, frozen, repo_root=REPO_ROOT, verify=False
    )

    assert invocation is not None
    assert invocation.config_path == str(frozen.stage_config_path)
    assert str(frozen.run_dir) in invocation.config_path
    assert invocation.argv.count(invocation.config_path) >= 1


@pytest.mark.parametrize("stage_id", EXPLICIT_ARGUMENT_STAGES)
def test_stage14_commands_carry_the_frozen_threshold_explicitly(
    frozen, stage_id
):
    """No configuration file, so every scientific value is an argument.

    The frozen run resolves tau to 0.005 A, i.e. 5 representation units.  If
    these commands were still built from a re-resolved configuration -- or
    from the script's own default -- this would read 10.
    """

    invocation = stage_invocation(
        stage_id, frozen, repo_root=REPO_ROOT, verify=False
    )

    argv = list(invocation.argv)

    assert "--threshold-mA" in argv
    assert argv[argv.index("--threshold-mA") + 1] == "5"


@pytest.mark.parametrize("stage_id", ORCHESTRATED)
def test_frozen_run_commands_carry_the_frozen_commit(frozen, stage_id):
    invocation = stage_invocation(
        stage_id, frozen, repo_root=REPO_ROOT, verify=False
    )

    assert invocation.pipeline_git_commit == frozen.git_commit

    commit_options = {
        "--pipeline-git-commit",
        "--brain-pipeline-git-commit",
        "--finalizer-pipeline-git-commit",
        "--length-bucket-pipeline-git-commit",
        "--producer-git-commit",
        "--finalizer-git-commit",
    }

    argv = list(invocation.argv)

    for index, token in enumerate(argv):
        if token in commit_options:
            assert argv[index + 1] == frozen.git_commit


def test_manual_reproduction_matches_the_submitted_invocation(frozen):
    """Advanced -> Reproduce manually must describe the same run.

    The orchestrator submits ``run_stage.sbatch <run-dir> <stage>``; the batch
    job then asks for exactly this invocation. Showing the operator anything
    else would be showing them a different experiment.
    """

    from pdbclean.stage_runner import run_stage

    displayed = stage_invocation(
        "candidate_filtering", frozen, repo_root=REPO_ROOT, verify=False
    )

    executed, record = run_stage(
        frozen,
        "candidate_filtering",
        repo_root=REPO_ROOT,
        allow_dirty=True,
        dry_run=True,
    )

    assert executed == 0
    assert record["argv"] == list(displayed.argv)
    assert record["command_text"] == displayed.command_text


def test_every_registered_stage_has_an_execution_description():
    for stage in STAGES:
        assert stage_execution(stage.stage_id).stage_id == stage.stage_id


def test_unknown_stage_has_no_execution_description():
    with pytest.raises(KeyError, match="no execution description"):
        stage_execution("not_a_stage")


def test_non_orchestrated_stages_generate_no_command(resolved, paths):
    for stage_id in ("snapshot", "silver_parse"):
        assert stage_id in STAGES_BY_ID
        assert _argv(stage_id, resolved, paths) is None
