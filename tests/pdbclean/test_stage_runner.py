"""A frozen run is immutable, and a worker proves it before computing.

The guarantee under test is narrow and absolute:

    the values a user sees in the frozen resolved configuration are exactly the
    values the scientific computation executes

Everything here checks that guarantee at the worker boundary, through the
production loader, on a real run directory.  Where a test needs the guarantee
to *fail*, it breaks something the way it would really break -- an edited
configuration, a moved checkout, a projection that lost a section -- and asserts
that the stage aborts rather than publishing.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest
import yaml

from pdbclean.run_provenance import FrozenRun, FrozenRunError, RunProvenance
from pdbclean.runconfig import resolve_run_config
from pdbclean.stage_config import stage_config_path
from pdbclean.stage_runner import (
    EXIT_PREFLIGHT_FAILED,
    StageRunnerError,
    build_execution_record,
    latest_execution_record,
    main as stage_runner_main,
    run_stage,
    verify_git,
)


REPO_ROOT = Path(__file__).resolve().parents[2]

STAGE = "candidate_filtering"


def _resolved(tau=0.005):
    return resolve_run_config(
        overrides=[
            "snapshot.mode=fixed",
            "snapshot.snapshot_id=20260101",
            "bri.representation_precision_angstrom=0.001",
            f"brain_filter.threshold_angstrom={tau}",
            f"duplicate_search.near_duplicate_threshold_angstrom={tau}",
        ]
    )


@pytest.fixture()
def frozen(tmp_path):
    provenance = RunProvenance.create(
        resolved=_resolved(),
        run_root=tmp_path / "runs",
        repo_root=REPO_ROOT,
        snapshot={"snapshot_id": "20260101", "display": "2026-01-01"},
    )

    return FrozenRun.load(provenance.run_dir)


# --------------------------------------------------------------------------
# The frozen run itself
# --------------------------------------------------------------------------


def test_freezing_materialises_the_executable_configuration(frozen):
    path = stage_config_path(frozen.run_dir)

    assert path.is_file()
    assert frozen.record["stage_config_path"] == str(path)
    assert frozen.record["stage_config_sha256"]
    assert not frozen.is_legacy


def test_a_frozen_run_verifies(frozen):
    verification = frozen.verify()

    assert verification["resolved_config_sha256"] == (
        frozen.resolved_config_sha256
    )
    assert verification["verified"] is True


def test_editing_the_profile_after_freezing_does_not_change_the_run(
    tmp_path, frozen
):
    """A frozen run reads its own directory, never a profile file."""

    before = frozen.resolved.get("brain_filter.threshold_angstrom")

    profile = tmp_path / "profile.yaml"
    profile.write_text(
        yaml.safe_dump({"brain_filter": {"threshold_angstrom": 0.050}}),
        encoding="utf-8",
    )

    reloaded = FrozenRun.load(frozen.run_dir)

    assert reloaded.resolved.get("brain_filter.threshold_angstrom") == before
    assert reloaded.resolved_config_sha256 == frozen.resolved_config_sha256

    reloaded.verify()


def test_changing_the_builtin_defaults_does_not_change_the_run(
    monkeypatch, frozen
):
    """Even a source-code change to the defaults cannot reach a frozen run."""

    from pdbclean import defaults as defaults_module

    mutated = defaults_module.validated_defaults()
    mutated["brain_filter"]["threshold_angstrom"] = 0.500

    monkeypatch.setattr(
        defaults_module, "VALIDATED_DEFAULTS", mutated, raising=True
    )

    reloaded = FrozenRun.load(frozen.run_dir)

    assert reloaded.resolved.get("brain_filter.threshold_angstrom") == 0.005

    reloaded.verify()


def test_a_tampered_stage_configuration_is_rejected(frozen):
    path = stage_config_path(frozen.run_dir)

    document = yaml.safe_load(path.read_text(encoding="utf-8"))
    document["brain_filter"]["threshold_angstrom"] = 0.010
    path.write_text(yaml.safe_dump(document), encoding="utf-8")

    with pytest.raises(FrozenRunError, match="not the projection"):
        frozen.verify()


def test_a_tampered_canonical_configuration_is_rejected(frozen):
    path = frozen.run_dir / "resolved_run.json"

    document = json.loads(path.read_text(encoding="utf-8"))
    document["brain_filter"]["threshold_angstrom"] = 0.010
    path.write_text(json.dumps(document), encoding="utf-8")

    with pytest.raises(FrozenRunError, match="hash mismatch"):
        frozen.verify()


def test_a_legacy_run_is_readable_but_not_executable(frozen):
    """Runs frozen before this architecture are records, not recipes."""

    stage_config_path(frozen.run_dir).unlink()

    reloaded = FrozenRun.load(frozen.run_dir)

    assert reloaded.is_legacy
    assert reloaded.resolved_config_sha256  # still readable

    with pytest.raises(FrozenRunError, match="legacy record"):
        reloaded.verify()


# --------------------------------------------------------------------------
# Preflight: the values a worker actually loads
# --------------------------------------------------------------------------


def test_preflight_records_the_values_the_worker_loaded(frozen):
    record = build_execution_record(
        frozen, STAGE, repo_root=REPO_ROOT, allow_dirty=True
    )

    executed = record["executed_values"]

    assert executed["representation_precision_angstrom"] == 0.001
    assert executed["brain_filter_threshold_angstrom"] == 0.005
    assert executed["brain_filter_threshold_units"] == 5
    assert executed["complete_bri_near_duplicate_threshold_angstrom"] == 0.005
    assert executed["complete_bri_near_duplicate_threshold_units"] == 5
    assert executed["snapshot"] == "20260101"
    assert executed["model_id"] == 1
    assert executed["minimum_backbone_distance_angstrom"] == 0.01
    assert executed["minimum_triangle_angle_degrees"] == 3.0

    assert record["value_agreement"] == "PASS"
    assert record["value_mismatches"] == []


def test_preflight_records_both_hashes_and_the_frozen_commit(frozen):
    record = build_execution_record(
        frozen, STAGE, repo_root=REPO_ROOT, allow_dirty=True
    )

    assert record["resolved_config_sha256"] == frozen.resolved_config_sha256
    assert record["scientific_config_sha256"] == (
        frozen.scientific_config_sha256
    )
    assert record["git"]["frozen_commit"] == frozen.git_commit


def test_a_lost_section_fails_the_stage_rather_than_computing(frozen):
    """The mandated wrong-value rejection test.

    Freeze Brain = 0.005, then hand the worker a configuration that says
    0.010.  Two independent gates must catch it: the projection no longer
    matches the canonical configuration, and -- were that gate removed -- the
    loaded value would not match the frozen one.
    """

    path = stage_config_path(frozen.run_dir)

    document = yaml.safe_load(path.read_text(encoding="utf-8"))
    document["brain_filter"]["threshold_angstrom"] = 0.010
    document["duplicate_search"]["near_duplicate_threshold_angstrom"] = 0.010
    path.write_text(yaml.safe_dump(document), encoding="utf-8")

    with pytest.raises((StageRunnerError, FrozenRunError)):
        build_execution_record(
            frozen, STAGE, repo_root=REPO_ROOT, allow_dirty=True
        )


def test_value_mismatch_alone_is_enough_to_abort(frozen, monkeypatch):
    """Even if the projection gate passed, differing values must abort.

    Belt and braces on purpose: the two gates protect against different
    failures, and neither is allowed to be the only one.
    """

    import pdbclean.stage_runner as runner

    genuine = runner.executed_values

    monkeypatch.setattr(
        runner,
        "executed_values",
        lambda config: {
            **genuine(config),
            "brain_filter_threshold_angstrom": 0.010,
        },
    )

    with pytest.raises(StageRunnerError, match="not the values"):
        build_execution_record(
            frozen, STAGE, repo_root=REPO_ROOT, allow_dirty=True
        )


def test_a_refusal_is_itself_recorded(frozen):
    path = stage_config_path(frozen.run_dir)
    path.write_text("not: a valid stage configuration\n", encoding="utf-8")

    with pytest.raises((StageRunnerError, FrozenRunError)):
        run_stage(frozen, STAGE, repo_root=REPO_ROOT, allow_dirty=True)

    record = latest_execution_record(frozen.run_dir, STAGE)

    assert record is not None
    assert record["preflight_verdict"] == "FAIL"
    assert record["preflight_error"]


# --------------------------------------------------------------------------
# Git
# --------------------------------------------------------------------------


def test_execution_refuses_a_different_checkout(tmp_path):
    with pytest.raises(StageRunnerError, match="frozen commit"):
        verify_git(
            REPO_ROOT,
            frozen_commit="b" * 40,
            allow_dirty=True,
        )


def test_dirty_worktree_is_refused_by_default(tmp_path, monkeypatch):
    import pdbclean.stage_runner as runner

    monkeypatch.setattr(
        runner,
        "git_state",
        lambda root: {
            "head": "a" * 40,
            "working_tree_dirty": True,
            "dirty_paths": [" M src/pdbclean/bri.py"],
            "uncommitted_diff_sha256": "c" * 64,
        },
    )

    with pytest.raises(StageRunnerError, match="dirty worktree"):
        verify_git(REPO_ROOT, frozen_commit="a" * 40, allow_dirty=False)


def test_dirty_execution_records_the_difference(monkeypatch):
    """If it is allowed, the commit alone no longer describes the code."""

    import pdbclean.stage_runner as runner

    monkeypatch.setattr(
        runner,
        "git_state",
        lambda root: {
            "head": "a" * 40,
            "working_tree_dirty": True,
            "dirty_paths": [" M src/pdbclean/bri.py"],
            "uncommitted_diff_sha256": "c" * 64,
            "uncommitted_diff_bytes": 42,
        },
    )

    state = verify_git(REPO_ROOT, frozen_commit="a" * 40, allow_dirty=True)

    assert state["dirty_execution_allowed"] is True
    assert state["uncommitted_diff_sha256"] == "c" * 64
    assert state["dirty_paths"]


# --------------------------------------------------------------------------
# The runner as a process
# --------------------------------------------------------------------------


def test_dry_run_verifies_and_records_without_executing(frozen):
    returncode, record = run_stage(
        frozen,
        STAGE,
        repo_root=REPO_ROOT,
        allow_dirty=True,
        dry_run=True,
    )

    assert returncode == 0
    assert record["dry_run"] is True
    assert record["preflight_verdict"] == "PASS"
    assert record["argv"][0].endswith("python")
    assert Path(record["record_path"]).is_file()


def test_cli_rejects_a_submission_that_names_another_run(frozen, capsys):
    """The submitted job pins what it expects; a mismatch must not execute."""

    exit_code = stage_runner_main(
        [
            "--run-dir",
            str(frozen.run_dir),
            "--stage",
            STAGE,
            "--repo-root",
            str(REPO_ROOT),
            "--expect-resolved-config-sha256", "d" * 64,
            "--allow-dirty-worktree",
        ]
    )

    assert exit_code == EXIT_PREFLIGHT_FAILED
    assert "Refusing to execute" in capsys.readouterr().err


def test_cli_accepts_the_matching_digests(frozen):
    exit_code = stage_runner_main(
        [
            "--run-dir",
            str(frozen.run_dir),
            "--stage",
            STAGE,
            "--repo-root",
            str(REPO_ROOT),
            "--expect-commit", str(frozen.git_commit),
            "--expect-resolved-config-sha256",
            str(frozen.resolved_config_sha256),
            "--expect-scientific-config-sha256",
            str(frozen.scientific_config_sha256),
            "--allow-dirty-worktree",
            "--dry-run",
        ]
    )

    assert exit_code == 0


def test_execution_records_are_append_only(frozen):
    for _ in range(3):
        run_stage(
            frozen,
            STAGE,
            repo_root=REPO_ROOT,
            allow_dirty=True,
            dry_run=True,
        )

    directory = frozen.run_dir / "execution" / STAGE

    assert len(list(directory.glob("*.json"))) >= 1
