"""Run-provenance regression tests.

Provenance is first class: it is written *before* work starts, it records the
resolved configuration rather than the input file, and it is append-only.
Historical provenance is never overwritten.
"""

from __future__ import annotations

import json

import pytest
import yaml

from pdbclean.run_provenance import (
    RunProvenance,
    RunProvenanceError,
    collect_environment,
    collect_git_state,
    file_sha256,
    list_runs,
    new_run_id,
)
from pdbclean.runconfig import config_sha256, resolve_run_config


@pytest.fixture()
def resolved():
    return resolve_run_config(
        overrides=["snapshot.mode=fixed", "snapshot.snapshot_id=20260101"]
    )


@pytest.fixture()
def provenance(tmp_path, resolved):
    return RunProvenance.create(
        resolved=resolved,
        run_root=tmp_path / "runs",
        repo_root=tmp_path,
        snapshot={"snapshot_id": "20260101", "display": "2026-01-01"},
        invocation={"argv": ["pdbclean", "run"]},
    )


# --------------------------------------------------------------------------
# Run identity
# --------------------------------------------------------------------------


def test_run_ids_are_unique_and_time_ordered():
    first = new_run_id(resolved_config_sha256="a" * 64)
    second = new_run_id(resolved_config_sha256="a" * 64)

    assert first != second
    assert first.startswith("run-")
    assert second.startswith("run-")


# --------------------------------------------------------------------------
# Creation
# --------------------------------------------------------------------------


def test_provenance_is_written_before_any_work(provenance):
    record = json.loads(provenance.record_path.read_text(encoding="utf-8"))

    assert record["status"] == "created"
    assert record["stages"] == []
    assert record["run_id"] == provenance.run_id


def test_record_carries_the_resolved_configuration(provenance, resolved):
    record = provenance.record

    assert record["resolved_config"] == resolved.to_dict()
    assert record["resolved_config_sha256"] == resolved.sha256
    assert record["scientific_config_sha256"] == resolved.scientific_sha256
    assert record["defaults_version"] == resolved.get("defaults_version")


def test_resolved_configuration_is_persisted_beside_the_record(provenance):
    record = provenance.record

    written = yaml.safe_load(
        open(record["resolved_config_yaml"], encoding="utf-8").read()
    )

    assert config_sha256(written) == record["resolved_config_sha256"]

    as_json = json.loads(
        open(record["resolved_config_json"], encoding="utf-8").read()
    )

    assert as_json == record["resolved_config"]


def test_record_carries_value_provenance(tmp_path, resolved):
    run = RunProvenance.create(
        resolved=resolved,
        run_root=tmp_path / "runs",
        repo_root=tmp_path,
    )

    sources = run.record["config_value_sources"]

    assert sources["brain.dimension"] == "builtin_default"
    assert sources["snapshot.snapshot_id"].startswith("override:")


def test_record_carries_git_and_environment(provenance):
    record = provenance.record

    assert "branch" in record["git"]
    assert "commit" in record["git"]
    assert "working_tree_dirty" in record["git"]

    environment = record["environment"]

    assert environment["python_version"]
    assert "numpy_version" in environment
    assert "pyarrow_version" in environment


def test_snapshot_identity_is_recorded(provenance):
    assert provenance.record["snapshot"]["snapshot_id"] == "20260101"


def test_existing_run_directory_is_never_overwritten(tmp_path, resolved):
    run = RunProvenance.create(
        resolved=resolved,
        run_root=tmp_path / "runs",
        repo_root=tmp_path,
    )

    with pytest.raises(RunProvenanceError, match="historical provenance"):
        RunProvenance.create(
            resolved=resolved,
            run_root=tmp_path / "runs",
            repo_root=tmp_path,
            run_id=run.run_id,
        )


# --------------------------------------------------------------------------
# Stage and event recording
# --------------------------------------------------------------------------


def test_stage_lifecycle_is_recorded(provenance):
    provenance.register_stage(
        stage_id="complete_bri",
        title="Complete BRI",
        layer="gold",
        scientific_parameters={"implementation_version": "1.2.2"},
    )

    provenance.update_stage(
        "complete_bri",
        status="execution_complete",
        slurm_job_ids=["10284723"],
    )
    provenance.update_stage(
        "complete_bri",
        status="complete",
        validation="validation_pass",
        output_count=578524,
    )
    provenance.flush()

    stage = provenance.stage("complete_bri")

    assert stage.status == "complete"
    assert stage.validation == "validation_pass"
    assert stage.output_count == 578524
    assert stage.slurm_job_ids == ["10284723"]
    assert stage.scientific_parameters["implementation_version"] == "1.2.2"


def test_unknown_stage_field_is_refused(provenance):
    provenance.register_stage(
        stage_id="complete_bri", title="Complete BRI", layer="gold"
    )

    with pytest.raises(RunProvenanceError):
        provenance.update_stage("complete_bri", not_a_field=1)


def test_events_are_append_only(provenance):
    before = provenance.events_path.read_text(encoding="utf-8").splitlines()

    provenance.append_event("stage_started", stage_id="complete_bri")
    provenance.append_event("stage_finished", stage_id="complete_bri")

    after = provenance.events_path.read_text(encoding="utf-8").splitlines()

    assert after[: len(before)] == before
    assert len(after) == len(before) + 2

    for line in after:
        payload = json.loads(line)
        assert "at" in payload
        assert "event" in payload


def test_run_creation_event_names_both_hashes(provenance, resolved):
    events = [
        json.loads(line)
        for line in provenance.events_path.read_text(
            encoding="utf-8"
        ).splitlines()
    ]

    created = [e for e in events if e["event"] == "run_created"]

    assert created
    assert created[0]["resolved_config_sha256"] == resolved.sha256
    assert created[0]["scientific_config_sha256"] == resolved.scientific_sha256


def test_validation_verdicts_are_recorded(provenance):
    provenance.record_validation(
        "gold_release", "PASS", retained_chain_count=499770
    )
    provenance.flush()

    verdict = provenance.record["validation"]["gold_release"]

    assert verdict["verdict"] == "PASS"
    assert verdict["retained_chain_count"] == 499770


def test_release_is_recorded_with_artefact_hashes(tmp_path, provenance):
    artefact = tmp_path / "retained_chains.parquet"
    artefact.write_bytes(b"not really parquet")

    provenance.record_release(
        release_path=str(tmp_path),
        artefacts=[
            {
                "path": str(artefact),
                "bytes": artefact.stat().st_size,
                "sha256": file_sha256(artefact),
            }
        ],
        retained_chain_count=499770,
    )
    provenance.flush()

    release = provenance.record["release"]

    assert release["retained_chain_count"] == 499770
    assert release["artefacts"][0]["sha256"] == file_sha256(artefact)


def test_flush_is_reloadable(provenance):
    provenance.register_stage(
        stage_id="gold_release", title="Gold release", layer="gold"
    )
    provenance.update_stage("gold_release", status="complete")
    provenance.set_status("complete")
    provenance.flush()

    reloaded = RunProvenance.load(provenance.run_dir)

    assert reloaded.run_id == provenance.run_id
    assert reloaded.record["status"] == "complete"
    assert reloaded.stage("gold_release").status == "complete"


# --------------------------------------------------------------------------
# Listing
# --------------------------------------------------------------------------


def test_runs_are_listed_with_both_hashes(tmp_path, resolved):
    run_root = tmp_path / "runs"

    first = RunProvenance.create(
        resolved=resolved, run_root=run_root, repo_root=tmp_path
    )
    second = RunProvenance.create(
        resolved=resolved, run_root=run_root, repo_root=tmp_path
    )

    listed = list_runs(run_root)

    assert {entry["run_id"] for entry in listed} == {
        first.run_id,
        second.run_id,
    }

    for entry in listed:
        assert entry["resolved_config_sha256"] == resolved.sha256
        assert entry["scientific_config_sha256"] == resolved.scientific_sha256


def test_listing_a_missing_root_is_empty(tmp_path):
    assert list_runs(tmp_path / "absent") == []


def test_git_state_is_reported_for_a_non_repository(tmp_path):
    state = collect_git_state(tmp_path)

    assert state["commit"] is None


def test_environment_reports_the_bri_version(tmp_path):
    version_file = tmp_path / "reproducibility" / "bri_version.txt"
    version_file.parent.mkdir(parents=True)
    version_file.write_text("1.2.2\n", encoding="utf-8")

    environment = collect_environment(repo_root=tmp_path)

    assert environment["bri_implementation"] == "1.2.2"


# --------------------------------------------------------------------------
# Canonical configuration vs runtime provenance
# --------------------------------------------------------------------------
#
# A run resolves its configuration ONCE, at creation, and persists it. Host
# facts -- $TMPDIR, hostname, Slurm identifiers -- are execution provenance and
# must never alter the canonical document or either of its hashes.


def test_created_run_config_is_immutable_across_hosts(tmp_path, monkeypatch):
    """Re-reading an existing run on another node yields the same identity."""

    monkeypatch.setenv("TMPDIR", "/tmp/login-node")

    run = RunProvenance.create(
        resolved=resolve_run_config(
            overrides=["snapshot.mode=fixed", "snapshot.snapshot_id=20260101"]
        ),
        run_root=tmp_path / "runs",
        repo_root=tmp_path,
    )

    created_sha = run.record["resolved_config_sha256"]
    created_science = run.record["scientific_config_sha256"]

    # The job lands on a compute node with a completely different scratch.
    monkeypatch.setenv("TMPDIR", "/scratch/node010/999999")
    monkeypatch.setenv("SLURM_JOB_ID", "10285026")

    reloaded = RunProvenance.load(run.run_dir)

    assert reloaded.record["resolved_config_sha256"] == created_sha
    assert reloaded.record["scientific_config_sha256"] == created_science

    # And the canonical document itself is loaded, not re-resolved.
    canonical = reloaded.resolved_config

    assert canonical.sha256 == created_sha
    assert canonical.scientific_sha256 == created_science
    assert canonical.get("storage.temporary_root") == "${TMPDIR}/pdbclean"


def test_loaded_config_matches_the_configuration_that_created_the_run(tmp_path):
    resolved = resolve_run_config(
        overrides=["snapshot.mode=fixed", "snapshot.snapshot_id=20260101"]
    )

    run = RunProvenance.create(
        resolved=resolved,
        run_root=tmp_path / "runs",
        repo_root=tmp_path,
    )

    assert run.resolved_config.to_dict() == resolved.to_dict()
    assert run.resolved_config.sha256 == resolved.sha256


def test_loading_a_run_without_a_canonical_config_is_an_error(tmp_path):
    from pdbclean.runconfig import RunConfigError, load_resolved_config

    with pytest.raises(RunConfigError):
        load_resolved_config(tmp_path)


def test_runtime_facts_are_recorded_separately(tmp_path, monkeypatch):
    monkeypatch.setenv("TMPDIR", "/scratch/node010/999999")
    monkeypatch.setenv("SLURM_JOB_ID", "10285026")

    resolved = resolve_run_config(
        overrides=["snapshot.mode=fixed", "snapshot.snapshot_id=20260101"]
    )

    run = RunProvenance.create(
        resolved=resolved,
        run_root=tmp_path / "runs",
        repo_root=tmp_path,
    )

    before = run.record["resolved_config_sha256"]

    observation = run.record_runtime(resolved, stage_id="complete_bri")
    run.flush()

    # The host facts are captured ...
    assert observation["storage"]["temporary_root"] == (
        "/scratch/node010/999999/pdbclean"
    )
    assert observation["environment"]["TMPDIR"] == "/scratch/node010/999999"
    assert observation["environment"]["SLURM_JOB_ID"] == "10285026"
    assert observation["stage_id"] == "complete_bri"
    assert observation["hostname"]

    # ... in the runtime block, not the configuration.
    runtime = run.record["runtime"]

    assert runtime["created_on"]["environment"]["TMPDIR"] == (
        "/scratch/node010/999999"
    )
    assert len(runtime["observed"]) == 1

    assert run.record["resolved_config_sha256"] == before
    assert run.record["resolved_config"]["storage"]["temporary_root"] == (
        "${TMPDIR}/pdbclean"
    )


def test_runtime_observations_accumulate_and_are_logged(tmp_path):
    resolved = resolve_run_config(
        overrides=["snapshot.mode=fixed", "snapshot.snapshot_id=20260101"]
    )

    run = RunProvenance.create(
        resolved=resolved,
        run_root=tmp_path / "runs",
        repo_root=tmp_path,
    )

    run.record_runtime(resolved, stage_id="complete_bri")
    run.record_runtime(resolved, stage_id="brain")
    run.flush()

    assert [o["stage_id"] for o in run.record["runtime"]["observed"]] == [
        "complete_bri",
        "brain",
    ]

    events = [
        json.loads(line)
        for line in run.events_path.read_text(encoding="utf-8").splitlines()
    ]

    assert len([e for e in events if e["event"] == "runtime_environment"]) == 2


def test_provenance_preserves_canonical_stage_identities(tmp_path, resolved):
    """The historical record carries Stage 1-14, not an ordering index."""

    from pdbclean.pipeline import plan_pipeline, record_plan_in_provenance

    run = RunProvenance.create(
        resolved=resolved,
        run_root=tmp_path / "runs",
        repo_root=tmp_path,
    )

    plan = plan_pipeline(
        resolve_run_config(
            overrides=[
                "snapshot.mode=fixed",
                "snapshot.snapshot_id=20260101",
                f"storage.output_root={tmp_path / 'outputs'}",
                f"storage.release_root={tmp_path / 'releases'}",
            ]
        ),
        repo_root=tmp_path,
    )

    record_plan_in_provenance(plan, run)
    run.flush()

    recorded = {
        stage["stage_id"]: stage["canonical_stage"]
        for stage in run.record["stages"]
    }

    assert recorded["structural_cleaning"] == "Stage 1"
    assert recorded["duplicate_classification"] == "Stage 10"
    assert recorded["redundancy_graph"] == "Stage 14a"
    assert recorded["representative_selection"] == "Stage 14b"
    assert recorded["gold_release"] == "Stage 14c"
    assert recorded["bronze_source_manifest"] == "Prerequisite B"
    assert recorded["snapshot"] == "Prerequisite A"

    # Canonical identities survive a reload.
    reloaded = RunProvenance.load(run.run_dir)

    assert reloaded.stage("redundancy_graph").canonical_stage == "Stage 14a"
