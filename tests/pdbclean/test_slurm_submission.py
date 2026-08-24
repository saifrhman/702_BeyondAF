"""Submission is safe, idempotent, gated, and fully recorded.

Slurm itself is mocked -- these tests must not put jobs on a real queue -- but
everything above it is the production code path: the same preflight, the same
argv construction, the same ledger, the same state machine.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from pdbclean import pipeline as pipeline_module
from pdbclean.pipeline import plan_pipeline
from pdbclean.run_provenance import FrozenRun, RunProvenance
from pdbclean.runconfig import resolve_run_config
from pdbclean.slurm import (
    JOB_ID_PATTERN,
    SlurmClient,
    SlurmError,
    SlurmJobState,
    SubmissionRequest,
)
from pdbclean.submission import (
    ACTIVE_STATES,
    STATE_BLOCKED,
    STATE_COMPLETE,
    STATE_FAILED,
    STATE_QUEUED,
    STATE_READY,
    STATE_RUNNING,
    STATE_VALIDATING,
    DuplicateSubmission,
    SubmissionError,
    SubmissionLedger,
    Submitter,
    cancel_stage,
    next_submittable_stage,
    run_status,
    stage_status,
    submittable_stage_ids,
)


REPO_ROOT = Path(__file__).resolve().parents[2]

STAGE = "candidate_filtering"


# --------------------------------------------------------------------------
# A Slurm client that records instead of submitting
# --------------------------------------------------------------------------


class FakeSlurm:
    """Records submissions; answers state queries from a script."""

    def __init__(self, *, states=None, fail_on=None):
        self.submissions: list[dict] = []
        self.cancelled: list[str] = []
        self.states = states or {}
        self.fail_on = fail_on or set()
        self._next = 1000

    available = True

    def describe(self):
        return {"available": True, "binaries": {"sbatch": "/usr/bin/sbatch"}}

    def submit(self, request: SubmissionRequest, *, cwd):
        # Build the real argv, so every validation in SubmissionRequest runs.
        argv = request.sbatch_argv("/usr/bin/sbatch")

        if request.job_name in self.fail_on:
            raise SlurmError("sbatch failed with exit code 1: queue is full")

        self._next += 1
        job_id = str(self._next)

        self.submissions.append(
            {"argv": argv, "request": request, "job_id": job_id, "cwd": cwd}
        )

        return {"job_id": job_id, "argv": argv, "stdout": job_id}

    def state(self, job_id):
        return SlurmJobState(
            job_id=job_id, state=self.states.get(job_id, "COMPLETED")
        )

    def cancel(self, job_id):
        self.cancelled.append(job_id)


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


@pytest.fixture()
def plan(frozen):
    return plan_pipeline(frozen.resolved, repo_root=REPO_ROOT)


def _submitter(client, **kwargs):
    return Submitter(
        repo_root=REPO_ROOT,
        client=client,
        allow_dirty_worktree=kwargs.pop("allow_dirty_worktree", True),
        **kwargs,
    )


# --------------------------------------------------------------------------
# Argument safety
# --------------------------------------------------------------------------


def test_submission_is_an_argument_vector_never_a_shell_string():
    request = SubmissionRequest(
        script=REPO_ROOT / "task_scripts/run_stage.sbatch",
        script_args=("/runs/run-1", STAGE),
    )

    argv = request.sbatch_argv("/usr/bin/sbatch")

    assert isinstance(argv, list)
    assert all(isinstance(part, str) for part in argv)
    assert argv[0] == "/usr/bin/sbatch"
    assert "--parsable" in argv


@pytest.mark.parametrize(
    "field,value",
    [
        ("partition", "nodes; rm -rf /"),
        ("time_limit", "$(whoami)"),
        ("memory", "16G && curl evil"),
        ("job_name", "a b"),
        ("array", "0-3%4; echo"),
        ("dependency_afterok", "123 || true"),
    ],
)
def test_hostile_values_are_refused_at_the_boundary(field, value):
    request = SubmissionRequest(
        script=REPO_ROOT / "task_scripts/run_stage.sbatch",
        **{field: value},
    )

    with pytest.raises(SlurmError, match="invalid"):
        request.sbatch_argv("/usr/bin/sbatch")


def test_a_missing_batch_script_is_refused():
    request = SubmissionRequest(script=Path("/no/such/script.sbatch"))

    with pytest.raises(SlurmError, match="does not exist"):
        request.sbatch_argv("/usr/bin/sbatch")


def test_only_registered_stages_are_submittable():
    ids = submittable_stage_ids()

    assert STAGE in ids
    assert "snapshot" not in ids
    assert "silver_parse" not in ids


def test_an_unknown_stage_is_refused(frozen, plan):
    submitter = _submitter(FakeSlurm())

    with pytest.raises(SubmissionError, match="Unknown stage"):
        submitter.submit(frozen, "rm -rf /", plan=plan)


# --------------------------------------------------------------------------
# What gets submitted
# --------------------------------------------------------------------------


def test_submission_carries_the_run_the_commit_and_both_digests(frozen, plan):
    client = FakeSlurm()

    entry = _submitter(client).submit(frozen, STAGE, plan=plan)

    assert entry["run_id"] == frozen.run_id
    assert entry["stage_id"] == STAGE
    assert entry["pipeline_git_commit"] == frozen.git_commit
    assert entry["resolved_config_sha256"] == frozen.resolved_config_sha256
    assert entry["scientific_config_sha256"] == (
        frozen.scientific_config_sha256
    )
    assert entry["stage_config_path"] == str(frozen.stage_config_path)
    assert entry["submitted_at"]

    job = entry["jobs"][0]

    assert JOB_ID_PATTERN.fullmatch(job["job_id"])
    assert job["stdout_path"] and job["stderr_path"]

    argv = client.submissions[0]["argv"]

    # The submitted job's own arguments pin the run it may execute.
    assert str(frozen.run_dir) in argv
    assert STAGE in argv
    assert frozen.git_commit in argv
    assert frozen.resolved_config_sha256 in argv
    assert frozen.scientific_config_sha256 in argv
    assert any(part.endswith("task_scripts/run_stage.sbatch") for part in argv)


def test_resource_request_comes_from_the_registry(frozen, plan):
    client = FakeSlurm()

    _submitter(client).submit(frozen, STAGE, plan=plan)

    request = client.submissions[0]["request"]

    assert request.partition == "nodes"
    assert request.memory == "32G"
    assert request.cpus == 2


def test_resources_are_overridable_without_changing_the_science(tmp_path):
    """Tuning a time limit must not change a run's scientific identity."""

    baseline = resolve_run_config(
        overrides=["snapshot.mode=fixed", "snapshot.snapshot_id=20260101"]
    )

    tuned = resolve_run_config(
        overrides=[
            "snapshot.mode=fixed",
            "snapshot.snapshot_id=20260101",
            f"execution.slurm.{STAGE}.time_limit=48:00:00",
            f"execution.slurm.{STAGE}.memory=256G",
        ]
    )

    from pdbclean.cli import slurm_resources_for

    resources = slurm_resources_for(STAGE, tuned)

    # YAML reads an unquoted 48:00:00 as base 60, i.e. 172800 seconds. That is
    # the same duration, so it is repaired into Slurm's canonical D-HH:MM:SS
    # rather than being passed on as a bare integer sbatch would reject.
    assert resources.time_limit == "2-00:00:00"

    from pdbclean.slurm import TIME_LIMIT_PATTERN

    assert TIME_LIMIT_PATTERN.fullmatch(resources.time_limit)
    assert resources.memory == "256G"
    assert tuned.scientific_sha256 == baseline.scientific_sha256
    assert tuned.sha256 != baseline.sha256


def test_the_ledger_is_append_only(frozen, plan):
    client = FakeSlurm(states={})

    _submitter(client).submit(frozen, STAGE, plan=plan)

    ledger = SubmissionLedger(frozen.run_dir)

    assert len(ledger.records()) == 1
    assert ledger.path.name == "submissions.jsonl"

    lines = ledger.path.read_text(encoding="utf-8").strip().splitlines()

    assert len(lines) == 1
    assert json.loads(lines[0])["stage_id"] == STAGE


# --------------------------------------------------------------------------
# Idempotency and retry
# --------------------------------------------------------------------------


def test_a_second_press_of_run_does_not_submit_a_duplicate(frozen, plan):
    client = FakeSlurm(states={"1001": "RUNNING"})

    submitter = _submitter(client)

    submitter.submit(frozen, STAGE, plan=plan)

    with pytest.raises(DuplicateSubmission, match="already RUNNING"):
        submitter.submit(frozen, STAGE, plan=plan)

    assert len(client.submissions) == 1


def test_a_completed_stage_is_not_resubmitted(frozen, plan, monkeypatch):
    client = FakeSlurm(states={"1001": "COMPLETED"})

    submitter = _submitter(client)
    submitter.submit(frozen, STAGE, plan=plan)

    _pretend_validated(frozen, plan, monkeypatch)

    with pytest.raises(DuplicateSubmission, match="already COMPLETE"):
        submitter.submit(frozen, STAGE, plan=plan)


def test_a_failed_stage_needs_an_explicit_retry(frozen, plan):
    client = FakeSlurm(states={"1001": "FAILED"})

    submitter = _submitter(client)
    submitter.submit(frozen, STAGE, plan=plan)

    with pytest.raises(SubmissionError, match="previously FAILED"):
        submitter.submit(frozen, STAGE, plan=plan)

    client.states["1002"] = "RUNNING"

    entry = submitter.submit(frozen, STAGE, plan=plan, retry=True)

    assert entry["retry"] is True
    assert len(client.submissions) == 2


def test_a_failed_submission_is_recorded_and_raises(frozen, plan):
    client = FakeSlurm(fail_on={f"pdbclean-{STAGE}".replace("_", "-")})

    with pytest.raises(SubmissionError, match="queue is full"):
        _submitter(client).submit(frozen, STAGE, plan=plan)

    record = SubmissionLedger(frozen.run_dir).latest(STAGE)

    assert record["outcome"] == "failed"
    assert "queue is full" in record["error"]
    assert record["jobs"] == []

    status = stage_status(frozen, STAGE, plan=plan, client=client)

    assert status.state == STATE_FAILED
    assert status.may_retry


# --------------------------------------------------------------------------
# The preflight gates
# --------------------------------------------------------------------------


def test_submission_refuses_a_run_frozen_at_another_commit(frozen, plan):
    frozen.record["git"]["commit"] = "f" * 40
    (frozen.run_dir / "run.json").write_text(
        json.dumps(frozen.record), encoding="utf-8"
    )

    reloaded = FrozenRun.load(frozen.run_dir)

    with pytest.raises(SubmissionError, match="frozen at commit"):
        _submitter(FakeSlurm()).submit(reloaded, STAGE, plan=plan)


def test_submission_refuses_a_tampered_configuration(frozen, plan):
    import yaml

    path = frozen.stage_config_path
    document = yaml.safe_load(path.read_text(encoding="utf-8"))
    document["brain_filter"]["threshold_angstrom"] = 0.010
    path.write_text(yaml.safe_dump(document), encoding="utf-8")

    with pytest.raises(SubmissionError, match="not the projection"):
        _submitter(FakeSlurm()).submit(frozen, STAGE, plan=plan)


def test_submission_refuses_a_legacy_run(frozen, plan):
    frozen.stage_config_path.unlink()

    reloaded = FrozenRun.load(frozen.run_dir)

    with pytest.raises(SubmissionError, match="historical record"):
        _submitter(FakeSlurm()).submit(reloaded, STAGE, plan=plan)


def test_submission_refuses_a_dirty_worktree_by_default(frozen, plan):
    """The production rule the pipeline already encodes, applied to submission."""

    submitter = Submitter(
        repo_root=REPO_ROOT,
        client=FakeSlurm(),
        allow_dirty_worktree=False,
    )

    import subprocess

    status = subprocess.run(
        ["git", "-C", str(REPO_ROOT), "status", "--porcelain=v1",
         "--untracked-files=all"],
        capture_output=True, text=True, check=False,
    ).stdout.strip()

    if not status:
        pytest.skip("worktree is clean; nothing to refuse")

    with pytest.raises(SubmissionError, match="not clean"):
        submitter.submit(frozen, STAGE, plan=plan)


# --------------------------------------------------------------------------
# Stage gates
# --------------------------------------------------------------------------


def _pretend_validated(frozen, plan, monkeypatch):
    """Make the plan report this stage as validated, as a real run would."""

    observation = plan.by_id[STAGE]
    observation.validation = pipeline_module.VALIDATION_PASS
    observation.action = pipeline_module.ACTION_RUN

    from pdbclean import submission as submission_module

    record = {
        "preflight_verdict": "PASS",
        "value_agreement": "PASS",
        "returncode": 0,
        "attempt_id": "test",
    }

    monkeypatch.setattr(
        submission_module,
        "stage_status",
        submission_module.stage_status,
    )

    directory = frozen.run_dir / "execution" / STAGE
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "test.json").write_text(json.dumps(record), encoding="utf-8")


def test_slurm_success_alone_is_not_completion(frozen, plan):
    """Exit code zero without validation must stay VALIDATING, never COMPLETE."""

    client = FakeSlurm(states={"1001": "COMPLETED"})

    _submitter(client).submit(frozen, STAGE, plan=plan)

    # No execution record and no validation pass yet.
    status = stage_status(frozen, STAGE, plan=plan, client=client)

    assert status.state == STATE_VALIDATING
    assert "no execution provenance" in status.reason


def test_a_value_mismatch_fails_the_stage_even_after_slurm_success(
    frozen, plan
):
    client = FakeSlurm(states={"1001": "COMPLETED"})

    _submitter(client).submit(frozen, STAGE, plan=plan)

    directory = frozen.run_dir / "execution" / STAGE
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "attempt.json").write_text(
        json.dumps(
            {
                "preflight_verdict": "PASS",
                "value_agreement": "FAIL",
                "value_mismatches": [
                    {
                        "key": "brain_filter_threshold_angstrom",
                        "executed": 0.010,
                        "frozen": 0.005,
                    }
                ],
                "returncode": 0,
            }
        ),
        encoding="utf-8",
    )

    status = stage_status(frozen, STAGE, plan=plan, client=client)

    assert status.state == STATE_FAILED
    assert "differ from the frozen" in status.reason


def test_downstream_is_blocked_while_upstream_runs(frozen, plan):
    client = FakeSlurm(states={"1001": "RUNNING"})

    _submitter(client).submit(frozen, STAGE, plan=plan)

    downstream = stage_status(
        frozen, "complete_bri_nn", plan=plan, client=client
    )

    assert downstream.state == STATE_BLOCKED
    assert not downstream.may_submit


def test_only_one_stage_may_start_at_a_time(frozen, plan):
    client = FakeSlurm(states={"1001": "PENDING"})

    statuses = run_status(frozen, plan=plan, client=client)

    first = next_submittable_stage(statuses)

    assert first is not None
    assert first.stage_id == STAGE

    _submitter(client).submit(frozen, first.stage_id, plan=plan)

    statuses = run_status(frozen, plan=plan, client=client)

    assert any(status.state in ACTIVE_STATES for status in statuses)
    assert next_submittable_stage(statuses) is None


def test_a_failure_stops_the_chain(frozen, plan):
    client = FakeSlurm(states={"1001": "FAILED"})

    _submitter(client).submit(frozen, STAGE, plan=plan)

    statuses = run_status(frozen, plan=plan, client=client)

    assert next_submittable_stage(statuses) is None


def test_reused_output_is_reported_complete(frozen, plan):
    """Stages 1-6 do not depend on tau, so a tau study reuses them."""

    client = FakeSlurm()

    status = stage_status(frozen, "brain", plan=plan, client=client)

    if plan.by_id["brain"].action != pipeline_module.ACTION_REUSE:
        pytest.skip("the frozen 20260101 outputs are not present here")

    assert status.state == STATE_COMPLETE
    assert "already exists" in status.reason


# --------------------------------------------------------------------------
# Cancellation
# --------------------------------------------------------------------------


def test_cancelling_a_stage_cancels_every_job_it_created(frozen, plan):
    client = FakeSlurm(states={"1001": "RUNNING"})

    _submitter(client).submit(frozen, STAGE, plan=plan)

    cancelled = cancel_stage(frozen, STAGE, client=client)

    assert cancelled == client.cancelled
    assert len(cancelled) == 1


def test_cancelling_an_unsubmitted_stage_is_a_no_op(frozen):
    client = FakeSlurm()

    assert cancel_stage(frozen, "gold_release", client=client) == []


# --------------------------------------------------------------------------
# Slurm state mapping
# --------------------------------------------------------------------------


@pytest.mark.parametrize(
    "state,queued,running,succeeded,failed",
    [
        ("PENDING", True, False, False, False),
        ("RUNNING", False, True, False, False),
        ("COMPLETED", False, False, True, False),
        ("FAILED", False, False, False, True),
        ("TIMEOUT", False, False, False, True),
        ("OUT_OF_MEMORY", False, False, False, True),
        ("CANCELLED", False, False, False, True),
    ],
)
def test_slurm_states_are_classified(state, queued, running, succeeded, failed):
    job = SlurmJobState(job_id="1", state=state)

    assert job.queued is queued
    assert job.running is running
    assert job.succeeded is succeeded
    assert job.failed is failed


def test_a_host_without_slurm_reports_so_rather_than_pretending():
    client = SlurmClient(binaries={"sbatch": None})

    assert client.available is False
    assert client.describe()["available"] is False
