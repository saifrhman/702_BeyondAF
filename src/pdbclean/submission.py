"""Submitting a frozen run's stages to Slurm, and reporting what they did.

The orchestrator's job is narrow and deliberately so:

    freeze -> verify -> sbatch -> record job id -> return

It never computes science, never waits for a job, and never runs anything long
inside the web server.  What it does own is the part that is easy to get wrong:

*the gate*
    a stage may only be submitted when every dependency has reached
    ``COMPLETE``, which means Slurm success **and** the stage's own validation
    **and** an execution record whose loaded values match the frozen ones.
    A process exiting zero is not scientific validation.

*idempotency*
    pressing Run twice must not create two identical jobs.  A stage that is
    SUBMITTING, QUEUED, RUNNING, VALIDATING or COMPLETE refuses resubmission
    unless it failed and the caller explicitly retries.

*provenance*
    every submission appends one immutable line to ``submissions.jsonl`` in the
    run directory: the argv, the job id, both configuration hashes, the frozen
    commit, the log paths, and -- for arrays -- the parent job and its tasks.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Sequence

from pdbclean import pipeline as pipeline_module
from pdbclean.pipeline import PipelinePaths, PipelinePlan
from pdbclean.provenance import ProvenanceError, resolve_clean_git_commit
from pdbclean.run_provenance import FrozenRun, FrozenRunError, utc_now_iso
from pdbclean.slurm import SlurmClient, SlurmError, SubmissionRequest
from pdbclean.stage_registry import (
    EXECUTION_ARRAY_THEN_FINALIZE,
    EXECUTION_NONE,
    EXECUTION_SINGLE,
    EXECUTION_SUBMITTER,
    STAGES_BY_ID,
    stage_execution,
)


#: The generic batch wrapper. It takes a run directory and a stage id and
#: nothing else, so a submission's entire variable surface is two validated
#: tokens plus three digests it re-checks itself.
STAGE_BATCH_SCRIPT = "task_scripts/run_stage.sbatch"

SUBMISSIONS_BASENAME = "submissions.jsonl"


# -- stage lifecycle --------------------------------------------------------

STATE_PLANNED = "PLANNED"
STATE_READY = "READY"
STATE_BLOCKED = "BLOCKED"
STATE_SUBMITTING = "SUBMITTING"
STATE_QUEUED = "QUEUED"
STATE_RUNNING = "RUNNING"
STATE_VALIDATING = "VALIDATING"
STATE_COMPLETE = "COMPLETE"
STATE_FAILED = "FAILED"
STATE_NOT_APPLICABLE = "NOT_APPLICABLE"

#: States in which a stage already has work in flight or finished.
ACTIVE_STATES = frozenset(
    {STATE_SUBMITTING, STATE_QUEUED, STATE_RUNNING, STATE_VALIDATING}
)

TERMINAL_STATES = frozenset({STATE_COMPLETE, STATE_FAILED})

#: Display order, weakest to strongest. The UI renders exactly these.
ALL_STATES: tuple[str, ...] = (
    STATE_PLANNED,
    STATE_BLOCKED,
    STATE_READY,
    STATE_SUBMITTING,
    STATE_QUEUED,
    STATE_RUNNING,
    STATE_VALIDATING,
    STATE_COMPLETE,
    STATE_FAILED,
    STATE_NOT_APPLICABLE,
)


class SubmissionError(RuntimeError):
    """Raised when a stage must not be submitted."""


class DuplicateSubmission(SubmissionError):
    """Raised when the same run/stage already has work in flight."""


# ---------------------------------------------------------------------------
# The ledger
# ---------------------------------------------------------------------------


@dataclass
class SubmissionLedger:
    """Append-only record of every submission made for one run."""

    run_dir: Path

    @property
    def path(self) -> Path:
        return Path(self.run_dir) / SUBMISSIONS_BASENAME

    def records(self) -> list[dict[str, Any]]:
        if not self.path.is_file():
            return []

        entries: list[dict[str, Any]] = []

        for line in self.path.read_text(encoding="utf-8").splitlines():
            line = line.strip()

            if not line:
                continue

            try:
                loaded = json.loads(line)
            except ValueError:
                continue

            if isinstance(loaded, dict):
                entries.append(loaded)

        return entries

    def for_stage(self, stage_id: str) -> list[dict[str, Any]]:
        return [
            entry
            for entry in self.records()
            if entry.get("stage_id") == stage_id
        ]

    def latest(self, stage_id: str) -> dict[str, Any] | None:
        entries = self.for_stage(stage_id)

        return entries[-1] if entries else None

    def append(self, record: dict[str, Any]) -> dict[str, Any]:
        self.path.parent.mkdir(parents=True, exist_ok=True)

        record.setdefault("recorded_at", utc_now_iso())
        record.setdefault(
            "attempt", len(self.for_stage(record.get("stage_id", ""))) + 1
        )

        with self.path.open("a", encoding="utf-8") as handle:
            handle.write(json.dumps(record, sort_keys=True, default=str) + "\n")

        return record


# ---------------------------------------------------------------------------
# State
# ---------------------------------------------------------------------------


@dataclass
class StageStatus:
    """What one stage of one run is doing right now."""

    stage_id: str
    state: str
    reason: str = ""
    submission: dict[str, Any] | None = None
    slurm: list[dict[str, Any]] = field(default_factory=list)
    validation: str = pipeline_module.PENDING
    execution: dict[str, Any] | None = None
    may_submit: bool = False
    may_retry: bool = False

    def to_dict(self) -> dict[str, Any]:
        submission = dict(self.submission or {})

        return {
            "stage_id": self.stage_id,
            "state": self.state,
            "reason": self.reason,
            "validation": self.validation,
            "may_submit": self.may_submit,
            "may_retry": self.may_retry,
            "slurm_jobs": list(self.slurm),
            "submission": submission,
            "execution": _execution_summary(self.execution),
        }


def _execution_summary(record: dict[str, Any] | None) -> dict[str, Any] | None:
    """The parts of an execution record a UI should show."""

    if not record:
        return None

    return {
        "attempt_id": record.get("attempt_id"),
        "started_at": record.get("started_at"),
        "finished_at": record.get("finished_at"),
        "returncode": record.get("returncode"),
        "preflight_verdict": record.get("preflight_verdict"),
        "preflight_error": record.get("preflight_error"),
        "value_agreement": record.get("value_agreement"),
        "value_mismatches": record.get("value_mismatches", []),
        "executed_values": record.get("executed_values", {}),
        "frozen_values": record.get("frozen_values", {}),
        "resolved_config_sha256": record.get("resolved_config_sha256"),
        "scientific_config_sha256": record.get("scientific_config_sha256"),
        "stage_config_path": record.get("stage_config_path"),
        "pipeline_git_commit": record.get("pipeline_git_commit"),
        "argv": record.get("argv", []),
        "record_path": record.get("record_path"),
        "child_jobs": record.get("child_jobs", []),
        "hostname": record.get("hostname"),
    }


def _slurm_states(
    client: SlurmClient | None,
    job_ids: Sequence[str],
) -> list[dict[str, Any]]:
    if client is None or not client.available or not job_ids:
        return [{"job_id": job_id, "state": "UNKNOWN"} for job_id in job_ids]

    from pdbclean.slurm import states_of

    states = states_of(client, job_ids)

    return [states[job_id].to_dict() for job_id in job_ids]


def _aggregate_slurm(states: Sequence[dict[str, Any]]) -> str:
    """One state for a submission that may have created several jobs."""

    from pdbclean.slurm import (
        FAILURE_STATES,
        QUEUED_STATES,
        RUNNING_STATES,
        SUCCESS_STATES,
    )

    names = {entry.get("state", "UNKNOWN") for entry in states}

    if names & FAILURE_STATES:
        return STATE_FAILED

    if names & RUNNING_STATES:
        return STATE_RUNNING

    if names & QUEUED_STATES:
        return STATE_QUEUED

    if names and names <= SUCCESS_STATES:
        return STATE_VALIDATING

    return STATE_SUBMITTING


def stage_status(
    run: FrozenRun,
    stage_id: str,
    *,
    plan: PipelinePlan,
    client: SlurmClient | None = None,
) -> StageStatus:
    """Combine the plan, the ledger, Slurm and the execution record."""

    from pdbclean.stage_runner import latest_execution_record

    observation = plan.by_id.get(stage_id)
    execution_spec = stage_execution(stage_id)

    ledger = SubmissionLedger(run.run_dir)
    submission = ledger.latest(stage_id)
    execution = latest_execution_record(run.run_dir, stage_id)

    validation = (
        observation.validation if observation else pipeline_module.PENDING
    )

    if execution_spec.mode == EXECUTION_NONE or (
        observation is not None
        and observation.action == pipeline_module.ACTION_NOT_APPLICABLE
    ):
        return StageStatus(
            stage_id=stage_id,
            state=STATE_NOT_APPLICABLE,
            reason=execution_spec.note or "No artefact and no job by design.",
            validation=validation,
        )

    # Already-valid output for this exact configuration is complete, whether
    # this run produced it or an earlier one did.
    if observation is not None and observation.action == (
        pipeline_module.ACTION_REUSE
    ):
        return StageStatus(
            stage_id=stage_id,
            state=STATE_COMPLETE,
            reason=(
                "Validated output already exists for this scientific "
                "configuration."
            ),
            validation=validation,
            submission=submission,
            execution=execution,
        )

    if submission is None:
        blocked = observation is not None and observation.action == (
            pipeline_module.ACTION_BLOCKED
        )

        return StageStatus(
            stage_id=stage_id,
            state=STATE_BLOCKED if blocked else STATE_READY,
            reason=(
                "Upstream validation has not passed."
                if blocked
                else "Eligible; nothing submitted yet."
            ),
            validation=validation,
            may_submit=not blocked,
        )

    job_ids = [
        str(job["job_id"])
        for job in submission.get("jobs", [])
        if job.get("job_id")
    ]

    states = _slurm_states(client, job_ids)
    slurm_state = _aggregate_slurm(states)

    if submission.get("outcome") == "failed":
        return StageStatus(
            stage_id=stage_id,
            state=STATE_FAILED,
            reason=submission.get("error", "Submission failed."),
            submission=submission,
            slurm=states,
            validation=validation,
            execution=execution,
            may_retry=True,
        )

    if slurm_state in {STATE_QUEUED, STATE_RUNNING, STATE_SUBMITTING}:
        return StageStatus(
            stage_id=stage_id,
            state=slurm_state,
            reason=f"Slurm job(s) {', '.join(job_ids)}.",
            submission=submission,
            slurm=states,
            validation=validation,
            execution=execution,
        )

    if slurm_state == STATE_FAILED:
        return StageStatus(
            stage_id=stage_id,
            state=STATE_FAILED,
            reason="Slurm reported a failed job; downstream stages stay blocked.",
            submission=submission,
            slurm=states,
            validation=validation,
            execution=execution,
            may_retry=True,
        )

    # Slurm says the job finished successfully. That is *not* completion.
    problems: list[str] = []

    if execution is None:
        problems.append("no execution provenance was written")
    else:
        if execution.get("preflight_verdict") == "FAIL":
            problems.append(
                "preflight failed: " + str(execution.get("preflight_error"))
            )

        if execution.get("value_agreement") == "FAIL":
            problems.append(
                "the worker loaded values that differ from the frozen "
                "configuration"
            )

        returncode = execution.get("returncode")

        if returncode not in (None, 0):
            problems.append(f"stage exited {returncode}")

    if validation != pipeline_module.VALIDATION_PASS:
        problems.append(
            f"stage validation is {validation!r}, not "
            f"{pipeline_module.VALIDATION_PASS!r}"
        )

    if problems:
        failed = any(
            "preflight failed" in problem
            or "differ from the frozen" in problem
            or problem.startswith("stage exited")
            for problem in problems
        )

        return StageStatus(
            stage_id=stage_id,
            state=STATE_FAILED if failed else STATE_VALIDATING,
            reason="; ".join(problems),
            submission=submission,
            slurm=states,
            validation=validation,
            execution=execution,
            may_retry=failed,
        )

    return StageStatus(
        stage_id=stage_id,
        state=STATE_COMPLETE,
        reason="Slurm succeeded, validation passed, executed values match.",
        submission=submission,
        slurm=states,
        validation=validation,
        execution=execution,
    )


def run_status(
    run: FrozenRun,
    *,
    plan: PipelinePlan,
    client: SlurmClient | None = None,
) -> list[StageStatus]:
    """Status of every stage of one run, in pipeline order."""

    return [
        stage_status(run, observation.stage.stage_id, plan=plan, client=client)
        for observation in plan.observations
    ]


def next_submittable_stage(
    statuses: Sequence[StageStatus],
) -> StageStatus | None:
    """The one stage that may start next, or None.

    Deliberately returns *one*.  "Run all eligible stages" means submit this
    one, wait for it to reach COMPLETE -- Slurm success and validation and
    matching executed values -- and only then ask again.  It never means
    submitting a dependency chain of everything at once, which would bypass
    every validation gate the pipeline has.
    """

    for status in statuses:
        if status.state in ACTIVE_STATES:
            # Something is already in flight; nothing else may start.
            return None

        if status.state == STATE_FAILED:
            # A failure blocks the chain until a human retries it.
            return None

        if status.state == STATE_READY and status.may_submit:
            return status

    return None


# ---------------------------------------------------------------------------
# Submitting
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Submitter:
    """Submits one frozen run's stages."""

    repo_root: Path
    client: SlurmClient
    allow_dirty_worktree: bool = False

    # -- gates ----------------------------------------------------------

    def preflight(self, run: FrozenRun, stage_id: str) -> dict[str, Any]:
        """Everything that must be true before sbatch is called."""

        if stage_id not in STAGES_BY_ID:
            raise SubmissionError(f"Unknown stage: {stage_id!r}")

        execution = stage_execution(stage_id)

        if execution.mode == EXECUTION_NONE:
            raise SubmissionError(
                f"Stage {stage_id!r} is not executed by the orchestrated "
                f"pipeline: {execution.note}"
            )

        if run.is_legacy:
            raise SubmissionError(
                f"Run {run.run_id} was frozen before runs carried an "
                "executable configuration. It is a readable historical record "
                "and cannot be executed; freeze a new run instead."
            )

        try:
            verification = run.verify()
        except FrozenRunError as exc:
            raise SubmissionError(str(exc)) from exc

        # The code that will run must be the code the run was frozen at.
        if self.allow_dirty_worktree:
            from pdbclean.stage_runner import git_state

            state = git_state(self.repo_root)
            head = state["head"]
        else:
            try:
                head = resolve_clean_git_commit(self.repo_root)
            except ProvenanceError as exc:
                raise SubmissionError(
                    f"{exc}. A production stage records a commit as its "
                    "producer, so the worktree must match that commit exactly "
                    "(untracked files count). Commit or stash the changes, or "
                    "submit with dirty-worktree execution explicitly enabled."
                ) from exc

        if run.git_commit and head and head != run.git_commit:
            raise SubmissionError(
                f"Refusing to submit run {run.run_id}: it was frozen at commit "
                f"{run.git_commit} but HEAD is {head}. Submitting now would "
                "execute different code while recording the frozen commit. "
                "Check out the frozen commit, or freeze a new run."
            )

        return {**verification, "head_commit": head}

    def may_submit(
        self,
        run: FrozenRun,
        stage_id: str,
        *,
        plan: PipelinePlan,
        retry: bool = False,
    ) -> StageStatus:
        """Raise unless ``stage_id`` may be submitted right now."""

        status = stage_status(run, stage_id, plan=plan, client=self.client)

        if status.state in ACTIVE_STATES:
            raise DuplicateSubmission(
                f"Stage {stage_id} of run {run.run_id} is already "
                f"{status.state}; refusing to submit a duplicate job. "
                f"({status.reason})"
            )

        if status.state == STATE_COMPLETE:
            raise DuplicateSubmission(
                f"Stage {stage_id} of run {run.run_id} is already COMPLETE. "
                f"({status.reason})"
            )

        if status.state == STATE_BLOCKED:
            raise SubmissionError(
                f"Stage {stage_id} is blocked: {status.reason}"
            )

        if status.state == STATE_NOT_APPLICABLE:
            raise SubmissionError(
                f"Stage {stage_id} is not applicable: {status.reason}"
            )

        if status.state == STATE_FAILED and not retry:
            raise SubmissionError(
                f"Stage {stage_id} of run {run.run_id} previously FAILED "
                f"({status.reason}). Resubmit explicitly as a retry."
            )

        return status

    # -- submission -----------------------------------------------------

    def submit(
        self,
        run: FrozenRun,
        stage_id: str,
        *,
        plan: PipelinePlan,
        retry: bool = False,
    ) -> dict[str, Any]:
        """Verify, submit, and record.  Returns the ledger entry."""

        # The stage id reaches an argument vector, so it is checked against the
        # registry before anything else looks at it.
        if stage_id not in STAGES_BY_ID:
            raise SubmissionError(f"Unknown stage: {stage_id!r}")

        self.may_submit(run, stage_id, plan=plan, retry=retry)

        verification = self.preflight(run, stage_id)

        execution = stage_execution(stage_id)
        ledger = SubmissionLedger(run.run_dir)

        base: dict[str, Any] = {
            "run_id": run.run_id,
            "run_directory": str(run.run_dir),
            "stage_id": stage_id,
            "submitted_at": utc_now_iso(),
            "pipeline_git_commit": run.git_commit,
            "head_commit": verification.get("head_commit"),
            "resolved_config_sha256": run.resolved_config_sha256,
            "scientific_config_sha256": run.scientific_config_sha256,
            "stage_config_path": verification["stage_config_path"],
            "stage_config_sha256": verification["stage_config_sha256"],
            "execution_mode": execution.mode,
            "retry": bool(retry),
            "allow_dirty_worktree": self.allow_dirty_worktree,
        }

        try:
            if execution.mode == EXECUTION_ARRAY_THEN_FINALIZE:
                jobs = self._submit_array_then_finalize(run, stage_id)
            elif execution.mode == EXECUTION_SUBMITTER:
                jobs = self._submit_via_submitter(run, stage_id)
            else:
                jobs = [self._submit_single(run, stage_id)]
        except (SlurmError, SubmissionError, FrozenRunError) as exc:
            record = ledger.append(
                {**base, "outcome": "failed", "error": str(exc), "jobs": []}
            )

            raise SubmissionError(
                f"Submission of stage {stage_id} failed: {exc}"
            ) from exc

        return ledger.append({**base, "outcome": "submitted", "jobs": jobs})

    # -- one job per mode -----------------------------------------------

    def _log_paths(
        self,
        run: FrozenRun,
        stage_id: str,
        role: str,
        *,
        array: bool = False,
    ) -> tuple[Path, Path]:
        directory = run.run_dir / "logs" / stage_id
        pattern = "%A_%a" if array else "%j"

        return (
            directory / f"{role}_{pattern}.out",
            directory / f"{role}_{pattern}.err",
        )

    def _stage_batch_request(
        self,
        run: FrozenRun,
        stage_id: str,
        *,
        role: str,
        dependency: str | None = None,
    ) -> SubmissionRequest:
        from pdbclean.cli import slurm_resources_for

        script = self.repo_root / STAGE_BATCH_SCRIPT
        resources = slurm_resources_for(stage_id, run.resolved)

        stdout, stderr = self._log_paths(run, stage_id, role)

        script_args = [
            str(run.run_dir),
            stage_id,
            str(run.git_commit or ""),
            str(run.resolved_config_sha256 or ""),
            str(run.scientific_config_sha256 or ""),
        ]

        if self.allow_dirty_worktree:
            script_args.append("--allow-dirty-worktree")

        return SubmissionRequest(
            script=script,
            script_args=tuple(script_args),
            job_name=f"pdbclean-{stage_id}"[:64].replace("_", "-"),
            partition=resources.partition,
            time_limit=resources.time_limit,
            memory=resources.memory,
            cpus=resources.cpus,
            dependency_afterok=dependency,
            stdout_path=stdout,
            stderr_path=stderr,
        )

    def _submit_single(self, run: FrozenRun, stage_id: str) -> dict[str, Any]:
        request = self._stage_batch_request(run, stage_id, role="stage")

        result = self.client.submit(request, cwd=self.repo_root)

        stdout, stderr = self._log_paths(run, stage_id, "stage")

        return {
            "role": "stage",
            "job_id": result["job_id"],
            "argv": result["argv"],
            "stdout_path": str(stdout),
            "stderr_path": str(stderr),
        }

    def _submit_array_then_finalize(
        self,
        run: FrozenRun,
        stage_id: str,
    ) -> list[dict[str, Any]]:
        """Submit the per-task array, then its finalizer with ``afterok``.

        ``afterok`` here is not a way around the validation gate: it chains the
        two halves of *one* stage.  The stage is still only COMPLETE when its
        finalizer's own validation passes, and no downstream stage is submitted
        by this dependency.
        """

        execution = stage_execution(stage_id)
        shape = array_shape(run, self.repo_root)

        from pdbclean.cli import slurm_resources_for

        array_resources = slurm_resources_for(stage_id, run.resolved, array=True)

        array_stdout, array_stderr = self._log_paths(
            run, stage_id, "array", array=True
        )

        array_request = SubmissionRequest(
            script=self.repo_root / str(execution.array_script),
            script_args=(
                str(run.stage_config_path),
                shape["manifest_path"],
                str(shape["task_count"]),
                str(shape["worker_count"]),
                str(self.repo_root),
            ),
            job_name=f"pdbclean-{stage_id}-array"[:64].replace("_", "-"),
            partition=array_resources.partition,
            time_limit=array_resources.time_limit,
            memory=array_resources.memory,
            cpus=array_resources.cpus,
            array=f"0-{shape['worker_count'] - 1}%{shape['concurrency']}",
            stdout_path=array_stdout,
            stderr_path=array_stderr,
        )

        array_result = self.client.submit(array_request, cwd=self.repo_root)

        finalize_request = self._stage_batch_request(
            run,
            stage_id,
            role="finalize",
            dependency=array_result["job_id"],
        )

        finalize_result = self.client.submit(
            finalize_request, cwd=self.repo_root
        )

        finalize_stdout, finalize_stderr = self._log_paths(
            run, stage_id, "finalize"
        )

        return [
            {
                "role": "array",
                "job_id": array_result["job_id"],
                "argv": array_result["argv"],
                "array_specification": array_request.array,
                "array_task_count": shape["task_count"],
                "array_worker_count": shape["worker_count"],
                "stdout_path": str(array_stdout),
                "stderr_path": str(array_stderr),
            },
            {
                "role": "finalize",
                "job_id": finalize_result["job_id"],
                "argv": finalize_result["argv"],
                "depends_on": array_result["job_id"],
                "stdout_path": str(finalize_stdout),
                "stderr_path": str(finalize_stderr),
            },
        ]

    def _submit_via_submitter(
        self,
        run: FrozenRun,
        stage_id: str,
    ) -> list[dict[str, Any]]:
        """Run a stage whose entry point is itself a submitter.

        ``submit_quality_pipeline.sh`` reads the manifest, works out the array
        shape and submits the array plus its merge.  It computes no science, so
        it runs here rather than occupying a compute node for the seconds it
        takes -- exactly what the documented manual workflow does today.
        """

        from pdbclean.stage_runner import StageRunnerError, run_stage

        try:
            returncode, record = run_stage(
                run,
                stage_id,
                repo_root=self.repo_root,
                allow_dirty=self.allow_dirty_worktree,
            )
        except (StageRunnerError, FrozenRunError) as exc:
            raise SubmissionError(str(exc)) from exc

        if returncode != 0:
            raise SubmissionError(
                f"{stage_id} submitter exited {returncode}: "
                f"{record.get('submitter_stderr_tail', '')[-500:]}"
            )

        children = record.get("child_jobs") or []

        if not children:
            raise SubmissionError(
                f"{stage_id} submitter reported no child jobs. Nothing is "
                "running; refusing to record a submission that did not happen."
            )

        return [
            {
                "role": child["role"],
                "job_id": child["job_id"],
                "argv": record.get("argv", []),
                "submitted_by": "submit_quality_pipeline.sh",
                "stdout_path": None,
                "stderr_path": None,
            }
            for child in children
        ]


# ---------------------------------------------------------------------------
# Array shape
# ---------------------------------------------------------------------------


def array_shape(run: FrozenRun, repo_root: Path) -> dict[str, Any]:
    """Derive the array shape exactly as submit_quality_pipeline.sh does.

    Logical tasks come from the manifest row count and the configured batch
    size; physical workers stride over logical tasks because Barkla limits how
    many array elements one user may have queued.
    """

    import pyarrow.parquet as pq

    from pdbclean.manifest import manifest_partition_count

    resolved = run.resolved
    paths = PipelinePaths.from_config(resolved, repo_root=repo_root)

    manifest = (
        paths.output_root / paths.snapshot / "bronze/source_manifest.parquet"
    )

    if not manifest.is_file():
        raise SubmissionError(
            f"Cannot size the array: no source manifest at {manifest}. "
            "The Bronze manifest stage must complete first."
        )

    row_count = pq.read_metadata(str(manifest)).num_rows

    batch_size = int(resolved.get("execution.batch_size") or 500)
    task_count = manifest_partition_count(row_count, batch_size)

    worker_limit = int(resolved.get("execution.quality_array_worker_count") or 64)
    concurrency_limit = int(
        resolved.get("execution.quality_array_concurrency") or 4
    )

    worker_count = max(1, min(worker_limit, task_count))
    concurrency = max(1, min(concurrency_limit, worker_count))

    return {
        "manifest_path": str(manifest),
        "manifest_rows": row_count,
        "batch_size": batch_size,
        "task_count": task_count,
        "worker_count": worker_count,
        "concurrency": concurrency,
    }


def submittable_stage_ids() -> list[str]:
    """Stage ids the orchestrator will submit. A strict allowlist."""

    return [
        stage.stage_id
        for stage in sorted(STAGES_BY_ID.values(), key=lambda s: s.ordinal)
        if stage_execution(stage.stage_id).mode != EXECUTION_NONE
    ]


def cancel_stage(
    run: FrozenRun,
    stage_id: str,
    *,
    client: SlurmClient,
) -> list[str]:
    """Cancel every job the latest submission for this stage created."""

    submission = SubmissionLedger(run.run_dir).latest(stage_id)

    if submission is None:
        return []

    cancelled: list[str] = []

    for job in submission.get("jobs", []):
        job_id = job.get("job_id")

        if not job_id:
            continue

        client.cancel(str(job_id))
        cancelled.append(str(job_id))

    return cancelled
