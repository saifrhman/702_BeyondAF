"""A small, safe Slurm client for the PDBClean orchestrator.

Design rules, in order of importance:

1. **No shell.**  Every command is an argument vector passed to
   ``subprocess.run`` with the binary resolved to an absolute path.  There is
   no string interpolation anywhere in this module, so there is nothing for a
   quoting mistake or a hostile value to escape from.

2. **Nothing user-supplied reaches an argument unvalidated.**  The only
   variable parts of a submission are a run directory, a stage id, a commit
   and two SHA256 digests.  Each is checked against a strict pattern here,
   *again*, even though the callers already validated them.

3. **The web server submits; it never computes.**  ``sbatch`` returns in
   milliseconds.  The scientific work happens on a compute node.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


class SlurmError(RuntimeError):
    """Raised when a Slurm operation fails."""


class SlurmUnavailable(SlurmError):
    """Raised when this host has no Slurm client binaries."""


#: A Slurm job id, or one array element of one.
JOB_ID_PATTERN = re.compile(r"\d+(?:_\d+)?")

#: Accepted forms for the values this module puts on an sbatch command line.
STAGE_ID_PATTERN = re.compile(r"[a-z][a-z0-9_]{0,63}")
COMMIT_PATTERN = re.compile(r"[0-9a-f]{40}")
SHA256_PATTERN = re.compile(r"[0-9a-f]{64}")
ARRAY_SPEC_PATTERN = re.compile(r"\d+(?:-\d+)?(?:%\d+)?")
PARTITION_PATTERN = re.compile(r"[A-Za-z0-9][A-Za-z0-9._-]{0,63}")
TIME_LIMIT_PATTERN = re.compile(r"(?:\d+-)?\d{1,2}:\d{2}:\d{2}")
MEMORY_PATTERN = re.compile(r"\d+[KMGT]?")
JOB_NAME_PATTERN = re.compile(r"[A-Za-z0-9][A-Za-z0-9._-]{0,63}")


#: Slurm states that mean the job has not started yet.
QUEUED_STATES = frozenset(
    {"PENDING", "CONFIGURING", "REQUEUED", "RESIZING", "SUSPENDED"}
)

#: Slurm states that mean the job is executing.
RUNNING_STATES = frozenset({"RUNNING", "COMPLETING", "STAGE_OUT"})

#: The only Slurm state that permits a stage to be *considered* for
#: validation.  Note "considered": exit code zero is not scientific validation.
SUCCESS_STATES = frozenset({"COMPLETED"})

#: Terminal failure.
FAILURE_STATES = frozenset(
    {
        "FAILED",
        "CANCELLED",
        "TIMEOUT",
        "NODE_FAIL",
        "OUT_OF_MEMORY",
        "PREEMPTED",
        "BOOT_FAIL",
        "DEADLINE",
        "REVOKED",
        "SPECIAL_EXIT",
    }
)


def _require(pattern: re.Pattern[str], value: str, what: str) -> str:
    text = str(value)

    if not pattern.fullmatch(text):
        raise SlurmError(f"Refusing to submit: invalid {what}: {text!r}")

    return text


@dataclass(frozen=True)
class SlurmJobState:
    """What Slurm currently says about one job."""

    job_id: str
    state: str = "UNKNOWN"
    exit_code: str | None = None
    reason: str | None = None
    submitted_at: str | None = None
    started_at: str | None = None
    finished_at: str | None = None
    elapsed: str | None = None
    node_list: str | None = None
    array_tasks: list[dict[str, Any]] = field(default_factory=list)

    @property
    def queued(self) -> bool:
        return self.state in QUEUED_STATES

    @property
    def running(self) -> bool:
        return self.state in RUNNING_STATES

    @property
    def succeeded(self) -> bool:
        return self.state in SUCCESS_STATES

    @property
    def failed(self) -> bool:
        return self.state in FAILURE_STATES

    @property
    def terminal(self) -> bool:
        return self.succeeded or self.failed

    def to_dict(self) -> dict[str, Any]:
        return {
            "job_id": self.job_id,
            "state": self.state,
            "exit_code": self.exit_code,
            "reason": self.reason,
            "submitted_at": self.submitted_at,
            "started_at": self.started_at,
            "finished_at": self.finished_at,
            "elapsed": self.elapsed,
            "node_list": self.node_list,
            "array_tasks": list(self.array_tasks),
        }


@dataclass(frozen=True)
class SubmissionRequest:
    """One validated sbatch submission."""

    script: Path
    script_args: tuple[str, ...] = ()
    job_name: str = "pdbclean"
    partition: str = "nodes"
    time_limit: str = "02:00:00"
    memory: str = "16G"
    cpus: int = 1
    array: str | None = None
    dependency_afterok: str | None = None
    stdout_path: Path | None = None
    stderr_path: Path | None = None

    def sbatch_argv(self, sbatch: str) -> list[str]:
        """Return the exact argv, with every value re-validated here."""

        argv: list[str] = [
            sbatch,
            "--parsable",
            "--job-name",
            _require(JOB_NAME_PATTERN, self.job_name, "job name"),
            "--partition",
            _require(PARTITION_PATTERN, self.partition, "partition"),
            "--time",
            _require(TIME_LIMIT_PATTERN, self.time_limit, "time limit"),
            "--mem",
            _require(MEMORY_PATTERN, self.memory, "memory request"),
            "--cpus-per-task",
            str(int(self.cpus)),
        ]

        if self.array is not None:
            argv += [
                "--array",
                _require(ARRAY_SPEC_PATTERN, self.array, "array specification"),
            ]

        if self.dependency_afterok is not None:
            argv += [
                "--dependency",
                "afterok:"
                + _require(
                    JOB_ID_PATTERN, self.dependency_afterok, "dependency job id"
                ),
            ]

        if self.stdout_path is not None:
            argv += ["--output", str(self.stdout_path)]

        if self.stderr_path is not None:
            argv += ["--error", str(self.stderr_path)]

        script = Path(self.script)

        if not script.is_file():
            raise SlurmError(f"Batch script does not exist: {script}")

        argv.append(str(script))
        argv.extend(str(part) for part in self.script_args)

        return argv


class SlurmClient:
    """Thin wrapper over ``sbatch`` / ``squeue`` / ``sacct`` / ``scancel``."""

    def __init__(
        self,
        *,
        binaries: Mapping[str, str] | None = None,
        timeout_seconds: int = 60,
    ) -> None:
        self.timeout_seconds = timeout_seconds

        supplied = dict(binaries or {})

        # An explicitly supplied key wins even when its value is None, so a
        # caller can model "this host has no sbatch" exactly.
        self._binaries = {
            name: (
                supplied[name] if name in supplied else shutil.which(name)
            )
            for name in ("sbatch", "squeue", "sacct", "scancel")
        }

    # -- availability ---------------------------------------------------

    @property
    def available(self) -> bool:
        return self._binaries.get("sbatch") is not None

    def binary(self, name: str) -> str:
        path = self._binaries.get(name)

        if not path:
            raise SlurmUnavailable(
                f"This host has no {name}; Slurm submission is not available "
                "here. Run the UI on a Barkla login node, or use the printed "
                "command from a host that has one."
            )

        return path

    def describe(self) -> dict[str, Any]:
        return {
            "available": self.available,
            "binaries": dict(self._binaries),
        }

    # -- running commands -----------------------------------------------

    def _run(self, argv: Sequence[str], *, cwd: Path | None = None):
        try:
            return subprocess.run(
                list(argv),
                cwd=None if cwd is None else str(cwd),
                capture_output=True,
                text=True,
                check=False,
                timeout=self.timeout_seconds,
                # A submitted job must not inherit an array task identity from
                # whatever happens to be in this process' environment.
                env={
                    key: value
                    for key, value in os.environ.items()
                    if not key.startswith("SLURM_")
                },
            )
        except (OSError, subprocess.SubprocessError) as exc:
            raise SlurmError(f"{argv[0]} could not be run: {exc}") from exc

    # -- operations -----------------------------------------------------

    def submit(
        self,
        request: SubmissionRequest,
        *,
        cwd: Path,
    ) -> dict[str, Any]:
        """Submit one batch job and return its identifier and argv."""

        argv = request.sbatch_argv(self.binary("sbatch"))

        for path in (request.stdout_path, request.stderr_path):
            if path is not None:
                Path(path).parent.mkdir(parents=True, exist_ok=True)

        completed = self._run(argv, cwd=cwd)

        if completed.returncode != 0:
            raise SlurmError(
                "sbatch failed with exit code "
                f"{completed.returncode}: {completed.stderr.strip() or completed.stdout.strip()}"
            )

        # `--parsable` prints "<jobid>" or "<jobid>;<cluster>".
        first = completed.stdout.strip().splitlines()
        job_id = first[0].split(";", 1)[0].strip() if first else ""

        if not JOB_ID_PATTERN.fullmatch(job_id):
            raise SlurmError(
                f"sbatch did not return a job id: {completed.stdout!r}"
            )

        return {
            "job_id": job_id,
            "argv": list(argv),
            "stdout": completed.stdout.strip(),
        }

    def state(self, job_id: str) -> SlurmJobState:
        """Return the job's current state.

        ``squeue`` is asked first because a job that is still queued or running
        may not be in the accounting database yet; ``sacct`` answers for jobs
        that have already finished.
        """

        identifier = _require(JOB_ID_PATTERN, job_id, "job id")

        live = self._squeue(identifier)

        if live is not None:
            return live

        return self._sacct(identifier)

    def _squeue(self, job_id: str) -> SlurmJobState | None:
        completed = self._run(
            [
                self.binary("squeue"),
                "--job",
                job_id,
                "--noheader",
                "--array",
                "--Format=JobID:|,State:|,Reason:|,StartTime:|,NodeList:|,TimeUsed:|",
            ]
        )

        if completed.returncode != 0 or not completed.stdout.strip():
            return None

        tasks: list[dict[str, Any]] = []
        parent: SlurmJobState | None = None

        for line in completed.stdout.strip().splitlines():
            fields = [part.strip() for part in line.split("|")]

            if len(fields) < 6:
                continue

            identifier, state, reason, start, nodes, used = fields[:6]

            tasks.append({"job_id": identifier, "state": state})

            if parent is None:
                parent = SlurmJobState(
                    job_id=job_id,
                    state=state,
                    reason=reason or None,
                    started_at=start or None,
                    node_list=nodes or None,
                    elapsed=used or None,
                )

        if parent is None:
            return None

        # An array is only "running" while any element still is.
        states = {task["state"] for task in tasks}

        aggregate = parent.state

        if states & RUNNING_STATES:
            aggregate = "RUNNING"
        elif states & QUEUED_STATES:
            aggregate = "PENDING"

        return SlurmJobState(
            job_id=parent.job_id,
            state=aggregate,
            reason=parent.reason,
            started_at=parent.started_at,
            node_list=parent.node_list,
            elapsed=parent.elapsed,
            array_tasks=tasks if len(tasks) > 1 else [],
        )

    def _sacct(self, job_id: str) -> SlurmJobState:
        completed = self._run(
            [
                self.binary("sacct"),
                "--jobs",
                job_id,
                "--noheader",
                "--parsable2",
                "--format=JobID,State,ExitCode,Submit,Start,End,Elapsed,NodeList",
            ]
        )

        if completed.returncode != 0:
            raise SlurmError(
                f"sacct failed for job {job_id}: {completed.stderr.strip()}"
            )

        rows = [
            [part.strip() for part in line.split("|")]
            for line in completed.stdout.strip().splitlines()
            if line.strip()
        ]

        if not rows:
            return SlurmJobState(job_id=job_id, state="UNKNOWN")

        # Rows for ".batch" / ".extern" steps describe the step, not the job.
        job_rows = [row for row in rows if "." not in row[0]]

        if not job_rows:
            job_rows = rows

        tasks = [
            {
                "job_id": row[0],
                "state": row[1].split()[0] if row[1] else "UNKNOWN",
                "exit_code": row[2] or None,
            }
            for row in job_rows
        ]

        head = job_rows[0]
        states = {task["state"] for task in tasks}

        # An array succeeds only when every element does.
        if states & FAILURE_STATES:
            aggregate = sorted(states & FAILURE_STATES)[0]
        elif states & RUNNING_STATES:
            aggregate = "RUNNING"
        elif states & QUEUED_STATES:
            aggregate = "PENDING"
        elif states <= SUCCESS_STATES:
            aggregate = "COMPLETED"
        else:
            aggregate = sorted(states)[0]

        return SlurmJobState(
            job_id=job_id,
            state=aggregate,
            exit_code=head[2] or None,
            submitted_at=head[3] or None,
            started_at=head[4] or None,
            finished_at=head[5] or None,
            elapsed=head[6] or None,
            node_list=head[7] if len(head) > 7 else None,
            array_tasks=tasks if len(tasks) > 1 else [],
        )

    def cancel(self, job_id: str) -> None:
        identifier = _require(JOB_ID_PATTERN, job_id, "job id")

        completed = self._run([self.binary("scancel"), identifier])

        if completed.returncode != 0:
            raise SlurmError(
                f"scancel failed for job {identifier}: "
                f"{completed.stderr.strip()}"
            )


def states_of(
    client: SlurmClient,
    job_ids: Iterable[str],
) -> dict[str, SlurmJobState]:
    """Return the state of several jobs, tolerating individual failures."""

    result: dict[str, SlurmJobState] = {}

    for job_id in job_ids:
        try:
            result[job_id] = client.state(job_id)
        except SlurmError as exc:
            result[job_id] = SlurmJobState(
                job_id=job_id, state="UNKNOWN", reason=str(exc)
            )

    return result
