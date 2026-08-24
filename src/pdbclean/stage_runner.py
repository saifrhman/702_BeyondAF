"""The worker-side gate: prove the frozen configuration, then execute it.

Every orchestrated stage of every run goes through this module, and nothing
scientific happens until it has established, on the executing host:

1. the run's canonical configuration still hashes to what the run recorded;
2. the executable projection on disk is still the projection of that canonical
   configuration, byte for byte;
3. the checkout is the frozen commit, and the worktree is clean;
4. the values the *production loader* extracts from that projection are the
   values the run was frozen with -- Brain threshold, complete-BRI threshold,
   representation precision, snapshot, model scope, Q005 distance, minimum
   N-CA-C angle.

Only then is the stage's own argv built and executed.  Any failure aborts
before the stage starts, and is recorded.

Point (4) is the one that matters most.  It is deliberately performed by
loading ``stage_config.yaml`` through :func:`pdbclean.config.load_config` --
the same call the stage itself makes -- and reading the thresholds back through
the same helpers and the same fallbacks the stage uses.  So this is not a
restatement of what the browser displayed; it is a measurement of what the
worker will compute with.  If the projection ever stopped carrying
``brain_filter``, this check would observe the module constant 0.010 A and fail
the stage instead of publishing science under the wrong threshold.

Usage::

    python -m pdbclean.stage_runner \\
        --run-dir outputs/runs/run-... \\
        --stage candidate_filtering \\
        --expect-commit <sha1> \\
        --expect-resolved-config-sha256 <sha256> \\
        --expect-scientific-config-sha256 <sha256>
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import re
import socket
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from pdbclean.config import ConfigError, load_config
from pdbclean.run_provenance import FrozenRun, FrozenRunError, utc_now_iso
from pdbclean.stage_config import (
    compare_values,
    executed_values,
    frozen_values,
)
from pdbclean.stage_registry import (
    EXECUTION_NONE,
    EXECUTION_SUBMITTER,
    STAGES_BY_ID,
    stage_execution,
)


EXECUTION_SCHEMA_NAME = "pdbclean_stage_execution"
EXECUTION_SCHEMA_VERSION = "1.0"

#: Exit code used when the frozen configuration cannot be trusted.  Distinct
#: from any exit code a scientific stage produces, so a preflight abort is
#: never mistaken for a scientific failure.
EXIT_PREFLIGHT_FAILED = 78

#: Machine-readable line a submitter stage prints for each job it creates.
CHILD_JOB_MARKER = re.compile(
    r"^PDBCLEAN_CHILD_JOB\s+(?P<role>[a-z_]+)\s+(?P<job_id>\d+)\s*$"
)


class StageRunnerError(RuntimeError):
    """Raised when a stage must not execute."""


def _git(repo_root: Path, *args: str) -> str | None:
    try:
        completed = subprocess.run(
            ["git", "-C", str(repo_root), *args],
            capture_output=True,
            text=True,
            check=False,
            timeout=60,
        )
    except (OSError, subprocess.SubprocessError):
        return None

    return completed.stdout if completed.returncode == 0 else None


def git_state(repo_root: Path) -> dict[str, Any]:
    """HEAD, dirtiness, and enough detail to reproduce a dirty tree."""

    head = (_git(repo_root, "rev-parse", "HEAD") or "").strip() or None
    status = _git(repo_root, "status", "--porcelain=v1", "--untracked-files=all")

    dirty_paths = (
        [line for line in status.splitlines() if line.strip()]
        if status is not None
        else []
    )

    payload: dict[str, Any] = {
        "head": head,
        "working_tree_dirty": None if status is None else bool(dirty_paths),
        "dirty_paths": dirty_paths,
    }

    if dirty_paths:
        # A dirty tree means the commit alone does not describe the code that
        # ran. Record a digest of the actual difference so the executed source
        # is at least identifiable.
        diff = _git(repo_root, "diff", "HEAD") or ""

        payload["uncommitted_diff_sha256"] = hashlib.sha256(
            diff.encode("utf-8")
        ).hexdigest()
        payload["uncommitted_diff_bytes"] = len(diff.encode("utf-8"))

    return payload


def verify_git(
    repo_root: Path,
    *,
    frozen_commit: str | None,
    allow_dirty: bool,
) -> dict[str, Any]:
    """Refuse to run different code while recording the frozen commit."""

    state = git_state(repo_root)

    if frozen_commit and state["head"] and state["head"] != frozen_commit:
        raise StageRunnerError(
            "Checkout does not match the run's frozen commit.\n"
            f"  frozen: {frozen_commit}\n"
            f"  HEAD:   {state['head']}\n"
            "Executing different code while recording the frozen commit would "
            "make the provenance false. Check out the frozen commit, or freeze "
            "a new run at the current commit."
        )

    if state["working_tree_dirty"] and not allow_dirty:
        listing = "\n".join(f"    {line}" for line in state["dirty_paths"][:20])

        raise StageRunnerError(
            "Refusing a production run from a dirty worktree: a commit alone "
            "would not describe the code that ran.\n"
            f"{listing}\n"
            "Commit or stash the changes (untracked files count), or pass "
            "--allow-dirty-worktree to record the difference explicitly."
        )

    state["dirty_execution_allowed"] = bool(
        allow_dirty and state["working_tree_dirty"]
    )

    return state


def build_execution_record(
    run: FrozenRun,
    stage_id: str,
    *,
    repo_root: Path,
    allow_dirty: bool,
) -> dict[str, Any]:
    """Perform every preflight check and return the execution record.

    Raises :class:`StageRunnerError` when the stage must not run.  The record
    is returned either way through the exception path's caller, so a refusal is
    itself recorded rather than being a silent non-event.
    """

    record: dict[str, Any] = {
        "schema_name": EXECUTION_SCHEMA_NAME,
        "schema_version": EXECUTION_SCHEMA_VERSION,
        "run_id": run.run_id,
        "run_directory": str(run.run_dir),
        "stage_id": stage_id,
        "started_at": utc_now_iso(),
        "hostname": socket.gethostname(),
        "python_executable": sys.executable,
        "python_version": platform.python_version(),
        "slurm": {
            key: os.environ[key]
            for key in (
                "SLURM_JOB_ID",
                "SLURM_ARRAY_JOB_ID",
                "SLURM_ARRAY_TASK_ID",
                "SLURM_JOB_PARTITION",
                "SLURM_JOB_NODELIST",
            )
            if key in os.environ
        },
        "preflight": {},
    }

    # (1) + (2): the frozen configuration is intact and the projection is its
    # projection. FrozenRun.verify() recomputes both canonical hashes and
    # re-projects the executable document.
    verification = run.verify()

    record["preflight"]["frozen_configuration"] = verification
    record["resolved_config_sha256"] = verification["resolved_config_sha256"]
    record["scientific_config_sha256"] = verification[
        "scientific_config_sha256"
    ]
    record["stage_config_path"] = verification["stage_config_path"]
    record["stage_config_sha256"] = verification["stage_config_sha256"]

    # (3): the code.
    record["git"] = verify_git(
        repo_root,
        frozen_commit=run.git_commit,
        allow_dirty=allow_dirty,
    )
    record["git"]["frozen_commit"] = run.git_commit

    # (4): the values the production loader gives the worker.
    try:
        loaded = load_config(verification["stage_config_path"])
    except ConfigError as exc:
        raise StageRunnerError(
            "The frozen stage configuration is not a valid Protocol 3.2 "
            f"configuration: {exc}"
        ) from exc

    resolved = run.resolved

    executed = executed_values(loaded.data)
    frozen = frozen_values(resolved)
    mismatches = compare_values(executed, frozen)

    record["executed_values"] = executed
    record["frozen_values"] = frozen
    record["value_mismatches"] = mismatches
    record["value_agreement"] = "FAIL" if mismatches else "PASS"
    record["stage_config_file_sha256"] = hashlib.sha256(
        Path(verification["stage_config_path"]).read_bytes()
    ).hexdigest()

    if mismatches:
        detail = "\n".join(
            f"    {item['key']}: worker loaded {item['executed']!r}, "
            f"run was frozen with {item['frozen']!r}"
            for item in mismatches
        )

        raise StageRunnerError(
            "The values this worker loaded are not the values the run was "
            f"frozen with:\n{detail}\n"
            "Refusing to compute science under a configuration the operator "
            "did not review."
        )

    return record


def write_execution_record(
    run: FrozenRun,
    stage_id: str,
    record: dict[str, Any],
) -> Path:
    """Append-only execution provenance, one file per attempt."""

    directory = run.run_dir / "execution" / stage_id
    directory.mkdir(parents=True, exist_ok=True)

    attempt = record.setdefault(
        "attempt_id",
        "{stamp}-{job}".format(
            stamp=datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ"),
            job=record.get("slurm", {}).get("SLURM_JOB_ID", "local"),
        ),
    )

    path = directory / f"{attempt}.json"

    temporary = path.with_name(f".{path.name}.tmp")
    temporary.write_text(
        json.dumps(record, indent=2, sort_keys=True, default=str) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)

    return path


def latest_execution_record(
    run_dir: Path,
    stage_id: str,
) -> dict[str, Any] | None:
    """The most recent execution record for one stage of one run."""

    directory = Path(run_dir) / "execution" / stage_id

    if not directory.is_dir():
        return None

    candidates = sorted(
        path for path in directory.glob("*.json") if path.is_file()
    )

    if not candidates:
        return None

    try:
        loaded = json.loads(candidates[-1].read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None

    return loaded if isinstance(loaded, dict) else None


def run_stage(
    run: FrozenRun,
    stage_id: str,
    *,
    repo_root: Path,
    allow_dirty: bool = False,
    python: str | None = None,
    dry_run: bool = False,
) -> tuple[int, dict[str, Any]]:
    """Verify, record, and execute one stage.  Returns (exit code, record)."""

    from pdbclean.cli import stage_invocation

    if stage_id not in STAGES_BY_ID:
        raise StageRunnerError(f"Unknown stage: {stage_id!r}")

    execution = stage_execution(stage_id)

    if execution.mode == EXECUTION_NONE:
        raise StageRunnerError(
            f"Stage {stage_id!r} is not executed by the orchestrated "
            f"pipeline: {execution.note}"
        )

    record: dict[str, Any] = {}

    try:
        record = build_execution_record(
            run,
            stage_id,
            repo_root=repo_root,
            allow_dirty=allow_dirty,
        )
    except (StageRunnerError, FrozenRunError) as exc:
        record = {
            "schema_name": EXECUTION_SCHEMA_NAME,
            "schema_version": EXECUTION_SCHEMA_VERSION,
            "run_id": run.run_id,
            "stage_id": stage_id,
            "started_at": utc_now_iso(),
            "hostname": socket.gethostname(),
            "preflight_verdict": "FAIL",
            "preflight_error": str(exc),
        }

        write_execution_record(run, stage_id, record)

        raise

    invocation = stage_invocation(
        stage_id,
        run,
        repo_root=repo_root,
        python=python,
        # Already verified above; verifying twice would only re-read files.
        verify=False,
    )

    if invocation is None:
        raise StageRunnerError(
            f"Stage {stage_id!r} has no executable command registered."
        )

    record["argv"] = list(invocation.argv)
    record["command_text"] = invocation.command_text
    record["pipeline_git_commit"] = invocation.pipeline_git_commit
    record["execution_mode"] = execution.mode
    record["preflight_verdict"] = "PASS"

    if dry_run:
        record["dry_run"] = True
        record["finished_at"] = utc_now_iso()
        record["record_path"] = str(write_execution_record(run, stage_id, record))

        return 0, record

    # Written before the work starts: a job that dies without warning still
    # leaves behind exactly what it was about to compute with.
    record_path = write_execution_record(run, stage_id, record)
    record["record_path"] = str(record_path)

    print("=" * 70, flush=True)
    print(f"PDBClean stage        {stage_id}", flush=True)
    print(f"Run                   {run.run_id}", flush=True)
    print(f"Frozen commit         {invocation.pipeline_git_commit}", flush=True)
    print(f"Resolved config sha   {invocation.resolved_config_sha256}", flush=True)
    print(f"Scientific sha        {invocation.scientific_config_sha256}", flush=True)
    print(f"Stage configuration   {invocation.config_path}", flush=True)
    print("Executed values:", flush=True)

    for key in sorted(record["executed_values"]):
        print(f"  {key:<52} {record['executed_values'][key]}", flush=True)

    print(f"Execution provenance  {record_path}", flush=True)
    print("=" * 70, flush=True)
    print(invocation.command_text, flush=True)
    print("", flush=True)

    capture = execution.mode == EXECUTION_SUBMITTER

    completed = subprocess.run(
        list(invocation.argv),
        cwd=str(repo_root),
        check=False,
        capture_output=capture,
        text=True,
    )

    record["returncode"] = completed.returncode
    record["finished_at"] = utc_now_iso()

    if capture:
        # A submitter creates the jobs that do the real work; the run has to
        # know their identifiers to report on them.
        print(completed.stdout or "", flush=True)
        print(completed.stderr or "", file=sys.stderr, flush=True)

        record["child_jobs"] = _child_jobs(completed.stdout or "")
        record["submitter_stdout_tail"] = (completed.stdout or "")[-8000:]
        record["submitter_stderr_tail"] = (completed.stderr or "")[-8000:]

    write_execution_record(run, stage_id, record)

    return completed.returncode, record


def _child_jobs(stdout: str) -> list[dict[str, str]]:
    jobs: list[dict[str, str]] = []

    for line in stdout.splitlines():
        match = CHILD_JOB_MARKER.match(line.strip())

        if match:
            jobs.append(
                {"role": match.group("role"), "job_id": match.group("job_id")}
            )

    return jobs


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m pdbclean.stage_runner",
        description=(
            "Verify a frozen run's configuration and execute one of its "
            "stages. Refuses to run if anything about the frozen "
            "configuration, the commit or the loaded values has changed."
        ),
    )

    parser.add_argument("--run-dir", required=True, type=Path)
    parser.add_argument("--stage", required=True)
    parser.add_argument("--repo-root", default=None, type=Path)
    parser.add_argument("--python", default=None)

    parser.add_argument(
        "--expect-commit",
        default=None,
        help="Abort unless the run's frozen commit is exactly this.",
    )
    parser.add_argument(
        "--expect-resolved-config-sha256",
        default=None,
        help="Abort unless the run's resolved-config SHA256 is exactly this.",
    )
    parser.add_argument(
        "--expect-scientific-config-sha256",
        default=None,
        help="Abort unless the run's scientific SHA256 is exactly this.",
    )

    parser.add_argument(
        "--allow-dirty-worktree",
        action="store_true",
        help=(
            "Permit execution from a dirty worktree. Provenance then records "
            "the changed paths and a SHA256 of the uncommitted diff, because "
            "the commit alone no longer describes the code that ran."
        ),
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Verify and record, but do not execute the stage.",
    )

    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    repo_root = (
        args.repo_root
        if args.repo_root is not None
        else Path(os.environ.get("PDBCLEAN_REPO_ROOT") or Path.cwd())
    ).resolve()

    try:
        run = FrozenRun.load(args.run_dir)
    except FrozenRunError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return EXIT_PREFLIGHT_FAILED

    # The submitted job itself pins what it expects. A mismatch here means the
    # run directory is not the run that was submitted.
    for expected, observed, what in (
        (args.expect_commit, run.git_commit, "frozen Git commit"),
        (
            args.expect_resolved_config_sha256,
            run.resolved_config_sha256,
            "resolved-config SHA256",
        ),
        (
            args.expect_scientific_config_sha256,
            run.scientific_config_sha256,
            "scientific SHA256",
        ),
    ):
        if expected is not None and expected != observed:
            print(
                f"error: submitted job expected {what} {expected}, but run "
                f"{run.run_id} records {observed}. Refusing to execute.",
                file=sys.stderr,
            )
            return EXIT_PREFLIGHT_FAILED

    try:
        returncode, _ = run_stage(
            run,
            args.stage,
            repo_root=repo_root,
            allow_dirty=args.allow_dirty_worktree,
            python=args.python,
            dry_run=args.dry_run,
        )
    except (StageRunnerError, FrozenRunError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return EXIT_PREFLIGHT_FAILED

    return returncode


if __name__ == "__main__":
    raise SystemExit(main())
