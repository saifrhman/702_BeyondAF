"""First-class run provenance for the COMP702 PDBClean pipeline.

Every run gets its own directory under ``storage.run_root`` containing:

``resolved_run.yaml`` / ``resolved_run.json``
    The final resolved configuration -- not the input file -- plus its SHA256.

``run.json``
    The run record: identity, snapshot, configuration hash, git state,
    environment and tool versions, per-stage input/output counts and hashes,
    Slurm job identifiers, validation verdicts and the final release.

``events.jsonl``
    An append-only event log.  Nothing in it is ever rewritten, so retries and
    failures stay visible even after a later attempt succeeds.

Historical provenance is never overwritten: a run directory is created with
``os.makedirs(exist_ok=False)`` and refuses to reuse an existing run ID.
"""

from __future__ import annotations

import hashlib
import json
import os
import platform
import socket
import subprocess
import uuid
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

from pdbclean.runconfig import (
    ResolvedRunConfig,
    runtime_environment,
    write_resolved_config,
)


RUN_RECORD_SCHEMA_NAME = "pdbclean_run_provenance"
RUN_RECORD_SCHEMA_VERSION = "1.0"


class RunProvenanceError(RuntimeError):
    """Raised when run provenance cannot be created or updated."""


def utc_now_iso() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def new_run_id(*, resolved_config_sha256: str, when: datetime | None = None) -> str:
    """Return a unique, sortable, configuration-linked run identifier."""

    moment = when or datetime.now(timezone.utc)
    stamp = moment.strftime("%Y%m%dT%H%M%SZ")

    entropy = hashlib.sha256(
        f"{resolved_config_sha256}:{uuid.uuid4()}".encode("utf-8")
    ).hexdigest()[:8]

    return f"run-{stamp}-{entropy}"


def file_sha256(path: str | Path) -> str:
    digest = hashlib.sha256()

    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)

    return digest.hexdigest()


def _git(repo_root: Path, *args: str) -> str | None:
    try:
        completed = subprocess.run(
            ["git", *args],
            cwd=str(repo_root),
            capture_output=True,
            text=True,
            check=False,
            timeout=30,
        )
    except (OSError, subprocess.SubprocessError):
        return None

    if completed.returncode != 0:
        return None

    return completed.stdout.strip()


def collect_git_state(repo_root: str | Path) -> dict[str, Any]:
    """Record branch, commit and whether the working tree was dirty."""

    root = Path(repo_root)

    commit = _git(root, "rev-parse", "HEAD")
    branch = _git(root, "rev-parse", "--abbrev-ref", "HEAD")
    status = _git(root, "status", "--porcelain")

    dirty: bool | None

    if status is None:
        dirty = None
    else:
        dirty = bool(status.strip())

    return {
        "repository_root": str(root),
        "branch": branch,
        "commit": commit,
        "working_tree_dirty": dirty,
        "dirty_path_count": (
            len([line for line in status.splitlines() if line.strip()])
            if status is not None
            else None
        ),
    }


def _module_version(module_name: str) -> str | None:
    try:
        module = __import__(module_name)
    except Exception:
        return None

    return getattr(module, "__version__", None)


def collect_environment(*, repo_root: str | Path | None = None) -> dict[str, Any]:
    """Record the interpreter, host and scientific library versions."""

    payload: dict[str, Any] = {
        "hostname": socket.gethostname(),
        "platform": platform.platform(),
        "python_version": platform.python_version(),
        "python_executable": os.sys.executable,
        "numpy_version": _module_version("numpy"),
        "scipy_version": _module_version("scipy"),
        "pyarrow_version": _module_version("pyarrow"),
        "gemmi_version": _module_version("gemmi"),
        "yaml_version": _module_version("yaml"),
    }

    slurm = {
        key: os.environ[key]
        for key in (
            "SLURM_JOB_ID",
            "SLURM_ARRAY_JOB_ID",
            "SLURM_ARRAY_TASK_ID",
            "SLURM_JOB_PARTITION",
            "SLURM_JOB_NODELIST",
        )
        if key in os.environ
    }

    if slurm:
        payload["slurm"] = slurm

    if repo_root is not None:
        bri_version = Path(repo_root) / "reproducibility" / "bri_version.txt"

        if bri_version.is_file():
            payload["bri_implementation"] = (
                bri_version.read_text(encoding="utf-8").strip()
            )

    return payload


def _atomic_write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)

    tmp = path.with_name(f".{path.name}.tmp.{os.getpid()}")

    try:
        with tmp.open("w", encoding="utf-8") as handle:
            handle.write(text)
            handle.flush()
            os.fsync(handle.fileno())

        os.replace(tmp, path)
    finally:
        if tmp.exists():
            tmp.unlink()


@dataclass
class StageProvenance:
    """Provenance for one pipeline stage within one run."""

    stage_id: str
    title: str
    layer: str

    #: The authoritative scientific stage this entry realises ("Stage 7",
    #: "Stage 14a", ...). Provenance records the canonical identity, not the
    #: registry's execution-order index, so a re-ordering of the runner can
    #: never renumber the science in the historical record.
    canonical_stage: str = "prerequisite"

    status: str = "pending"
    validation: str = "pending"
    started_at: str | None = None
    finished_at: str | None = None
    runtime_seconds: float | None = None
    input_count: int | None = None
    output_count: int | None = None
    output_path: str | None = None
    manifest_path: str | None = None
    summary_path: str | None = None
    checksums: dict[str, str] = field(default_factory=dict)
    slurm_job_ids: list[str] = field(default_factory=list)
    attempts: int = 0
    reused: bool = False
    scientific_parameters: dict[str, Any] = field(default_factory=dict)
    messages: list[str] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return {
            "stage_id": self.stage_id,
            "title": self.title,
            "layer": self.layer,
            "canonical_stage": self.canonical_stage,
            "status": self.status,
            "validation": self.validation,
            "started_at": self.started_at,
            "finished_at": self.finished_at,
            "runtime_seconds": self.runtime_seconds,
            "input_count": self.input_count,
            "output_count": self.output_count,
            "output_path": self.output_path,
            "manifest_path": self.manifest_path,
            "summary_path": self.summary_path,
            "checksums": dict(self.checksums),
            "slurm_job_ids": list(self.slurm_job_ids),
            "attempts": self.attempts,
            "reused": self.reused,
            "scientific_parameters": dict(self.scientific_parameters),
            "messages": list(self.messages),
        }


class RunProvenance:
    """Owns one run directory and its provenance artefacts."""

    def __init__(
        self,
        *,
        run_id: str,
        run_dir: Path,
        record: dict[str, Any],
    ) -> None:
        self.run_id = run_id
        self.run_dir = run_dir
        self._record = record
        self._stages: dict[str, StageProvenance] = {}

    # -- construction ---------------------------------------------------

    @classmethod
    def create(
        cls,
        *,
        resolved: ResolvedRunConfig,
        run_root: str | Path,
        repo_root: str | Path,
        snapshot: dict[str, Any] | None = None,
        run_id: str | None = None,
        invocation: dict[str, Any] | None = None,
    ) -> "RunProvenance":
        """Create a new run directory and write provenance *before* work."""

        identifier = run_id or new_run_id(
            resolved_config_sha256=resolved.sha256
        )

        run_dir = Path(run_root) / identifier

        if run_dir.exists():
            raise RunProvenanceError(
                f"Run directory already exists, refusing to overwrite "
                f"historical provenance: {run_dir}"
            )

        run_dir.mkdir(parents=True, exist_ok=False)

        config_paths = write_resolved_config(resolved, run_dir)

        record: dict[str, Any] = {
            "schema_name": RUN_RECORD_SCHEMA_NAME,
            "schema_version": RUN_RECORD_SCHEMA_VERSION,
            "run_id": identifier,
            "created_at": utc_now_iso(),
            "status": "created",
            "run_directory": str(run_dir),
            "snapshot": snapshot or {},
            "resolved_config": resolved.to_dict(),
            "resolved_config_sha256": resolved.sha256,
            "scientific_config_sha256": resolved.scientific_sha256,
            "resolved_config_yaml": config_paths["resolved_config_yaml"],
            "resolved_config_json": config_paths["resolved_config_json"],
            "config_layers": list(resolved.layers),
            "config_file": resolved.config_path,
            "config_overrides": list(resolved.override_items),
            "config_value_sources": dict(resolved.sources),
            "defaults_version": resolved.get("defaults_version"),
            "git": collect_git_state(repo_root),
            "environment": collect_environment(repo_root=repo_root),
            # Host-specific values are execution provenance, never part of the
            # canonical configuration. Recording them here is what lets the
            # canonical configuration stay identical on every host.
            "runtime": {
                "created_on": runtime_environment(resolved),
                "observed": [],
            },
            "invocation": invocation or {},
            "inputs": {},
            "stages": [],
            "validation": {},
            "release": {},
        }

        instance = cls(run_id=identifier, run_dir=run_dir, record=record)
        instance.flush()
        instance.append_event(
            "run_created",
            resolved_config_sha256=resolved.sha256,
            scientific_config_sha256=resolved.scientific_sha256,
            snapshot_id=(snapshot or {}).get("snapshot_id"),
        )

        return instance

    @classmethod
    def load(cls, run_dir: str | Path) -> "RunProvenance":
        directory = Path(run_dir)
        record_path = directory / "run.json"

        if not record_path.is_file():
            raise RunProvenanceError(f"No run.json in {directory}")

        record = json.loads(record_path.read_text(encoding="utf-8"))

        instance = cls(
            run_id=record["run_id"],
            run_dir=directory,
            record=record,
        )

        for stage in record.get("stages", []):
            instance._stages[stage["stage_id"]] = StageProvenance(
                stage_id=stage["stage_id"],
                title=stage.get("title", stage["stage_id"]),
                layer=stage.get("layer", "gold"),
                canonical_stage=stage.get("canonical_stage", "prerequisite"),
                status=stage.get("status", "pending"),
                validation=stage.get("validation", "pending"),
                started_at=stage.get("started_at"),
                finished_at=stage.get("finished_at"),
                runtime_seconds=stage.get("runtime_seconds"),
                input_count=stage.get("input_count"),
                output_count=stage.get("output_count"),
                output_path=stage.get("output_path"),
                manifest_path=stage.get("manifest_path"),
                summary_path=stage.get("summary_path"),
                checksums=dict(stage.get("checksums", {})),
                slurm_job_ids=list(stage.get("slurm_job_ids", [])),
                attempts=stage.get("attempts", 0),
                reused=stage.get("reused", False),
                scientific_parameters=dict(
                    stage.get("scientific_parameters", {})
                ),
                messages=list(stage.get("messages", [])),
            )

        return instance

    # -- record access --------------------------------------------------

    @property
    def record(self) -> dict[str, Any]:
        return self._record

    @property
    def record_path(self) -> Path:
        return self.run_dir / "run.json"

    @property
    def events_path(self) -> Path:
        return self.run_dir / "events.jsonl"

    def stage(self, stage_id: str) -> StageProvenance | None:
        return self._stages.get(stage_id)

    # -- mutation -------------------------------------------------------

    def register_stage(
        self,
        *,
        stage_id: str,
        title: str,
        layer: str,
        canonical_stage: str = "prerequisite",
        scientific_parameters: dict[str, Any] | None = None,
    ) -> StageProvenance:
        stage = self._stages.get(stage_id)

        if stage is None:
            stage = StageProvenance(
                stage_id=stage_id,
                title=title,
                layer=layer,
                canonical_stage=canonical_stage,
                scientific_parameters=scientific_parameters or {},
            )
            self._stages[stage_id] = stage
        elif scientific_parameters:
            stage.scientific_parameters.update(scientific_parameters)

        return stage

    def update_stage(self, stage_id: str, **fields: Any) -> StageProvenance:
        stage = self._stages.get(stage_id)

        if stage is None:
            raise RunProvenanceError(f"Unknown stage: {stage_id}")

        for key, value in fields.items():
            if not hasattr(stage, key):
                raise RunProvenanceError(
                    f"Unknown stage provenance field: {key}"
                )

            setattr(stage, key, value)

        self.append_event("stage_updated", stage_id=stage_id, **fields)

        return stage

    @property
    def resolved_config(self) -> ResolvedRunConfig:
        """Return the canonical configuration this run was created with.

        Loaded from the run's own persisted document, never re-resolved, so it
        is identical on every host that touches the run.
        """

        from pdbclean.runconfig import load_resolved_config

        return load_resolved_config(self.run_dir)

    def record_runtime(
        self,
        resolved: ResolvedRunConfig,
        *,
        stage_id: str | None = None,
    ) -> dict[str, Any]:
        """Record the host-specific values one execution actually ran with.

        This is where ``$TMPDIR``, the node name and the Slurm identifiers
        live.  They are execution facts, not configuration: the canonical
        resolved configuration stays byte-identical whichever host runs it.
        """

        observation = runtime_environment(resolved)
        observation["observed_at"] = utc_now_iso()

        if stage_id is not None:
            observation["stage_id"] = stage_id

        runtime = self._record.setdefault(
            "runtime", {"created_on": {}, "observed": []}
        )
        runtime.setdefault("observed", []).append(observation)

        self.append_event(
            "runtime_environment",
            stage_id=stage_id,
            hostname=observation["hostname"],
            storage=observation["storage"],
        )

        return observation

    def record_input(self, name: str, **fields: Any) -> None:
        self._record.setdefault("inputs", {})[name] = fields

    def record_validation(self, name: str, verdict: str, **fields: Any) -> None:
        self._record.setdefault("validation", {})[name] = {
            "verdict": verdict,
            "recorded_at": utc_now_iso(),
            **fields,
        }
        self.append_event("validation", name=name, verdict=verdict, **fields)

    def record_release(
        self,
        *,
        release_path: str,
        artefacts: Iterable[dict[str, Any]] | None = None,
        **fields: Any,
    ) -> None:
        self._record["release"] = {
            "release_path": release_path,
            "recorded_at": utc_now_iso(),
            "artefacts": list(artefacts or []),
            **fields,
        }
        self.append_event("release", release_path=release_path, **fields)

    def set_status(self, status: str, **fields: Any) -> None:
        self._record["status"] = status
        self._record["status_updated_at"] = utc_now_iso()

        for key, value in fields.items():
            self._record[key] = value

        self.append_event("run_status", status=status, **fields)

    # -- persistence ----------------------------------------------------

    def append_event(self, event: str, **fields: Any) -> None:
        """Append one immutable line to the run event log."""

        payload = {
            "at": utc_now_iso(),
            "event": event,
            **{
                key: value
                for key, value in fields.items()
                if _json_safe(value)
            },
        }

        self.run_dir.mkdir(parents=True, exist_ok=True)

        with self.events_path.open("a", encoding="utf-8") as handle:
            handle.write(json.dumps(payload, sort_keys=True) + "\n")

    def flush(self) -> Path:
        """Write the current run record atomically."""

        self._record["stages"] = [
            stage.to_dict()
            for stage in sorted(
                self._stages.values(),
                key=lambda item: item.stage_id,
            )
        ]
        self._record["updated_at"] = utc_now_iso()

        _atomic_write_text(
            self.record_path,
            json.dumps(self._record, indent=2, sort_keys=True) + "\n",
        )

        return self.record_path


def _json_safe(value: Any) -> bool:
    try:
        json.dumps(value)
    except (TypeError, ValueError):
        return False

    return True


def list_runs(run_root: str | Path) -> list[dict[str, Any]]:
    """Return a newest-first summary of every recorded run."""

    root = Path(run_root)

    if not root.is_dir():
        return []

    summaries: list[dict[str, Any]] = []

    for entry in sorted(root.iterdir(), reverse=True):
        record_path = entry / "run.json"

        if not record_path.is_file():
            continue

        try:
            record = json.loads(record_path.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            continue

        summaries.append(
            {
                "run_id": record.get("run_id", entry.name),
                "created_at": record.get("created_at"),
                "updated_at": record.get("updated_at"),
                "status": record.get("status"),
                "snapshot_id": (record.get("snapshot") or {}).get(
                    "snapshot_id"
                ),
                "resolved_config_sha256": record.get(
                    "resolved_config_sha256"
                ),
                "scientific_config_sha256": record.get(
                    "scientific_config_sha256"
                ),
                "run_directory": str(entry),
                "stage_count": len(record.get("stages", [])),
            }
        )

    return summaries
