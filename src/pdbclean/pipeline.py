"""Top-level orchestration for the COMP702 PDBClean pipeline.

This module is an *orchestrator*, not a second scientific implementation.  It
inspects the artefacts the existing stage entry points produce, decides which
stages may be reused, enforces the validation gates, and submits or reports the
stages that still need to run.  No duplicate-detection, BRI, Brain, filtering or
representative-selection logic lives here.

Validation philosophy, preserved from the existing pipeline:

    stage runs -> output materialises -> validation runs -> PASS ->
    downstream stage may start

A process exiting zero is not scientific validation.  A stage is only
``validation_pass`` when its own success marker and machine-readable summary
agree with the resolved configuration.

Restartability is decided from configuration and provenance, never from
directory names: an existing output is only reused when every compatibility key
declared in :mod:`pdbclean.stage_registry` matches the resolved configuration.
"""

from __future__ import annotations

import json
import shlex
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable

from pdbclean.run_provenance import RunProvenance, file_sha256
from pdbclean.defaults import (
    IMPLEMENTED_REPRESENTATION_PRECISION_ANGSTROM,
    precision_is_implemented,
)
from pdbclean.runconfig import ResolvedRunConfig
from pdbclean.stage_registry import (
    PRECISION_DEPENDENT_STAGES,
    STAGES,
    LAYER_GOLD,
    StageSpec,
    release_name,
    resolve_directory,
)


# Stage lifecycle states.  The UI renders exactly these.
PENDING = "pending"
RUNNING = "running"
EXECUTION_COMPLETE = "execution_complete"
VALIDATING = "validating"
VALIDATION_PASS = "validation_pass"
VALIDATION_FAIL = "validation_fail"
PARTIAL = "partial"
COMPLETE = "complete"
BLOCKED = "blocked"
NOT_APPLICABLE = "not_applicable"

STATE_ORDER: tuple[str, ...] = (
    PENDING,
    BLOCKED,
    RUNNING,
    EXECUTION_COMPLETE,
    VALIDATING,
    PARTIAL,
    VALIDATION_FAIL,
    VALIDATION_PASS,
    COMPLETE,
    NOT_APPLICABLE,
)

#: Actions the planner can assign to a stage.
ACTION_REUSE = "reuse"
ACTION_RUN = "run"
ACTION_BLOCKED = "blocked"
ACTION_NOT_APPLICABLE = "not_applicable"

#: Outputs above this size are recorded without a checksum unless the caller
#: explicitly asks for one, so that inspecting a pipeline stays cheap enough to
#: run on a login node.
DEFAULT_MAX_CHECKSUM_BYTES = 256 * 1024 * 1024


class PipelineError(RuntimeError):
    """Raised when the pipeline cannot be planned or executed."""


@dataclass(frozen=True)
class PipelinePaths:
    """Concrete filesystem roots for one resolved run."""

    repo_root: Path
    output_root: Path
    release_root: Path
    run_root: Path
    snapshot: str
    protocol: str
    release: str

    @classmethod
    def from_config(
        cls,
        resolved: ResolvedRunConfig,
        *,
        repo_root: str | Path,
    ) -> "PipelinePaths":
        root = Path(repo_root)

        snapshot = resolved.get("snapshot.snapshot_id")

        if not snapshot:
            raise PipelineError(
                "Cannot derive pipeline paths before the snapshot has been "
                "resolved to a concrete identity"
            )

        protocol = resolved.get("release.protocol_version")
        dataset = resolved.get("release.dataset_name")

        def _root(dotted: str, fallback: str) -> Path:
            value = resolved.get(dotted) or fallback
            candidate = Path(value)

            return candidate if candidate.is_absolute() else root / candidate

        return cls(
            repo_root=root,
            output_root=_root("storage.output_root", "outputs/pdbclean"),
            release_root=_root("storage.release_root", "outputs/releases"),
            run_root=_root("storage.run_root", "outputs/runs"),
            snapshot=str(snapshot),
            protocol=str(protocol),
            release=release_name(
                dataset_name=str(dataset),
                snapshot=str(snapshot),
                protocol=str(protocol),
            ),
        )

    def stage_directory(self, stage: StageSpec) -> Path | None:
        template = resolve_directory(
            stage,
            output_root=str(self.output_root),
            snapshot=self.snapshot,
            protocol=self.protocol,
            release_root=str(self.release_root),
            release=self.release,
        )

        return None if template is None else Path(template)


@dataclass
class StageObservation:
    """What the filesystem says about one stage of one configuration."""

    stage: StageSpec
    status: str = PENDING
    validation: str = PENDING
    action: str = ACTION_RUN
    directory: str | None = None
    exists: bool = False
    success_marker_present: bool = False
    summary_path: str | None = None
    summary: dict[str, Any] = field(default_factory=dict)
    primary_output_path: str | None = None
    primary_output_bytes: int | None = None
    primary_output_sha256: str | None = None
    input_count: int | None = None
    output_count: int | None = None
    incompatibilities: list[dict[str, Any]] = field(default_factory=list)
    messages: list[str] = field(default_factory=list)
    scientific_parameters: dict[str, Any] = field(default_factory=dict)
    slurm_job_ids: list[str] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return {
            "stage_id": self.stage.stage_id,
            "ordinal": self.stage.ordinal,
            "canonical_stage": self.stage.canonical_stage,
            "title": self.stage.title,
            "layer": self.stage.layer,
            "purpose": self.stage.purpose,
            "persisted": self.stage.persisted,
            "depends_on": list(self.stage.depends_on),
            "validation_description": self.stage.validation,
            "entry_point": self.stage.entry_point,
            "status": self.status,
            "validation": self.validation,
            "action": self.action,
            "directory": self.directory,
            "exists": self.exists,
            "success_marker_present": self.success_marker_present,
            "summary_path": self.summary_path,
            "summary": self.summary,
            "manifest_path": self.primary_output_path,
            "output_path": self.directory,
            "primary_output_bytes": self.primary_output_bytes,
            "primary_output_sha256": self.primary_output_sha256,
            "input_count": self.input_count,
            "output_count": self.output_count,
            "incompatibilities": list(self.incompatibilities),
            "messages": list(self.messages),
            "scientific_parameters": dict(self.scientific_parameters),
            "slurm_job_ids": list(self.slurm_job_ids),
        }


def _read_json(path: Path) -> dict[str, Any] | None:
    try:
        loaded = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None

    return loaded if isinstance(loaded, dict) else None


def _first_count(summary: dict[str, Any], keys: Iterable[str]) -> int | None:
    for key in keys:
        value = summary.get(key)

        if isinstance(value, bool):
            continue

        if isinstance(value, int):
            return value

    return None


def _expected_compatibility_value(
    key: str,
    resolved: ResolvedRunConfig,
) -> Any:
    """Resolve a compatibility expectation, including pseudo-keys."""

    if key == "__near_duplicate_threshold_text__":
        return f"d_bri_mA <= {resolved.near_duplicate_threshold_mA}"

    return resolved.get(key)


def _values_agree(observed: Any, expected: Any) -> bool:
    if observed is None or expected is None:
        return observed == expected

    if isinstance(observed, (int, float)) and isinstance(
        expected, (int, float)
    ):
        if isinstance(observed, bool) or isinstance(expected, bool):
            return observed == expected

        return abs(float(observed) - float(expected)) <= 1.0e-12

    return str(observed) == str(expected)


def inspect_stage(
    stage: StageSpec,
    *,
    paths: PipelinePaths,
    resolved: ResolvedRunConfig,
    compute_checksums: bool = False,
    max_checksum_bytes: int = DEFAULT_MAX_CHECKSUM_BYTES,
) -> StageObservation:
    """Inspect one stage's artefacts against the resolved configuration."""

    observation = StageObservation(stage=stage)

    observation.scientific_parameters = {
        key: resolved.get(key) for key in stage.scientific_parameters
    }

    if stage.stage_id == "snapshot":
        # The snapshot stage has no artefact of its own: it is satisfied when
        # the run configuration carries a concrete, pinned snapshot identity.
        snapshot_id = resolved.get("snapshot.snapshot_id")

        if snapshot_id:
            observation.status = COMPLETE
            observation.validation = VALIDATION_PASS
            observation.action = ACTION_REUSE
            observation.output_count = None
            observation.messages.append(
                f"Snapshot pinned to {snapshot_id} "
                f"(selection mode: "
                f"{resolved.get('snapshot.resolved_selection_mode') or resolved.get('snapshot.mode')})."
            )
        else:
            observation.status = PENDING
            observation.validation = PENDING
            observation.action = ACTION_RUN
            observation.messages.append(
                "Snapshot has not yet been resolved to a concrete identity."
            )

        return observation

    if not stage.persisted:
        observation.status = NOT_APPLICABLE
        observation.validation = NOT_APPLICABLE
        observation.action = ACTION_NOT_APPLICABLE
        observation.messages.append(
            "Deterministic stage with no persisted artefact by design."
        )
        return observation

    directory = paths.stage_directory(stage)

    if directory is None:
        observation.status = NOT_APPLICABLE
        observation.validation = NOT_APPLICABLE
        observation.action = ACTION_NOT_APPLICABLE
        return observation

    observation.directory = str(directory)
    observation.exists = directory.is_dir()

    if not observation.exists:
        observation.status = PENDING
        observation.validation = PENDING
        observation.action = ACTION_RUN
        return observation

    if stage.success_marker:
        marker = directory / stage.success_marker
        observation.success_marker_present = marker.exists()
    else:
        observation.success_marker_present = True

    summary: dict[str, Any] = {}

    if stage.summary:
        summary_path = directory / stage.summary
        observation.summary_path = str(summary_path)

        loaded = _read_json(summary_path)

        if loaded is not None:
            summary = loaded

    observation.summary = summary

    if stage.primary_output:
        output_path = directory / stage.primary_output
        observation.primary_output_path = str(output_path)

        if output_path.is_file():
            observation.primary_output_bytes = output_path.stat().st_size

            if compute_checksums and (
                observation.primary_output_bytes <= max_checksum_bytes
            ):
                observation.primary_output_sha256 = file_sha256(output_path)
        else:
            observation.messages.append(
                f"Primary output missing: {stage.primary_output}"
            )

    observation.input_count = _first_count(summary, stage.input_count_keys)
    observation.output_count = _first_count(summary, stage.output_count_keys)

    for job_key in ("slurm_job_id", "analysis_job_id", "job_id"):
        value = summary.get(job_key)

        if value is not None:
            observation.slurm_job_ids.append(str(value))

    # Representation precision is not recorded in the frozen stage summaries,
    # because no grid other than the implemented one has ever been executable.
    # An existing artefact was therefore necessarily produced on that grid, so
    # a run configured for a different precision must never reuse it.
    if stage.stage_id in PRECISION_DEPENDENT_STAGES:
        try:
            implemented = precision_is_implemented(resolved.data)
        except ValueError:
            implemented = True

        if not implemented:
            observation.incompatibilities.append(
                {
                    "summary_key": "bri_representation_precision_angstrom",
                    "config_key": "bri.representation_precision_angstrom",
                    "observed": IMPLEMENTED_REPRESENTATION_PRECISION_ANGSTROM,
                    "expected": resolved.get(
                        "bri.representation_precision_angstrom"
                    ),
                }
            )

    # Compatibility: configuration and provenance, never directory names.
    for summary_key, config_key in stage.compatibility.items():
        expected = _expected_compatibility_value(config_key, resolved)
        observed = summary.get(summary_key)

        if observed is None:
            observation.messages.append(
                f"Summary does not record {summary_key!r}; compatibility with "
                f"{config_key} could not be confirmed."
            )
            continue

        if not _values_agree(observed, expected):
            observation.incompatibilities.append(
                {
                    "summary_key": summary_key,
                    "config_key": config_key,
                    "observed": observed,
                    "expected": expected,
                }
            )

    has_output = (
        observation.primary_output_path is None
        or Path(observation.primary_output_path).is_file()
    )

    if observation.incompatibilities:
        observation.status = EXECUTION_COMPLETE
        observation.validation = VALIDATION_FAIL
        observation.action = ACTION_RUN
        observation.messages.append(
            "Existing output was produced under a different scientific "
            "configuration and will not be reused."
        )
        return observation

    if observation.success_marker_present and has_output:
        observation.status = COMPLETE
        observation.validation = VALIDATION_PASS
        observation.action = ACTION_REUSE
        return observation

    if has_output and not observation.success_marker_present:
        observation.status = PARTIAL
        observation.validation = PENDING
        observation.action = ACTION_RUN
        observation.messages.append(
            "Output present without a success marker; the stage did not "
            "complete its own validation."
        )
        return observation

    observation.status = PARTIAL
    observation.validation = PENDING
    observation.action = ACTION_RUN

    return observation


@dataclass
class PipelinePlan:
    """The full stage plan for one resolved configuration."""

    resolved: ResolvedRunConfig
    paths: PipelinePaths
    observations: list[StageObservation]

    @property
    def by_id(self) -> dict[str, StageObservation]:
        return {obs.stage.stage_id: obs for obs in self.observations}

    @property
    def blocked(self) -> list[StageObservation]:
        return [obs for obs in self.observations if obs.action == ACTION_BLOCKED]

    @property
    def to_run(self) -> list[StageObservation]:
        return [obs for obs in self.observations if obs.action == ACTION_RUN]

    @property
    def reusable(self) -> list[StageObservation]:
        return [obs for obs in self.observations if obs.action == ACTION_REUSE]

    @property
    def complete(self) -> bool:
        return all(
            obs.action in {ACTION_REUSE, ACTION_NOT_APPLICABLE}
            for obs in self.observations
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "snapshot": self.resolved.get("snapshot.snapshot_id"),
            "protocol": self.resolved.get("release.protocol_version"),
            "resolved_config_sha256": self.resolved.sha256,
            "scientific_config_sha256": self.resolved.scientific_sha256,
            "release": self.paths.release,
            "release_root": str(self.paths.release_root),
            "output_root": str(self.paths.output_root),
            "complete": self.complete,
            "stages": [obs.to_dict() for obs in self.observations],
        }


def plan_pipeline(
    resolved: ResolvedRunConfig,
    *,
    repo_root: str | Path,
    compute_checksums: bool = False,
    stages: Iterable[StageSpec] | None = None,
) -> PipelinePlan:
    """Inspect every stage and decide reuse / run / blocked."""

    paths = PipelinePaths.from_config(resolved, repo_root=repo_root)

    specs = sorted(
        stages if stages is not None else STAGES,
        key=lambda item: item.ordinal,
    )

    observations: list[StageObservation] = []
    by_id: dict[str, StageObservation] = {}

    for spec in specs:
        observation = inspect_stage(
            spec,
            paths=paths,
            resolved=resolved,
            compute_checksums=compute_checksums,
        )

        # A stage may not start until every dependency has passed validation.
        blocking = [
            dependency
            for dependency in spec.depends_on
            if dependency in by_id
            and by_id[dependency].validation
            not in {VALIDATION_PASS, NOT_APPLICABLE}
        ]

        if blocking and observation.action == ACTION_RUN:
            observation.action = ACTION_BLOCKED
            observation.status = BLOCKED
            observation.messages.append(
                "Blocked: upstream validation has not passed for "
                + ", ".join(blocking)
            )

        observations.append(observation)
        by_id[spec.stage_id] = observation

    return PipelinePlan(
        resolved=resolved,
        paths=paths,
        observations=observations,
    )


def record_plan_in_provenance(
    plan: PipelinePlan,
    provenance: RunProvenance,
) -> None:
    """Copy the plan into the run provenance record."""

    for observation in plan.observations:
        stage = observation.stage

        provenance.register_stage(
            stage_id=stage.stage_id,
            title=stage.title,
            layer=stage.layer,
            canonical_stage=stage.canonical_stage,
            scientific_parameters=observation.scientific_parameters,
        )

        provenance.update_stage(
            stage.stage_id,
            status=observation.status,
            validation=observation.validation,
            input_count=observation.input_count,
            output_count=observation.output_count,
            output_path=observation.directory,
            manifest_path=observation.primary_output_path,
            summary_path=observation.summary_path,
            reused=(observation.action == ACTION_REUSE),
            slurm_job_ids=list(observation.slurm_job_ids),
            messages=list(observation.messages),
            checksums=(
                {"primary_output": observation.primary_output_sha256}
                if observation.primary_output_sha256
                else {}
            ),
        )

        if observation.incompatibilities:
            provenance.record_validation(
                f"{stage.stage_id}:configuration_compatibility",
                "FAIL",
                incompatibilities=observation.incompatibilities,
            )
        elif observation.validation == VALIDATION_PASS:
            provenance.record_validation(
                f"{stage.stage_id}:stage_gate",
                "PASS",
                output_count=observation.output_count,
            )

    provenance.flush()


# ----------------------------------------------------------------------
# Executors
# ----------------------------------------------------------------------


@dataclass
class ExecutionResult:
    stage_id: str
    command: list[str]
    executed: bool
    returncode: int | None = None
    stdout: str = ""
    stderr: str = ""
    slurm_job_id: str | None = None

    def to_dict(self) -> dict[str, Any]:
        return {
            "stage_id": self.stage_id,
            "command": list(self.command),
            "command_text": " ".join(shlex.quote(part) for part in self.command),
            "executed": self.executed,
            "returncode": self.returncode,
            "slurm_job_id": self.slurm_job_id,
            "stdout_tail": self.stdout[-4000:],
            "stderr_tail": self.stderr[-4000:],
        }


class Executor:
    """Base executor. Subclasses decide how a stage command is realised."""

    name = "base"

    def run(self, stage_id: str, command: list[str], *, cwd: Path) -> ExecutionResult:
        raise NotImplementedError


class DryRunExecutor(Executor):
    """Print what would run. The default, and always safe."""

    name = "dry-run"

    def run(self, stage_id: str, command: list[str], *, cwd: Path) -> ExecutionResult:
        return ExecutionResult(
            stage_id=stage_id,
            command=command,
            executed=False,
        )


class LocalExecutor(Executor):
    """Run a stage in the current process' machine.

    Intended for small stages and for fixtures.  Heavy stages should use the
    Slurm executor; Barkla login nodes are for orchestration only.
    """

    name = "local"

    def __init__(self, *, timeout: int | None = None) -> None:
        self.timeout = timeout

    def run(self, stage_id: str, command: list[str], *, cwd: Path) -> ExecutionResult:
        completed = subprocess.run(
            command,
            cwd=str(cwd),
            capture_output=True,
            text=True,
            check=False,
            timeout=self.timeout,
        )

        return ExecutionResult(
            stage_id=stage_id,
            command=command,
            executed=True,
            returncode=completed.returncode,
            stdout=completed.stdout,
            stderr=completed.stderr,
        )


class SlurmExecutor(Executor):
    """Submit a stage with ``sbatch --parsable``.

    Batch scripts are always submitted, never executed directly.
    """

    name = "slurm"

    def __init__(self, *, extra_sbatch_args: Iterable[str] = ()) -> None:
        self.extra_sbatch_args = list(extra_sbatch_args)

    def run(self, stage_id: str, command: list[str], *, cwd: Path) -> ExecutionResult:
        submission = ["sbatch", "--parsable", *self.extra_sbatch_args, *command]

        completed = subprocess.run(
            submission,
            cwd=str(cwd),
            capture_output=True,
            text=True,
            check=False,
        )

        job_id: str | None = None

        if completed.returncode == 0:
            first = completed.stdout.strip().splitlines()
            if first:
                job_id = first[0].split(";", 1)[0].strip() or None

        return ExecutionResult(
            stage_id=stage_id,
            command=submission,
            executed=True,
            returncode=completed.returncode,
            stdout=completed.stdout,
            stderr=completed.stderr,
            slurm_job_id=job_id,
        )


def build_executor(name: str, **kwargs: Any) -> Executor:
    if name in {"dry-run", "dry_run", "plan"}:
        return DryRunExecutor()

    if name == "local":
        return LocalExecutor(**kwargs)

    if name == "slurm":
        return SlurmExecutor(**kwargs)

    raise PipelineError(
        f"Unknown executor {name!r}; expected dry-run, local or slurm"
    )


def gold_release_summary(plan: PipelinePlan) -> dict[str, Any]:
    """Return the Gold release page payload for a plan, or an empty dict.

    Nothing is fabricated: if the release has not been published for this
    configuration, the caller gets ``{}`` and the UI shows an explicit
    "not yet produced" state rather than placeholder numbers.
    """

    observation = plan.by_id.get("gold_release")

    if observation is None or observation.validation != VALIDATION_PASS:
        return {}

    manifest = observation.summary

    if not manifest:
        return {}

    payload = dict(manifest)
    payload["release_directory"] = observation.directory
    payload["retained_manifest"] = observation.primary_output_path

    return payload


def gold_stage_observations(plan: PipelinePlan) -> list[StageObservation]:
    return [obs for obs in plan.observations if obs.stage.layer == LAYER_GOLD]
