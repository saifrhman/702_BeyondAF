"""Read-only inspection of a historical run.

Everything here opens a completed run's own records and artefacts and returns
what they say.  Nothing in this module writes: inspecting a run must never
modify ``run.json``, append an event, re-resolve a snapshot, recompute a hash
or launch anything.

The timeline it returns is always in **canonical scientific order** -- the
order of :func:`pdbclean.stage_registry.canonical_timeline` -- never
alphabetical and never the orchestrator's execution index.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from pdbclean import artefacts
from pdbclean.stage_registry import (
    method_references,
    ROLE_INVESTIGATION,
    ROLE_VALIDATION,
    STAGES_BY_ID,
    canonical_timeline,
)


#: Status shown for a canonical identity that the orchestrated pipeline does
#: not produce (Stages 11-13).  These are real scientific stages, so they are
#: displayed rather than omitted -- with an honest status.
STATUS_NOT_EXECUTED = "not_executed"
STATUS_NOT_RECORDED = "not_recorded"
STATUS_OPTIONAL = "optional_investigation"

#: Shown wherever a historical run simply did not record a field.  Never
#: invented, never back-filled.
NOT_RECORDED = "not recorded"


def _stage_records(record: dict[str, Any]) -> dict[str, dict[str, Any]]:
    return {
        entry["stage_id"]: entry
        for entry in record.get("stages", [])
        if isinstance(entry, dict) and "stage_id" in entry
    }


def _runtime_for_stage(
    record: dict[str, Any],
    stage_id: str,
) -> list[dict[str, Any]]:
    runtime = record.get("runtime") or {}

    return [
        observation
        for observation in runtime.get("observed", [])
        if isinstance(observation, dict)
        and observation.get("stage_id") == stage_id
    ]


def run_timeline(record: dict[str, Any]) -> list[dict[str, Any]]:
    """Return the canonical timeline for one run, in canonical order.

    One row per canonical scientific identity.  Where two identities share a
    producer (Stage 3/4, Stage 8/9) both rows appear and both reference the
    same producer, which is stated explicitly rather than hidden.
    """

    stages = _stage_records(record)

    rows: list[dict[str, Any]] = []

    for entry in canonical_timeline():
        producer = entry.producer
        recorded = stages.get(producer) if producer else None

        if recorded is not None:
            status = recorded.get("status") or STATUS_NOT_RECORDED
            validation = recorded.get("validation") or STATUS_NOT_RECORDED
            reused = bool(recorded.get("reused"))
        elif producer is None:
            # A canonical stage with no producer on the release path.
            status = (
                STATUS_OPTIONAL
                if entry.role in {ROLE_INVESTIGATION, ROLE_VALIDATION}
                else STATUS_NOT_EXECUTED
            )
            validation = "not_applicable"
            reused = False
        else:
            status = STATUS_NOT_RECORDED
            validation = STATUS_NOT_RECORDED
            reused = False

        spec = STAGES_BY_ID.get(producer) if producer else None

        rows.append(
            {
                "key": entry.key,
                "label": entry.label,
                "title": entry.title,
                "display": entry.display,
                "position": entry.position,
                "role": entry.role,
                "parent": entry.parent,
                "layer": entry.layer,
                "purpose": entry.purpose,
                "note": entry.note,
                "producer": producer,
                "shared_producer": entry.shared_producer,
                "frozen_output": entry.frozen_output,
                "status": status,
                "validation": validation,
                "reused": reused,
                "input_count": (
                    recorded.get("input_count") if recorded else None
                ),
                "output_count": (
                    recorded.get("output_count") if recorded else None
                ),
                "slurm_job_ids": (
                    list(recorded.get("slurm_job_ids", [])) if recorded else []
                ),
                "entry_point": spec.entry_point if spec else None,
                "has_detail": recorded is not None,
            }
        )

    return rows


def stage_detail(
    record: dict[str, Any],
    canonical_key: str,
    *,
    repo_root: str | Path | None = None,
    list_artefacts: bool = True,
) -> dict[str, Any]:
    """Return everything one run recorded about one canonical stage.

    Read-only.  Fields the run did not record are reported as
    :data:`NOT_RECORDED` rather than invented.
    """

    entry = next(
        (item for item in canonical_timeline() if item.key == canonical_key),
        None,
    )

    if entry is None:
        raise KeyError(f"Unknown canonical stage: {canonical_key}")

    stages = _stage_records(record)
    producer = entry.producer
    recorded = stages.get(producer) if producer else None
    spec = STAGES_BY_ID.get(producer) if producer else None

    identity = {
        "canonical_key": entry.key,
        "canonical_label": entry.label,
        "canonical_title": entry.title,
        "display": entry.display,
        "position": entry.position,
        "role": entry.role,
        "parent_stage": entry.parent or NOT_RECORDED,
        "layer": entry.layer,
        "purpose": entry.purpose,
        "note": entry.note or NOT_RECORDED,
        "producer": producer or NOT_RECORDED,
        "shared_producer": entry.shared_producer,
        "sibling_identities": [
            item.display
            for item in canonical_timeline()
            if producer is not None
            and item.producer == producer
            and item.key != entry.key
        ],
        "implementation": (spec.entry_point if spec else NOT_RECORDED),
        "frozen_output": entry.frozen_output or NOT_RECORDED,
    }

    # The scientific explanation is shared across runs; the values shown
    # alongside it come from this run's own record.
    description = {
        "rationale": entry.rationale or NOT_RECORDED,
        "scientific_method": entry.scientific_method or NOT_RECORDED,
        "stage_input": entry.stage_input or NOT_RECORDED,
        "stage_output": entry.stage_output or NOT_RECORDED,
        "downstream_role": entry.downstream_role or NOT_RECORDED,
        "implementation_note": entry.implementation_note or NOT_RECORDED,
        "references": method_references(entry.references),
    }

    if recorded is None:
        return {
            "identity": identity,
            "description": description,
            "status": {
                "status": (
                    STATUS_OPTIONAL
                    if entry.role in {ROLE_INVESTIGATION, ROLE_VALIDATION}
                    else STATUS_NOT_EXECUTED
                ),
                "validation": "not_applicable",
                "reused": False,
                "explanation": (
                    entry.note
                    or "This canonical stage was not produced by this run."
                ),
            },
            "configuration": {},
            "inputs": {},
            "outputs": {},
            "validation": {},
            "execution": {},
            "reuse": {},
            "artefacts": [],
        }

    def _or_not_recorded(value: Any) -> Any:
        return NOT_RECORDED if value in (None, "", [], {}) else value

    run_config = record.get("resolved_config") or {}
    snapshot = record.get("snapshot") or {}
    git = record.get("git") or {}
    environment = record.get("environment") or {}

    configuration = dict(recorded.get("scientific_parameters") or {})
    configuration.setdefault(
        "snapshot.snapshot_id", snapshot.get("snapshot_id", NOT_RECORDED)
    )

    # Stage-relevant scientific context, taken from the run's own frozen
    # configuration -- never re-resolved.
    for dotted in (
        "bri.representation_precision_angstrom",
        "brain_filter.threshold_angstrom",
        "duplicate_search.near_duplicate_threshold_angstrom",
        "selection.models.model_id",
    ):
        section, _, leaf = dotted.partition(".")
        cursor: Any = run_config.get(section)

        for part in leaf.split("."):
            cursor = cursor.get(part) if isinstance(cursor, dict) else None

        if cursor is not None:
            configuration.setdefault(dotted, cursor)

    upstream = [
        item.display
        for item in canonical_timeline()
        if spec is not None and item.producer in set(spec.depends_on)
    ]

    directory = recorded.get("output_path")

    artefact_rows: list[dict[str, Any]] = []

    if list_artefacts and directory:
        artefact_rows = artefacts.list_directory(directory)

    return {
        "identity": identity,
        "description": description,
        "status": {
            "status": _or_not_recorded(recorded.get("status")),
            "validation": _or_not_recorded(recorded.get("validation")),
            "reused": bool(recorded.get("reused")),
            "attempts": recorded.get("attempts", 0),
            "messages": recorded.get("messages") or [],
        },
        "configuration": configuration,
        "inputs": {
            "input_count": _or_not_recorded(recorded.get("input_count")),
            "upstream_stages": upstream or [NOT_RECORDED],
            "manifest_path": _or_not_recorded(recorded.get("manifest_path")),
            "run_inputs": record.get("inputs") or {},
        },
        "outputs": {
            "output_count": _or_not_recorded(recorded.get("output_count")),
            "output_path": _or_not_recorded(directory),
            "summary_path": _or_not_recorded(recorded.get("summary_path")),
            "checksums": recorded.get("checksums") or {},
        },
        "validation": {
            "verdict": _or_not_recorded(recorded.get("validation")),
            "gate": spec.validation if spec else NOT_RECORDED,
            "recorded": {
                name: verdict
                for name, verdict in (record.get("validation") or {}).items()
                if producer is not None and name.startswith(producer)
            },
        },
        "execution": {
            "started_at": _or_not_recorded(recorded.get("started_at")),
            "finished_at": _or_not_recorded(recorded.get("finished_at")),
            "runtime_seconds": _or_not_recorded(
                recorded.get("runtime_seconds")
            ),
            "slurm_job_ids": recorded.get("slurm_job_ids") or [],
            "entry_point": _or_not_recorded(
                spec.entry_point if spec else None
            ),
            "git_branch": _or_not_recorded(git.get("branch")),
            "git_commit": _or_not_recorded(git.get("commit")),
            "working_tree_dirty": git.get("working_tree_dirty"),
            "python_version": _or_not_recorded(
                environment.get("python_version")
            ),
            "bri_implementation": _or_not_recorded(
                environment.get("bri_implementation")
            ),
            "resolved_config_sha256": _or_not_recorded(
                record.get("resolved_config_sha256")
            ),
            "scientific_config_sha256": _or_not_recorded(
                record.get("scientific_config_sha256")
            ),
            "runtime_observations": _runtime_for_stage(record, producer),
        },
        "reuse": {
            "reused": bool(recorded.get("reused")),
            "explanation": (
                "Existing validated output was reused: every compatibility "
                "key recorded in the stage summary agreed with this run's "
                "resolved configuration."
                if recorded.get("reused")
                else "Output was produced by this run."
            ),
        },
        "artefacts": artefact_rows,
    }


#: Canonical stages whose outputs can be opened in the Duplicate Explorer.
DUPLICATE_STAGE_KEYS = frozenset(
    {"stage_8", "stage_9", "stage_10", "stage_14a", "stage_14b"}
)


def duplicate_navigation(canonical_key: str) -> dict[str, Any] | None:
    """Return the Duplicate Explorer entry point for a duplicate-bearing stage.

    Mol* inspection is reached from the Explorer.  Neither ever reclassifies a
    pair: the classification shown is the one the pipeline recorded.
    """

    if canonical_key not in DUPLICATE_STAGE_KEYS:
        return None

    filters = {
        "stage_8": {},
        "stage_9": {},
        "stage_10": {},
        "stage_14a": {"relationship": "removed"},
        "stage_14b": {"relationship": "removed"},
    }[canonical_key]

    return {
        "available": True,
        "label": "Open in Duplicate Explorer",
        "filters": filters,
        "caveat": (
            "The Duplicate Explorer filters and displays recorded results. "
            "Mol* is for human inspection only and never determines whether "
            "two chains are duplicates."
        ),
    }
