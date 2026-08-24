"""The one configuration document a frozen run's workers execute.

A run resolves its configuration once, from validated defaults + profile +
explicit overrides, and freezes it as ``resolved_run.json``.  That document is
canonical and hashed, but the production stage entry points cannot read it:
they consume the Protocol 3.2 schema validated by
:func:`pdbclean.config.load_config`.

This module closes that gap.  It projects the frozen canonical configuration
onto the stage-execution schema, writes it into the run directory as
``stage_config.yaml``, and provides the verification a worker performs before
any scientific work begins.

Why it exists
-------------

Before this module, every generated stage command pointed at
``config/pdbclean/protocol_3_2_comp702_v1.yaml`` -- the byte-frozen *base*
file.  That file has no ``brain_filter`` section, so a UI-entered Brain
threshold of 0.005 A was silently replaced by the module constant
``brain_prefilter_production.BRAIN_PREFILTER_TAU_ANGSTROM`` (0.010 A).  It also
pins ``duplicate_search.near_duplicate_threshold_angstrom: 0.010`` and
``snapshot.mode: latest_complete``, so an entered tau and an entered snapshot
were both discarded too.  Nothing failed; the run simply computed different
science from the one the operator reviewed and froze.

The document this module writes is:

*deterministic*
    a pure function of the canonical configuration, so the same frozen run
    always projects to the same bytes and a tampered copy is detectable;

*complete*
    it carries every scientific section, so no stage can fall back to a
    source-code default;

*self-identifying*
    it embeds the canonical configuration's two SHA256 digests, so a worker
    holding only this file can still prove which run it belongs to.
"""

from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Any, Mapping

import yaml

from pdbclean.runconfig import ResolvedRunConfig, RunConfigError


#: Filename of the executable projection inside a run directory.
STAGE_CONFIG_BASENAME = "stage_config.yaml"

STAGE_CONFIG_SCHEMA_NAME = "pdbclean_stage_config"
STAGE_CONFIG_SCHEMA_VERSION = "1.0"


class StageConfigError(RuntimeError):
    """Raised when the executable configuration is missing or inconsistent."""


_HEADER = """\
# PDBClean stage configuration -- GENERATED, DO NOT EDIT
#
# This is the configuration the scientific workers of one frozen run execute.
# It is a deterministic projection of that run's canonical `resolved_run.json`
# (defaults -> profile -> explicit overrides), written once when the run was
# frozen and never regenerated afterwards.
#
# Editing this file does not change the run: `pdbclean.stage_runner` recomputes
# this projection from `resolved_run.json` and aborts the stage before any
# scientific work if the bytes differ.
#
# Canonical resolved configuration SHA256: {resolved_sha}
# Scientific configuration SHA256:         {scientific_sha}
# Schema:                                  {schema} v{version}
"""


def stage_config_document(resolved: ResolvedRunConfig) -> dict[str, Any]:
    """Return the projected stage configuration as a mapping.

    Pure: no clock, no environment, no filesystem.  ``${TMPDIR}`` and friends
    stay as templates and are expanded by ``load_config`` on the executing
    host.
    """

    try:
        document = resolved.to_protocol_config()
    except RunConfigError as exc:
        raise StageConfigError(str(exc)) from exc

    document["schema_name"] = STAGE_CONFIG_SCHEMA_NAME
    document["schema_version"] = STAGE_CONFIG_SCHEMA_VERSION

    return document


def stage_config_text(resolved: ResolvedRunConfig) -> str:
    """Return the exact bytes (as text) written to ``stage_config.yaml``."""

    header = _HEADER.format(
        resolved_sha=resolved.sha256,
        scientific_sha=resolved.scientific_sha256,
        schema=STAGE_CONFIG_SCHEMA_NAME,
        version=STAGE_CONFIG_SCHEMA_VERSION,
    )

    body = yaml.safe_dump(
        stage_config_document(resolved),
        sort_keys=True,
        default_flow_style=False,
        allow_unicode=True,
    )

    return header + "\n" + body


def stage_config_sha256(resolved: ResolvedRunConfig) -> str:
    """SHA256 of the projected document's file bytes."""

    return hashlib.sha256(
        stage_config_text(resolved).encode("utf-8")
    ).hexdigest()


def stage_config_path(directory: str | Path) -> Path:
    return Path(directory) / STAGE_CONFIG_BASENAME


def write_stage_config(
    resolved: ResolvedRunConfig,
    directory: str | Path,
) -> dict[str, str]:
    """Write ``stage_config.yaml`` into ``directory``.

    Writing is idempotent for the same configuration and refuses to change an
    existing document that differs -- a frozen run's executable configuration
    is written once, at freeze time, and is immutable thereafter.
    """

    target = Path(directory)
    target.mkdir(parents=True, exist_ok=True)

    path = stage_config_path(target)
    text = stage_config_text(resolved)

    if path.is_file():
        existing = path.read_text(encoding="utf-8")

        if existing != text:
            raise StageConfigError(
                f"Refusing to overwrite an existing frozen stage "
                f"configuration with different content: {path}"
            )
    else:
        path.write_text(text, encoding="utf-8")

    return {
        "stage_config_path": str(path),
        "stage_config_sha256": hashlib.sha256(
            text.encode("utf-8")
        ).hexdigest(),
    }


def cached_stage_config(
    resolved: ResolvedRunConfig,
    cache_root: str | Path,
) -> Path:
    """Materialise a content-addressed projection outside any run.

    ``pdbclean stage-command`` and ``pdbclean plan`` show what *would* run
    before a run exists.  The command they print has to be runnable and has to
    carry the operator's own resolved values, so it cannot point at the frozen
    base YAML.  It points here instead: a file named by the SHA256 of the
    canonical configuration it was projected from, which makes it immutable,
    self-verifying and shared between every preview of the same configuration.

    This is a *preview* artefact.  Execution goes through a frozen run, where
    the same bytes live in the run directory next to the provenance that
    explains them.
    """

    root = Path(cache_root)
    root.mkdir(parents=True, exist_ok=True)

    path = root / f"{resolved.sha256}.yaml"
    text = stage_config_text(resolved)

    if not path.is_file() or path.read_text(encoding="utf-8") != text:
        temporary = path.with_name(f".{path.name}.tmp")
        temporary.write_text(text, encoding="utf-8")
        temporary.replace(path)

    return path


def verify_stage_config(
    directory: str | Path,
    resolved: ResolvedRunConfig,
) -> dict[str, Any]:
    """Confirm the on-disk projection is the one this configuration produces.

    Returns the verification detail.  Raises :class:`StageConfigError` when the
    file is missing or has been modified, so a stage aborts *before* computing
    anything rather than publishing science from an edited configuration.
    """

    path = stage_config_path(directory)

    if not path.is_file():
        raise StageConfigError(
            f"Run has no executable stage configuration: {path}. Runs frozen "
            "before this file existed are legacy records; re-freeze the "
            "configuration to execute it."
        )

    observed = path.read_text(encoding="utf-8")
    expected = stage_config_text(resolved)

    observed_sha = hashlib.sha256(observed.encode("utf-8")).hexdigest()
    expected_sha = hashlib.sha256(expected.encode("utf-8")).hexdigest()

    if observed_sha != expected_sha:
        raise StageConfigError(
            "The stage configuration on disk is not the projection of this "
            f"run's canonical configuration.\n  file:     {path}\n"
            f"  observed: {observed_sha}\n  expected: {expected_sha}\n"
            "A frozen run is immutable; refusing to execute."
        )

    return {
        "stage_config_path": str(path),
        "stage_config_sha256": observed_sha,
        "verified": True,
    }


# ---------------------------------------------------------------------------
# The scientific values a worker actually loaded
# ---------------------------------------------------------------------------
#
# These are read back through the *production* loader and the *production*
# helpers, including the very fallbacks that used to hide a lost override.  If
# a projection ever stopped carrying `brain_filter`, `executed_values` would
# report the module constant 0.010 A and the comparison against the frozen
# configuration would fail the stage.


def executed_values(config: Mapping[str, Any]) -> dict[str, Any]:
    """Extract the scientific values a stage will compute with.

    ``config`` must be the mapping returned by
    :func:`pdbclean.config.load_config`, i.e. what the worker itself holds.
    """

    from pdbclean.brain_prefilter_production import BRAIN_PREFILTER_TAU_ANGSTROM
    from pdbclean.defaults import (
        brain_filter_threshold_units,
        near_duplicate_threshold_units,
        representation_precision_angstrom,
    )

    snapshot = config.get("snapshot") or {}
    selection = config.get("selection") or {}
    quality_rules = config.get("quality_rules") or {}
    geometry = config.get("post_cleaning_geometric_validation") or {}
    release = config.get("release") or {}

    # Reproduce brain_prefilter_production.main() exactly, fallback included.
    brain_section = config.get(
        "brain_filter",
        {"threshold_angstrom": BRAIN_PREFILTER_TAU_ANGSTROM},
    )

    values: dict[str, Any] = {
        "snapshot": snapshot.get("snapshot_id"),
        "snapshot_mode": snapshot.get("mode"),
        "protocol_version": release.get("protocol_version"),
        "model_policy": (selection.get("models") or {}).get("policy"),
        "model_id": (selection.get("models") or {}).get("model_id"),
        "minimum_backbone_distance_angstrom": (
            (quality_rules.get("backbone_distance") or {}).get(
                "minimum_distance_angstrom"
            )
        ),
        "minimum_triangle_angle_degrees": geometry.get(
            "minimum_triangle_angle_degrees"
        ),
        "representation_precision_angstrom": _safe(
            representation_precision_angstrom, config
        ),
        "brain_filter_threshold_angstrom": brain_section.get(
            "threshold_angstrom"
        ),
        "brain_filter_threshold_units": _safe(
            brain_filter_threshold_units,
            {
                "brain_filter": brain_section,
                "bri": config.get("bri", {}),
            },
        ),
        "complete_bri_near_duplicate_threshold_angstrom": (
            (config.get("duplicate_search") or {}).get(
                "near_duplicate_threshold_angstrom"
            )
        ),
        "complete_bri_near_duplicate_threshold_units": _safe(
            near_duplicate_threshold_units, config
        ),
    }

    return values


def _safe(function, config) -> Any:
    """Return the helper's value, or the error text it raised.

    A worker must be able to *record* that a value could not be derived; that
    is itself evidence, and it is what a mismatch check needs in order to fail
    rather than crash while writing provenance.
    """

    try:
        return function(config)
    except Exception as exc:  # noqa: BLE001 - recorded, then compared
        return f"ERROR: {type(exc).__name__}: {exc}"


def frozen_values(resolved: ResolvedRunConfig) -> dict[str, Any]:
    """The same scientific values, read from the frozen canonical config."""

    from pdbclean.defaults import (
        brain_filter_threshold_units,
        near_duplicate_threshold_units,
        representation_precision_angstrom,
    )

    return {
        "snapshot": resolved.get("snapshot.snapshot_id"),
        "snapshot_mode": "fixed",
        "protocol_version": resolved.get("release.protocol_version"),
        "model_policy": resolved.get("selection.models.policy"),
        "model_id": resolved.get("selection.models.model_id"),
        "minimum_backbone_distance_angstrom": resolved.get(
            "quality_rules.backbone_distance.minimum_distance_angstrom"
        ),
        "minimum_triangle_angle_degrees": resolved.get(
            "post_cleaning_geometric_validation."
            "minimum_triangle_angle_degrees"
        ),
        "representation_precision_angstrom": _safe(
            representation_precision_angstrom, resolved.data
        ),
        "brain_filter_threshold_angstrom": resolved.get(
            "brain_filter.threshold_angstrom"
        ),
        "brain_filter_threshold_units": _safe(
            brain_filter_threshold_units, resolved.data
        ),
        "complete_bri_near_duplicate_threshold_angstrom": resolved.get(
            "duplicate_search.near_duplicate_threshold_angstrom"
        ),
        "complete_bri_near_duplicate_threshold_units": _safe(
            near_duplicate_threshold_units, resolved.data
        ),
    }


def compare_values(
    executed: Mapping[str, Any],
    frozen: Mapping[str, Any],
) -> list[dict[str, Any]]:
    """Return every scientific value the worker did not load as frozen."""

    mismatches: list[dict[str, Any]] = []

    for key in sorted(set(frozen) | set(executed)):
        want = frozen.get(key)
        got = executed.get(key)

        if _agree(got, want):
            continue

        mismatches.append({"key": key, "executed": got, "frozen": want})

    return mismatches


def _agree(observed: Any, expected: Any) -> bool:
    if isinstance(observed, bool) or isinstance(expected, bool):
        return observed == expected

    if isinstance(observed, (int, float)) and isinstance(
        expected, (int, float)
    ):
        return abs(float(observed) - float(expected)) <= 1.0e-12

    return observed == expected
