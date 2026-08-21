"""Deterministic run-configuration resolution for the COMP702 pipeline.

Precedence, from weakest to strongest:

1. ``pdbclean.defaults.VALIDATED_DEFAULTS`` -- the built-in validated defaults
   that reproduce the frozen COMP702 scientific pipeline;
2. a configuration file (``--config`` / UI profile selection);
3. explicit key overrides (``--set a.b=c`` from the CLI, or UI form fields).

A later layer replaces a value supplied by an earlier layer.  Mappings merge
recursively; scalars and sequences are replaced wholesale.  The resolver records
which layer supplied every leaf, so provenance can answer "where did this value
come from" for any scientific parameter.

The resolved configuration -- not the input file -- is what a run executes and
what is written to provenance.
"""

from __future__ import annotations

import copy
import hashlib
import json
import os
import re
import socket
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Mapping

import yaml

from pdbclean.config import StrictSafeLoader
from pdbclean.defaults import (
    DEFAULTS_VERSION,
    brain_filter_threshold_milliangstrom,
    brain_filter_threshold_units,
    near_duplicate_threshold_milliangstrom,
    near_duplicate_threshold_units,
    precision_is_implemented,
    representation_precision_angstrom,
    representation_unit_label,
    validated_defaults,
)


RESOLVED_CONFIG_SCHEMA_NAME = "pdbclean_resolved_run_config"
RESOLVED_CONFIG_SCHEMA_VERSION = "1.0"

LAYER_BUILTIN = "builtin_default"
LAYER_CONFIG_FILE = "config_file"
LAYER_OVERRIDE = "override"

#: Scientific parameters that must always be visible in a run summary.  These
#: are display keys only; they never affect resolution.
SCIENTIFIC_SUMMARY_KEYS: tuple[tuple[str, str], ...] = (
    ("snapshot.mode", "Snapshot selection mode"),
    ("snapshot.snapshot_id", "Snapshot identity"),
    ("selection.models.policy", "Model policy"),
    ("selection.models.model_id", "Model scope"),
    ("selection.chain_namespace.canonical", "Canonical chain namespace"),
    ("quality_rules.backbone_distance.minimum_distance_angstrom",
     "Q005 minimum backbone distance (A)"),
    ("post_cleaning_geometric_validation.minimum_triangle_angle_degrees",
     "Minimum N-CA-C angle (deg)"),
    ("bri.implementation_version", "BRI implementation version"),
    ("bri.representation_precision_angstrom", "BRI representation precision (A)"),
    ("brain.dimension", "Brain dimension"),
    ("brain_filter.threshold_angstrom", "Brain filter threshold (A)"),
    ("brain_filter.operator", "Brain filter operator"),
    ("duplicate_search.near_duplicate_threshold_angstrom",
     "Complete-BRI near-duplicate threshold (A)"),
    ("duplicate_search.operator", "Near-duplicate operator"),
    ("duplicate_search.metric", "Distance metric"),
    ("duplicate_search.search_engine", "Production NN engine"),
    ("representative_selection.policy_name", "Representative policy"),
    ("representative_selection.m1_policy", "m=1 policy"),
    ("graph.require_direct_edge_for_removal", "Direct-edge removal required"),
    ("graph.connected_component_is_duplicate_equivalence",
     "Components treated as equivalence"),
)


class RunConfigError(ValueError):
    """Raised when a run configuration cannot be resolved."""


def _is_mapping(value: Any) -> bool:
    return isinstance(value, Mapping)


def _deep_merge(
    base: dict[str, Any],
    incoming: Mapping[str, Any],
    *,
    layer: str,
    sources: dict[str, str],
    prefix: str = "",
) -> dict[str, Any]:
    """Recursively merge ``incoming`` into ``base``, recording leaf sources."""

    result = base

    for key, value in incoming.items():
        path = f"{prefix}{key}"

        if _is_mapping(value) and _is_mapping(result.get(key)):
            _deep_merge(
                result[key],
                value,
                layer=layer,
                sources=sources,
                prefix=f"{path}.",
            )
            continue

        if _is_mapping(value):
            result[key] = {}
            _deep_merge(
                result[key],
                value,
                layer=layer,
                sources=sources,
                prefix=f"{path}.",
            )
            continue

        result[key] = copy.deepcopy(value)
        sources[path] = layer

    return result


def _record_leaf_sources(
    value: Any,
    *,
    layer: str,
    sources: dict[str, str],
    prefix: str = "",
) -> None:
    if _is_mapping(value):
        for key, child in value.items():
            _record_leaf_sources(
                child,
                layer=layer,
                sources=sources,
                prefix=f"{prefix}{key}.",
            )
        return

    sources[prefix.rstrip(".")] = layer


def parse_override(item: str) -> tuple[str, Any]:
    """Parse one ``dotted.key=value`` override.

    The value is parsed with YAML scalar rules, so ``0.010`` becomes a float,
    ``10`` an int, ``true`` a bool, ``null`` None, and anything else a string.
    """

    if not isinstance(item, str) or "=" not in item:
        raise RunConfigError(
            f"Override must look like dotted.key=value, got {item!r}"
        )

    key, _, raw = item.partition("=")
    key = key.strip()

    if not key:
        raise RunConfigError(f"Override has an empty key: {item!r}")

    for part in key.split("."):
        if not part:
            raise RunConfigError(
                f"Override key has an empty path segment: {item!r}"
            )

    try:
        value = yaml.safe_load(raw)
    except yaml.YAMLError as exc:
        raise RunConfigError(
            f"Override value is not valid YAML: {item!r}"
        ) from exc

    return key, value


def overrides_to_mapping(
    overrides: Iterable[str] | Mapping[str, Any] | None,
) -> dict[str, Any]:
    """Turn CLI/UI overrides into a nested mapping."""

    if overrides is None:
        return {}

    pairs: list[tuple[str, Any]]

    if isinstance(overrides, Mapping):
        pairs = [(str(k), v) for k, v in overrides.items()]
    else:
        pairs = [parse_override(item) for item in overrides]

    nested: dict[str, Any] = {}

    for key, value in pairs:
        cursor = nested

        parts = key.split(".")

        for part in parts[:-1]:
            existing = cursor.get(part)

            if not isinstance(existing, dict):
                existing = {}
                cursor[part] = existing

            cursor = existing

        cursor[parts[-1]] = value

    return nested


def canonical_config_bytes(config: Mapping[str, Any]) -> bytes:
    """Return the canonical byte serialisation used for config hashing."""

    return json.dumps(
        config,
        sort_keys=True,
        ensure_ascii=False,
        separators=(",", ":"),
        allow_nan=False,
    ).encode("utf-8")


def config_sha256(config: Mapping[str, Any]) -> str:
    """Return the SHA256 of the canonical configuration serialisation."""

    return hashlib.sha256(canonical_config_bytes(config)).hexdigest()


#: Sections that define the scientific method and the dataset it is applied to.
#: Infrastructure sections (storage, execution, observability) are excluded so
#: that the scientific identity of a run does not change merely because it ran
#: on a node with a different ``$TMPDIR``.
SCIENTIFIC_SECTIONS: tuple[str, ...] = (
    "release",
    "selection",
    "quality_rules",
    "post_cleaning_geometric_validation",
    "bri",
    "brain",
    "brain_filter",
    "duplicate_search",
    "graph",
    "representative_selection",
)

#: Snapshot keys that carry scientific identity.  The selection *mode* and the
#: archive layout details do not: two runs that resolve to the same snapshot
#: identity are scientifically the same run whether the snapshot was named
#: explicitly or discovered as the latest complete one.
SCIENTIFIC_SNAPSHOT_KEYS: tuple[str, ...] = ("snapshot_id", "bucket_url")


def scientific_projection(config: Mapping[str, Any]) -> dict[str, Any]:
    """Return only the parts of a configuration that define the science."""

    projection: dict[str, Any] = {
        section: copy.deepcopy(config[section])
        for section in SCIENTIFIC_SECTIONS
        if section in config
    }

    snapshot = config.get("snapshot")

    if isinstance(snapshot, Mapping):
        projection["snapshot"] = {
            key: snapshot.get(key)
            for key in SCIENTIFIC_SNAPSHOT_KEYS
            if key in snapshot
        }

    projection["defaults_version"] = config.get("defaults_version")

    return projection


def scientific_sha256(config: Mapping[str, Any]) -> str:
    """Return the SHA256 of the scientific projection of a configuration."""

    return config_sha256(scientific_projection(config))


#: Environment variables whose value belongs to the *host a job happens to run
#: on*, not to the run.  They are deliberately left unexpanded in the canonical
#: resolved configuration.
#:
#: ``${TMPDIR}/pdbclean`` is a policy -- "use this node's own scratch" -- and
#: that policy is the same statement on every node.  Baking a login node's
#: ``$TMPDIR`` into the configuration that a compute node will execute would be
#: both wrong and a source of spurious hash drift, so the template is preserved
#: and resolved at execution time instead.  The concrete value the job actually
#: used is recorded as runtime provenance (see :func:`runtime_environment`).
RUNTIME_ENVIRONMENT_VARIABLES: frozenset[str] = frozenset(
    {
        "TMPDIR",
        "TMP",
        "TEMP",
        "HOSTNAME",
        "SLURMD_NODENAME",
        "SLURM_JOB_ID",
        "SLURM_JOBID",
        "SLURM_ARRAY_JOB_ID",
        "SLURM_ARRAY_TASK_ID",
        "SLURM_JOB_NODELIST",
        "SLURM_JOB_PARTITION",
        "SLURM_SUBMIT_HOST",
        "SLURM_SUBMIT_DIR",
    }
)

_ENV_REFERENCE = re.compile(r"\$\{([A-Za-z_][A-Za-z0-9_]*)\}|\$([A-Za-z_][A-Za-z0-9_]*)")


def _expand_environment(value: Any, *, preserve_runtime: bool = True) -> Any:
    """Expand ``${VAR}`` references in string leaves.

    Variables in :data:`RUNTIME_ENVIRONMENT_VARIABLES` are preserved verbatim
    when ``preserve_runtime`` is true, so the canonical configuration stays
    identical on every host.  Pass ``preserve_runtime=False`` at execution time
    to obtain the concrete host-specific paths.
    """

    if isinstance(value, str):
        def _substitute(match: "re.Match[str]") -> str:
            name = match.group(1) or match.group(2)

            if preserve_runtime and name in RUNTIME_ENVIRONMENT_VARIABLES:
                return match.group(0)

            return os.environ.get(name, match.group(0))

        return _ENV_REFERENCE.sub(_substitute, value)

    if isinstance(value, dict):
        return {
            k: _expand_environment(v, preserve_runtime=preserve_runtime)
            for k, v in value.items()
        }

    if isinstance(value, list):
        return [
            _expand_environment(v, preserve_runtime=preserve_runtime)
            for v in value
        ]

    return value


def runtime_environment(resolved: "ResolvedRunConfig") -> dict[str, Any]:
    """Resolve the host-specific values a job actually runs with.

    The canonical configuration states the *policy* (``${TMPDIR}/pdbclean``);
    this returns what that policy means on this host, right now.  The result is
    recorded as execution provenance and never folded back into the canonical
    configuration or into either configuration hash.
    """

    storage = resolved.get("storage") or {}

    payload: dict[str, Any] = {
        "hostname": socket.gethostname(),
        "environment": {
            name: os.environ[name]
            for name in sorted(RUNTIME_ENVIRONMENT_VARIABLES)
            if name in os.environ
        },
        "storage": {
            key: _expand_environment(value, preserve_runtime=False)
            for key, value in storage.items()
            if isinstance(value, str)
        },
    }

    return payload


def load_config_file(path: str | Path) -> dict[str, Any]:
    """Load a YAML configuration layer without legacy schema validation.

    ``pdbclean.config.load_config`` validates the *complete* legacy protocol
    schema.  A run-configuration layer is allowed to be partial, so this loader
    performs syntax and type checks only.  The merged result is validated by
    :func:`validate_resolved_config`.
    """

    config_path = Path(path)

    if not config_path.is_file():
        raise RunConfigError(f"Configuration file not found: {config_path}")

    try:
        loaded = yaml.load(
            config_path.read_text(encoding="utf-8"),
            Loader=StrictSafeLoader,
        )
    except yaml.YAMLError as exc:
        raise RunConfigError(
            f"Invalid YAML in {config_path}: {exc}"
        ) from exc

    if loaded is None:
        return {}

    if not isinstance(loaded, dict):
        raise RunConfigError(
            f"Configuration root must be a mapping: {config_path}"
        )

    return loaded


@dataclass(frozen=True)
class ResolvedRunConfig:
    """A fully resolved, hashable run configuration."""

    data: dict[str, Any]
    sources: dict[str, str] = field(default_factory=dict)
    layers: tuple[str, ...] = ()
    config_path: str | None = None
    override_items: tuple[str, ...] = ()

    @property
    def sha256(self) -> str:
        """SHA256 of the complete resolved configuration."""

        return config_sha256(self.data)

    @property
    def scientific_sha256(self) -> str:
        """SHA256 of the scientific projection only.

        Two runs with the same value here are scientifically the same run.
        This is the identity used to decide whether existing stage output may
        be reused, and it is stable across machines and ``$TMPDIR`` values.
        """

        return scientific_sha256(self.data)

    def get(self, dotted: str, default: Any = None) -> Any:
        cursor: Any = self.data

        for part in dotted.split("."):
            if not isinstance(cursor, Mapping) or part not in cursor:
                return default

            cursor = cursor[part]

        return cursor

    def source_of(self, dotted: str) -> str:
        return self.sources.get(dotted, LAYER_BUILTIN)

    @property
    def near_duplicate_threshold_mA(self) -> int:
        return near_duplicate_threshold_units(self.data)

    @property
    def brain_threshold_mA(self) -> int:
        return brain_filter_threshold_units(self.data)

    def scientific_summary(self) -> list[dict[str, Any]]:
        """Return the display rows shown before a run starts."""

        rows: list[dict[str, Any]] = []

        for dotted, label in SCIENTIFIC_SUMMARY_KEYS:
            rows.append(
                {
                    "key": dotted,
                    "label": label,
                    "value": self.get(dotted),
                    "source": self.source_of(dotted),
                }
            )

        return rows

    def to_dict(self) -> dict[str, Any]:
        return copy.deepcopy(self.data)

    def to_yaml(self) -> str:
        return yaml.safe_dump(
            self.data,
            sort_keys=True,
            default_flow_style=False,
            allow_unicode=True,
        )

    def to_protocol_config(self) -> dict[str, Any]:
        """Project the resolved run config onto the legacy protocol schema.

        The existing stage entry points consume the Protocol 3.2 configuration
        schema validated by :func:`pdbclean.config.load_config`.  This
        projection carries the resolved scientific values into that schema so
        that CLI, UI and Slurm execution all drive the same stage code with the
        same numbers.

        The snapshot is always emitted as ``fixed`` with its concrete resolved
        identity, because a run must be pinned to an immutable snapshot before
        any processing begins.

        This projection is the *execution* handoff, so runtime templates such
        as ``${TMPDIR}`` are resolved to their concrete value on the executing
        host here -- and only here.  The canonical configuration itself keeps
        the template, so it stays identical on every host.
        """

        data = self.to_dict()

        snapshot_id = self.get("snapshot.snapshot_id")

        if not snapshot_id:
            raise RunConfigError(
                "Cannot project a protocol configuration before the "
                "snapshot has been resolved to a concrete identity"
            )

        snapshot_section: dict[str, Any] = {
            "mode": "fixed",
            "snapshot_id": snapshot_id,
            "bucket_url": self.get("snapshot.bucket_url"),
        }

        expected_count = self.get("snapshot.expected_mmcif_count")

        if expected_count is not None:
            snapshot_section["expected_mmcif_count"] = expected_count

        expected_bytes = self.get("snapshot.expected_total_bytes")

        if expected_bytes is not None:
            snapshot_section["expected_total_bytes"] = expected_bytes

        projected: dict[str, Any] = {
            "release": data["release"],
            "snapshot": snapshot_section,
            "selection": data["selection"],
            "quality_rules": data["quality_rules"],
            "post_cleaning_geometric_validation": data[
                "post_cleaning_geometric_validation"
            ],
            "bri": {
                "enabled": self.get("bri.enabled"),
                "implementation": self.get("bri.implementation"),
                "vector_partition_key": self.get("bri.vector_partition_key"),
                "require_full_accounting": self.get(
                    "bri.require_full_accounting"
                ),
            },
            "duplicate_search": {
                "near_duplicate_threshold_angstrom": self.get(
                    "duplicate_search.near_duplicate_threshold_angstrom"
                ),
            },
            "graph": {
                "build_connected_components": self.get(
                    "graph.build_connected_components"
                ),
                "calculate_edge_density": self.get(
                    "graph.calculate_edge_density"
                ),
                "calculate_clique_status": self.get(
                    "graph.calculate_clique_status"
                ),
                "require_direct_edge_for_removal": self.get(
                    "graph.require_direct_edge_for_removal"
                ),
            },
            "execution": {
                key: value
                for key, value in data["execution"].items()
                if key != "executor"
            },
            "storage": {
                "temporary_root": _expand_environment(
                    self.get("storage.temporary_root"),
                    preserve_runtime=False,
                ),
                "output_root": _expand_environment(
                    self.get("storage.output_root"),
                    preserve_runtime=False,
                ),
                "retain_downloaded_mmcif": self.get(
                    "storage.retain_downloaded_mmcif"
                ),
            },
            "observability": data["observability"],
            "automation": {
                "snapshot_watcher_enabled": False,
                "require_stable_manifest": True,
                "prevent_concurrent_release_runs": True,
            },
        }

        return projected


def normalise_resolved_config(config: dict[str, Any]) -> dict[str, Any]:
    """Repair representations that YAML scalar rules get wrong.

    A snapshot identity is a date-shaped *string*.  Both ``--set
    snapshot.snapshot_id=20260101`` and an unquoted ``snapshot_id: 20260101``
    in a YAML file arrive here as an integer, which is the same snapshot --
    not a different one.  Coercing it back is a representation repair; it
    cannot select a different snapshot than the user asked for.
    """

    snapshot = config.get("snapshot")

    if _is_mapping(snapshot):
        snapshot_id = snapshot.get("snapshot_id")

        if isinstance(snapshot_id, int) and not isinstance(snapshot_id, bool):
            snapshot["snapshot_id"] = f"{snapshot_id:08d}"

    return config


def validate_resolved_config(config: Mapping[str, Any]) -> None:
    """Validate the merged run configuration.

    Only checks that would change a scientific result are enforced here.  The
    thresholds are re-validated through the same helpers the production stages
    use, so an override that cannot be represented exactly is rejected before a
    run identity is created rather than deep inside a Slurm array.
    """

    for section in (
        "release",
        "snapshot",
        "selection",
        "quality_rules",
        "post_cleaning_geometric_validation",
        "bri",
        "brain",
        "brain_filter",
        "duplicate_search",
        "graph",
        "representative_selection",
        "execution",
        "storage",
    ):
        if section not in config:
            raise RunConfigError(
                f"Resolved configuration is missing section: {section}"
            )

        if not _is_mapping(config[section]):
            raise RunConfigError(
                f"Resolved configuration section must be a mapping: {section}"
            )

    mode = config["snapshot"].get("mode")

    if mode not in {"fixed", "latest_complete", "explicit"}:
        raise RunConfigError(
            "snapshot.mode must be 'latest_complete', 'fixed' or 'explicit'"
        )

    bucket_url = config["snapshot"].get("bucket_url")

    if not isinstance(bucket_url, str) or not bucket_url.strip():
        raise RunConfigError("snapshot.bucket_url must be a non-empty string")

    snapshot_id = config["snapshot"].get("snapshot_id")

    if mode in {"fixed", "explicit"}:
        if (
            not isinstance(snapshot_id, str)
            or len(snapshot_id) != 8
            or not snapshot_id.isdigit()
        ):
            raise RunConfigError(
                "An explicit snapshot requires snapshot.snapshot_id in "
                "YYYYMMDD format"
            )
    elif snapshot_id is not None:
        if (
            not isinstance(snapshot_id, str)
            or len(snapshot_id) != 8
            or not snapshot_id.isdigit()
        ):
            raise RunConfigError(
                "snapshot.snapshot_id must be null or YYYYMMDD"
            )

    model_id = config["selection"].get("models", {}).get("model_id")

    if isinstance(model_id, bool) or not isinstance(model_id, int):
        raise RunConfigError("selection.models.model_id must be an integer")

    # Re-use the exact validators the production stages use, so an override
    # that cannot be represented exactly fails here rather than deep inside a
    # Slurm array.
    # Thresholds are validated against the CONFIGURED precision grid, not
    # against a hard-coded milliangstrom grid, so that a precision study is
    # rejected cleanly rather than silently rounding a threshold.
    try:
        precision = representation_precision_angstrom(config)
        near_duplicate_milliangstrom = near_duplicate_threshold_units(config)
        brain_milliangstrom = brain_filter_threshold_units(config)
    except ValueError as exc:
        raise RunConfigError(str(exc)) from exc

    for name, value in (
        ("duplicate_search.operator", config["duplicate_search"].get("operator")),
        ("brain_filter.operator", config["brain_filter"].get("operator")),
    ):
        if value != "less_than_or_equal":
            raise RunConfigError(
                f"{name} must be 'less_than_or_equal'; the validated "
                "COMP702 near-duplicate criterion is inclusive. Changing it "
                "is a scientific change and is not available as an override."
            )

    if config["brain"].get("dimension") != 9:
        raise RunConfigError(
            "brain.dimension is fixed at 9 by Definition 5.1 and is not "
            "configurable"
        )

    if config["duplicate_search"].get("metric") != "L_infinity":
        raise RunConfigError(
            "duplicate_search.metric is fixed at 'L_infinity' and is not "
            "configurable"
        )

    if config["duplicate_search"].get("final_classification_basis") != (
        "complete_BRI"
    ):
        raise RunConfigError(
            "duplicate_search.final_classification_basis is fixed at "
            "'complete_BRI'; Brain is a filtering layer only"
        )

    if brain_milliangstrom < near_duplicate_milliangstrom:
        raise RunConfigError(
            "brain_filter.threshold_angstrom "
            f"({brain_milliangstrom} mA) must not be smaller than "
            "duplicate_search.near_duplicate_threshold_angstrom "
            f"({near_duplicate_milliangstrom} mA); the Brain filter is a "
            "lossless prefilter and a smaller radius would discard true "
            "complete-BRI near duplicates"
        )


def resolve_run_config(
    *,
    config_path: str | Path | None = None,
    overrides: Iterable[str] | Mapping[str, Any] | None = None,
    override_origin: str = "cli",
    defaults: Mapping[str, Any] | None = None,
    expand_environment: bool = True,
    validate: bool = True,
) -> ResolvedRunConfig:
    """Resolve a run configuration from defaults, file and overrides."""

    base = (
        copy.deepcopy(dict(defaults))
        if defaults is not None
        else validated_defaults()
    )

    sources: dict[str, str] = {}
    _record_leaf_sources(base, layer=LAYER_BUILTIN, sources=sources)

    layers: list[str] = [LAYER_BUILTIN]

    resolved_path: str | None = None

    if config_path is not None:
        resolved_path = str(Path(config_path))
        file_layer = f"{LAYER_CONFIG_FILE}:{resolved_path}"
        file_data = load_config_file(config_path)

        # A profile may point at a base configuration file; that base is
        # applied first so the profile itself remains the stronger layer.
        base_config = file_data.pop("base_config", None)

        if base_config is not None:
            base_path = Path(base_config)

            if not base_path.is_absolute():
                base_path = Path(config_path).parent / base_path

            base_layer = f"{LAYER_CONFIG_FILE}:{base_path}"
            _deep_merge(
                base,
                load_config_file(base_path),
                layer=base_layer,
                sources=sources,
            )
            layers.append(base_layer)

        _deep_merge(base, file_data, layer=file_layer, sources=sources)
        layers.append(file_layer)

    override_mapping = overrides_to_mapping(overrides)

    override_items: tuple[str, ...] = ()

    if override_mapping:
        override_layer = f"{LAYER_OVERRIDE}:{override_origin}"
        _deep_merge(
            base,
            override_mapping,
            layer=override_layer,
            sources=sources,
        )
        layers.append(override_layer)

        if isinstance(overrides, Mapping):
            override_items = tuple(
                f"{key}={value}" for key, value in sorted(overrides.items())
            )
        elif overrides is not None:
            override_items = tuple(str(item) for item in overrides)

    if expand_environment:
        base = _expand_environment(base)

    base = normalise_resolved_config(base)

    base["schema_name"] = RESOLVED_CONFIG_SCHEMA_NAME
    base["schema_version"] = RESOLVED_CONFIG_SCHEMA_VERSION
    base["defaults_version"] = DEFAULTS_VERSION

    if validate:
        validate_resolved_config(base)

    return ResolvedRunConfig(
        data=base,
        sources=sources,
        layers=tuple(layers),
        config_path=resolved_path,
        override_items=override_items,
    )


def with_resolved_snapshot(
    resolved: ResolvedRunConfig,
    *,
    snapshot_id: str,
    layout: str | None = None,
    source_prefix: str | None = None,
    sample_mmcif_key: str | None = None,
    selection_mode: str | None = None,
) -> ResolvedRunConfig:
    """Return a copy of ``resolved`` pinned to a concrete snapshot identity."""

    data = resolved.to_dict()

    snapshot = data["snapshot"]
    snapshot["snapshot_id"] = snapshot_id
    snapshot["resolved_selection_mode"] = (
        selection_mode
        if selection_mode is not None
        else snapshot.get("mode")
    )

    if layout is not None:
        snapshot["layout"] = layout

    if source_prefix is not None:
        snapshot["source_prefix"] = source_prefix

    if sample_mmcif_key is not None:
        snapshot["sample_mmcif_key"] = sample_mmcif_key

    sources = dict(resolved.sources)
    sources["snapshot.snapshot_id"] = "snapshot_resolution"

    for key in (
        "resolved_selection_mode",
        "layout",
        "source_prefix",
        "sample_mmcif_key",
    ):
        if key in snapshot:
            sources[f"snapshot.{key}"] = "snapshot_resolution"

    validate_resolved_config(data)

    return ResolvedRunConfig(
        data=data,
        sources=sources,
        layers=resolved.layers + ("snapshot_resolution",),
        config_path=resolved.config_path,
        override_items=resolved.override_items,
    )


def write_resolved_config(
    resolved: ResolvedRunConfig,
    directory: str | Path,
    *,
    basename: str = "resolved_run",
) -> dict[str, str]:
    """Write ``resolved_run.yaml`` / ``.json`` and return their paths+hash."""

    target = Path(directory)
    target.mkdir(parents=True, exist_ok=True)

    yaml_path = target / f"{basename}.yaml"
    json_path = target / f"{basename}.json"

    payload = resolved.to_dict()

    yaml_path.write_text(resolved.to_yaml(), encoding="utf-8")
    json_path.write_bytes(canonical_config_bytes(payload))

    return {
        "resolved_config_yaml": str(yaml_path),
        "resolved_config_json": str(json_path),
        "resolved_config_sha256": resolved.sha256,
        "scientific_config_sha256": resolved.scientific_sha256,
    }


def load_resolved_config(
    directory: str | Path,
    *,
    basename: str = "resolved_run",
) -> ResolvedRunConfig:
    """Load the canonical configuration a run was created with.

    A run resolves its configuration **once**, at creation, and persists it.
    Downstream work -- a later stage, a Slurm job on another node, a UI
    inspecting the run -- must *load* that document rather than re-resolve it,
    so the run's identity cannot drift with the host it happens to execute on.

    The loaded configuration is returned verbatim: no defaults are re-layered
    and no environment is expanded, so its hashes are exactly the hashes that
    were recorded when the run was created.
    """

    json_path = Path(directory) / f"{basename}.json"

    if not json_path.is_file():
        raise RunConfigError(
            f"No canonical resolved configuration in {directory}"
        )

    data = json.loads(json_path.read_text(encoding="utf-8"))

    if not isinstance(data, dict):
        raise RunConfigError(
            f"Canonical resolved configuration is not a mapping: {json_path}"
        )

    return ResolvedRunConfig(
        data=data,
        sources={},
        layers=("loaded_from_run",),
        config_path=str(json_path),
    )
