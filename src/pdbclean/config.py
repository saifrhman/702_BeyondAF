"""Loading and validation for versioned PDBClean configuration files."""

from __future__ import annotations

import hashlib
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml


class ConfigError(ValueError):
    """Raised when a PDBClean configuration is invalid."""


class StrictSafeLoader(yaml.SafeLoader):
    """YAML loader that rejects duplicate mapping keys."""


def _construct_unique_mapping(
    loader: StrictSafeLoader,
    node: yaml.MappingNode,
    deep: bool = False,
) -> dict[str, Any]:
    mapping: dict[str, Any] = {}

    for key_node, value_node in node.value:
        key = loader.construct_object(key_node, deep=deep)

        if key in mapping:
            raise ConfigError(f"Duplicate YAML key: {key!r}")

        mapping[key] = loader.construct_object(value_node, deep=deep)

    return mapping


StrictSafeLoader.add_constructor(
    yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG,
    _construct_unique_mapping,
)


@dataclass(frozen=True)
class LoadedConfig:
    """Validated configuration plus provenance information."""

    path: Path
    data: dict[str, Any]
    sha256: str


REQUIRED_TOP_LEVEL_SECTIONS = (
    "release",
    "snapshot",
    "selection",
    "quality_rules",
    "bri",
    "geometric_search",
    "graph",
    "execution",
    "storage",
    "observability",
    "automation",
)


def _expand_environment(value: Any) -> Any:
    if isinstance(value, str):
        return os.path.expandvars(value)

    if isinstance(value, list):
        return [_expand_environment(item) for item in value]

    if isinstance(value, dict):
        return {
            key: _expand_environment(item)
            for key, item in value.items()
        }

    return value


def _validate_required_sections(config: dict[str, Any]) -> None:
    missing = [
        section
        for section in REQUIRED_TOP_LEVEL_SECTIONS
        if section not in config
    ]

    if missing:
        raise ConfigError(
            "Missing required top-level sections: "
            + ", ".join(missing)
        )


def _validate_release(config: dict[str, Any]) -> None:
    release = config["release"]

    if not isinstance(release, dict):
        raise ConfigError("release must be a mapping")

    dataset_name = release.get("dataset_name")

    if not isinstance(dataset_name, str) or not dataset_name.strip():
        raise ConfigError(
            "release.dataset_name must be a non-empty string"
        )

    protocol_version = release.get("protocol_version")

    if (
        not isinstance(protocol_version, str)
        or not protocol_version.strip()
    ):
        raise ConfigError(
            "release.protocol_version must be a non-empty string"
        )


def _validate_snapshot(config: dict[str, Any]) -> None:
    snapshot = config["snapshot"]

    if not isinstance(snapshot, dict):
        raise ConfigError("snapshot must be a mapping")

    mode = snapshot.get("mode")

    if mode not in {"fixed", "latest_complete"}:
        raise ConfigError(
            "snapshot.mode must be either 'fixed' or 'latest_complete'"
        )

    bucket_url = snapshot.get("bucket_url")

    if not isinstance(bucket_url, str) or not bucket_url.strip():
        raise ConfigError(
            "snapshot.bucket_url must be a non-empty string"
        )

    if mode == "fixed":
        snapshot_id = snapshot.get("snapshot_id")

        if (
            not isinstance(snapshot_id, str)
            or len(snapshot_id) != 8
            or not snapshot_id.isdigit()
        ):
            raise ConfigError(
                "fixed snapshot mode requires snapshot.snapshot_id "
                "in YYYYMMDD format"
            )

        expected_count = snapshot.get("expected_mmcif_count")
        expected_bytes = snapshot.get("expected_total_bytes")

        if expected_count is not None:
            if not isinstance(expected_count, int) or expected_count <= 0:
                raise ConfigError(
                    "snapshot.expected_mmcif_count must be a "
                    "positive integer"
                )

        if expected_bytes is not None:
            if not isinstance(expected_bytes, int) or expected_bytes <= 0:
                raise ConfigError(
                    "snapshot.expected_total_bytes must be a "
                    "positive integer"
                )

    if mode == "latest_complete":
        if "snapshot_id" in snapshot:
            raise ConfigError(
                "latest_complete mode must not define snapshot_id"
            )


def _validate_selection(config: dict[str, Any]) -> None:
    selection = config["selection"]

    if not isinstance(selection, dict):
        raise ConfigError("selection must be a mapping")

    models = selection.get("models")

    if not isinstance(models, dict):
        raise ConfigError("selection.models must be a mapping")

    policy = models.get("policy")

    if policy != "first_model":
        raise ConfigError(
            "selection.models.policy must be 'first_model'"
        )

    model_id = models.get("model_id")

    if (
        not isinstance(model_id, int)
        or isinstance(model_id, bool)
        or model_id <= 0
    ):
        raise ConfigError(
            "selection.models.model_id must be a positive integer"
        )


def _validate_execution(config: dict[str, Any]) -> None:
    execution = config["execution"]

    if not isinstance(execution, dict):
        raise ConfigError("execution must be a mapping")

    positive_integer_fields = (
        "batch_size",
        "download_concurrency",
        "connection_timeout_seconds",
    )

    for field in positive_integer_fields:
        value = execution.get(field)

        if (
            not isinstance(value, int)
            or isinstance(value, bool)
            or value <= 0
        ):
            raise ConfigError(
                f"execution.{field} must be a positive integer"
            )

    max_retries = execution.get("max_retries")

    if (
        not isinstance(max_retries, int)
        or isinstance(max_retries, bool)
        or max_retries < 0
    ):
        raise ConfigError(
            "execution.max_retries must be a non-negative integer"
        )

    for field in (
        "atomic_writes",
        "write_success_markers",
    ):
        if not isinstance(execution.get(field), bool):
            raise ConfigError(
                f"execution.{field} must be a boolean"
            )


def load_config(path: str | Path) -> LoadedConfig:
    """Load, expand, validate, and checksum a YAML configuration."""

    config_path = Path(path).resolve()

    if not config_path.is_file():
        raise ConfigError(
            f"Configuration file does not exist: {config_path}"
        )

    raw_bytes = config_path.read_bytes()
    sha256 = hashlib.sha256(raw_bytes).hexdigest()

    try:
        loaded = yaml.load(
            raw_bytes.decode("utf-8"),
            Loader=StrictSafeLoader,
        )
    except yaml.YAMLError as exc:
        raise ConfigError(f"Invalid YAML: {exc}") from exc

    if not isinstance(loaded, dict):
        raise ConfigError("Configuration root must be a mapping")

    expanded = _expand_environment(loaded)

    _validate_required_sections(expanded)
    _validate_release(expanded)
    _validate_snapshot(expanded)
    _validate_selection(expanded)
    _validate_execution(expanded)

    return LoadedConfig(
        path=config_path,
        data=expanded,
        sha256=sha256,
    )
