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


def _validate_snapshot(config: dict[str, Any]) -> None:
    snapshot = config["snapshot"]

    if not isinstance(snapshot, dict):
        raise ConfigError("snapshot must be a mapping")

    snapshot_id = snapshot.get("source_prefix", "").split("/", 1)[0]
    release_snapshot = config["release"].get("snapshot")

    if snapshot_id != release_snapshot:
        raise ConfigError(
            "release.snapshot does not match snapshot.source_prefix"
        )

    expected_count = snapshot.get("expected_mmcif_count")
    expected_bytes = snapshot.get("expected_total_bytes")

    if not isinstance(expected_count, int) or expected_count <= 0:
        raise ConfigError(
            "snapshot.expected_mmcif_count must be a positive integer"
        )

    if not isinstance(expected_bytes, int) or expected_bytes <= 0:
        raise ConfigError(
            "snapshot.expected_total_bytes must be a positive integer"
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
    _validate_snapshot(expanded)

    return LoadedConfig(
        path=config_path,
        data=expanded,
        sha256=sha256,
    )
