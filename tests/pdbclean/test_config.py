"""Tests for strict PDBClean configuration loading."""

from pathlib import Path

import pytest
import yaml

from pdbclean.config import ConfigError, load_config


CONFIG_PATH = Path(
    "config/pdbclean/protocol_3_2_comp702_v1.yaml"
)


def test_load_valid_project_config(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setenv("TMPDIR", "/tmp/test-user")

    loaded = load_config(CONFIG_PATH)

    assert loaded.data["snapshot"]["mode"] == "fixed"
    assert loaded.data["snapshot"]["snapshot_id"] == "20260101"
    assert (
        loaded.data["storage"]["temporary_root"]
        == "/tmp/test-user/pdbclean"
    )
    assert len(loaded.sha256) == 64
    assert loaded.path == CONFIG_PATH.resolve()


def test_missing_required_section_is_rejected(
    tmp_path: Path,
) -> None:
    config_path = tmp_path / "missing.yaml"
    config_path.write_text(
        """
release:
  snapshot: "20260101"
""".strip()
    )

    with pytest.raises(
        ConfigError,
        match="Missing required top-level sections",
    ):
        load_config(config_path)


def test_duplicate_yaml_key_is_rejected(
    tmp_path: Path,
) -> None:
    config_path = tmp_path / "duplicate.yaml"
    config_path.write_text(
        """
release:
  snapshot: "20260101"
release:
  snapshot: "20260101"
""".strip()
    )

    with pytest.raises(
        ConfigError,
        match="Duplicate YAML key",
    ):
        load_config(config_path)


def test_invalid_snapshot_mode_is_rejected(
    tmp_path: Path,
) -> None:
    original = CONFIG_PATH.read_text()
    modified = original.replace(
        "mode: fixed",
        "mode: invalid",
        1,
    )

    config_path = tmp_path / "invalid_mode.yaml"
    config_path.write_text(modified)

    with pytest.raises(
        ConfigError,
        match="snapshot.mode must be either",
    ):
        load_config(config_path)


def test_fixed_mode_requires_snapshot_id(
    tmp_path: Path,
) -> None:
    data = yaml.safe_load(CONFIG_PATH.read_text())
    data["snapshot"].pop("snapshot_id")

    config_path = tmp_path / "missing_snapshot_id.yaml"
    config_path.write_text(
        yaml.safe_dump(data, sort_keys=False)
    )

    with pytest.raises(
        ConfigError,
        match="fixed snapshot mode requires",
    ):
        load_config(config_path)


def test_latest_complete_rejects_snapshot_id(
    tmp_path: Path,
) -> None:
    data = yaml.safe_load(CONFIG_PATH.read_text())
    data["snapshot"]["mode"] = "latest_complete"

    config_path = tmp_path / "latest_with_id.yaml"
    config_path.write_text(
        yaml.safe_dump(data, sort_keys=False)
    )

    with pytest.raises(
        ConfigError,
        match="latest_complete mode must not define snapshot_id",
    ):
        load_config(config_path)
