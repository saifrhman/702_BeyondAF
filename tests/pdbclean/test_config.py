"""Tests for strict PDBClean configuration loading."""

from pathlib import Path

import pytest

from pdbclean.config import ConfigError, load_config


CONFIG_PATH = Path(
    "config/pdbclean/protocol_3_2_comp702_v1.yaml"
)


def test_load_valid_project_config(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setenv("TMPDIR", "/tmp/test-user")

    loaded = load_config(CONFIG_PATH)

    assert loaded.data["release"]["snapshot"] == "20260101"
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


def test_snapshot_prefix_must_match_release(
    tmp_path: Path,
) -> None:
    original = CONFIG_PATH.read_text()
    modified = original.replace(
        'snapshot: "20260101"',
        'snapshot: "20250101"',
        1,
    )

    config_path = tmp_path / "mismatch.yaml"
    config_path.write_text(modified)

    with pytest.raises(
        ConfigError,
        match="release.snapshot does not match",
    ):
        load_config(config_path)
