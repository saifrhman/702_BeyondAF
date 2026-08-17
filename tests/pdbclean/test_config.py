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

    assert loaded.data["snapshot"]["mode"] == "latest_complete"
    assert "snapshot_id" not in loaded.data["snapshot"]
    assert "expected_mmcif_count" not in loaded.data["snapshot"]
    assert "expected_total_bytes" not in loaded.data["snapshot"]
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
    data = yaml.safe_load(CONFIG_PATH.read_text())
    data["snapshot"]["mode"] = "invalid"

    config_path = tmp_path / "invalid_mode.yaml"
    config_path.write_text(
        yaml.safe_dump(data, sort_keys=False)
    )

    with pytest.raises(
        ConfigError,
        match="snapshot.mode must be either",
    ):
        load_config(config_path)


def test_fixed_mode_requires_snapshot_id(
    tmp_path: Path,
) -> None:
    data = yaml.safe_load(CONFIG_PATH.read_text())
    data["snapshot"]["mode"] = "fixed"
    data["snapshot"].pop("snapshot_id", None)

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
    data["snapshot"]["snapshot_id"] = "20260101"

    config_path = tmp_path / "latest_with_id.yaml"
    config_path.write_text(
        yaml.safe_dump(data, sort_keys=False)
    )

    with pytest.raises(
        ConfigError,
        match="latest_complete mode must not define snapshot_id",
    ):
        load_config(config_path)


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("expected_mmcif_count", 246905),
        ("expected_total_bytes", 85079649893),
    ],
)
def test_latest_complete_rejects_snapshot_specific_expectations(
    tmp_path: Path,
    field: str,
    value: int,
) -> None:
    data = yaml.safe_load(CONFIG_PATH.read_text())
    data["snapshot"]["mode"] = "latest_complete"
    data["snapshot"].pop("snapshot_id", None)
    data["snapshot"][field] = value

    config_path = tmp_path / f"latest_with_{field}.yaml"
    config_path.write_text(
        yaml.safe_dump(data, sort_keys=False)
    )

    with pytest.raises(
        ConfigError,
        match=(
            "latest_complete mode must not define "
            + field
        ),
    ):
        load_config(config_path)


@pytest.mark.parametrize(
    "model_id",
    [0, -1, True, 1.5, "1"],
)
def test_invalid_model_id_is_rejected(
    tmp_path: Path,
    model_id,
) -> None:
    data = yaml.safe_load(CONFIG_PATH.read_text())
    data["selection"]["models"]["model_id"] = model_id

    config_path = tmp_path / "invalid_model_id.yaml"
    config_path.write_text(
        yaml.safe_dump(data, sort_keys=False)
    )

    with pytest.raises(
        ConfigError,
        match="selection.models.model_id must be a positive integer",
    ):
        load_config(config_path)


def test_all_models_policy_is_accepted(
    tmp_path: Path,
) -> None:
    data = yaml.safe_load(CONFIG_PATH.read_text())
    data["selection"]["models"] = {
        "policy": "all_models",
    }

    config_path = tmp_path / "all_models.yaml"
    config_path.write_text(
        yaml.safe_dump(data, sort_keys=False)
    )

    load_config(config_path)


def test_unsupported_model_policy_is_rejected(
    tmp_path: Path,
) -> None:
    data = yaml.safe_load(CONFIG_PATH.read_text())
    data["selection"]["models"]["policy"] = "unsupported"

    config_path = tmp_path / "invalid_model_policy.yaml"
    config_path.write_text(
        yaml.safe_dump(data, sort_keys=False)
    )

    with pytest.raises(
        ConfigError,
        match="selection.models.policy",
    ):
        load_config(config_path)



@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("batch_size", 0),
        ("batch_size", -1),
        ("batch_size", True),
        ("download_concurrency", 0),
        ("quality_array_worker_count", 0),
        ("quality_array_worker_count", -1),
        ("quality_array_worker_count", True),
        ("quality_array_concurrency", 0),
        ("quality_array_concurrency", -1),
        ("quality_array_concurrency", True),
        ("connection_timeout_seconds", 0),
    ],
)
def test_invalid_positive_execution_integer_is_rejected(
    tmp_path: Path,
    field: str,
    value,
) -> None:
    data = yaml.safe_load(CONFIG_PATH.read_text())
    data["execution"][field] = value

    config_path = tmp_path / f"invalid_{field}.yaml"
    config_path.write_text(
        yaml.safe_dump(data, sort_keys=False)
    )

    with pytest.raises(
        ConfigError,
        match=rf"execution\.{field} must be a positive integer",
    ):
        load_config(config_path)


@pytest.mark.parametrize(
    "value",
    [-1, True, 1.5, "3"],
)
def test_invalid_max_retries_is_rejected(
    tmp_path: Path,
    value,
) -> None:
    data = yaml.safe_load(CONFIG_PATH.read_text())
    data["execution"]["max_retries"] = value

    config_path = tmp_path / "invalid_max_retries.yaml"
    config_path.write_text(
        yaml.safe_dump(data, sort_keys=False)
    )

    with pytest.raises(
        ConfigError,
        match="execution.max_retries must be a non-negative integer",
    ):
        load_config(config_path)


@pytest.mark.parametrize(
    "field",
    [
        "atomic_writes",
        "write_success_markers",
    ],
)
def test_execution_boolean_fields_are_validated(
    tmp_path: Path,
    field: str,
) -> None:
    data = yaml.safe_load(CONFIG_PATH.read_text())
    data["execution"][field] = "true"

    config_path = tmp_path / f"invalid_{field}.yaml"
    config_path.write_text(
        yaml.safe_dump(data, sort_keys=False)
    )

    with pytest.raises(
        ConfigError,
        match=rf"execution\.{field} must be a boolean",
    ):
        load_config(config_path)


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("dataset_name", ""),
        ("dataset_name", None),
        ("protocol_version", ""),
        ("protocol_version", None),
    ],
)
def test_invalid_release_provenance_is_rejected(
    tmp_path: Path,
    field: str,
    value,
) -> None:
    data = yaml.safe_load(CONFIG_PATH.read_text())
    data["release"][field] = value

    config_path = tmp_path / f"invalid_release_{field}.yaml"
    config_path.write_text(
        yaml.safe_dump(data, sort_keys=False)
    )

    with pytest.raises(
        ConfigError,
        match=rf"release\.{field} must be a non-empty string",
    ):
        load_config(config_path)


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("temporary_root", ""),
        ("temporary_root", None),
        ("output_root", ""),
        ("output_root", None),
    ],
)
def test_invalid_storage_path_is_rejected(
    tmp_path: Path,
    field: str,
    value,
) -> None:
    data = yaml.safe_load(CONFIG_PATH.read_text())
    data["storage"][field] = value

    config_path = tmp_path / f"invalid_storage_{field}.yaml"
    config_path.write_text(
        yaml.safe_dump(data, sort_keys=False)
    )

    with pytest.raises(
        ConfigError,
        match=rf"storage\.{field} must be a non-empty string",
    ):
        load_config(config_path)


@pytest.mark.parametrize(
    "value",
    [None, "false", 0, 1],
)
def test_retain_downloaded_mmcif_must_be_boolean(
    tmp_path: Path,
    value,
) -> None:
    data = yaml.safe_load(CONFIG_PATH.read_text())
    data["storage"]["retain_downloaded_mmcif"] = value

    config_path = tmp_path / "invalid_retain_downloaded.yaml"
    config_path.write_text(
        yaml.safe_dump(data, sort_keys=False)
    )

    with pytest.raises(
        ConfigError,
        match="storage.retain_downloaded_mmcif must be a boolean",
    ):
        load_config(config_path)


def test_post_cleaning_geometric_validation_config_is_loaded(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setenv("TMPDIR", "/tmp/test-user")

    loaded = load_config(CONFIG_PATH)

    assert (
        loaded.data[
            "post_cleaning_geometric_validation"
        ]["enabled"]
        is True
    )

    assert (
        loaded.data[
            "post_cleaning_geometric_validation"
        ]["minimum_triangle_angle_degrees"]
        == 3.0
    )

    assert (
        loaded.data["quality_rules"]
        ["backbone_distance"]
        ["minimum_distance_angstrom"]
        == 0.01
    )


@pytest.mark.parametrize(
    "value",
    [-1, 181, True, "3.0", None],
)
def test_invalid_post_cleaning_triangle_angle_is_rejected(
    tmp_path: Path,
    value,
) -> None:
    data = yaml.safe_load(CONFIG_PATH.read_text())

    data[
        "post_cleaning_geometric_validation"
    ]["minimum_triangle_angle_degrees"] = value

    config_path = tmp_path / "invalid_geometry_angle.yaml"
    config_path.write_text(
        yaml.safe_dump(data, sort_keys=False)
    )

    with pytest.raises(
        ConfigError,
        match="minimum_triangle_angle_degrees",
    ):
        load_config(config_path)


@pytest.mark.parametrize(
    "value",
    [-0.001, True, "0.01", None],
)
def test_invalid_backbone_distance_threshold_is_rejected(
    tmp_path: Path,
    value,
) -> None:
    data = yaml.safe_load(CONFIG_PATH.read_text())

    data["quality_rules"]["backbone_distance"][
        "minimum_distance_angstrom"
    ] = value

    config_path = tmp_path / "invalid_backbone_distance.yaml"
    config_path.write_text(
        yaml.safe_dump(data, sort_keys=False)
    )

    with pytest.raises(
        ConfigError,
        match="minimum_distance_angstrom",
    ):
        load_config(config_path)


def test_geometric_thresholds_can_be_changed_in_config(
    tmp_path: Path,
) -> None:
    data = yaml.safe_load(CONFIG_PATH.read_text())

    data["quality_rules"]["backbone_distance"][
        "minimum_distance_angstrom"
    ] = 0.02

    data[
        "post_cleaning_geometric_validation"
    ]["minimum_triangle_angle_degrees"] = 5.0

    config_path = tmp_path / "changed_geometry_thresholds.yaml"
    config_path.write_text(
        yaml.safe_dump(data, sort_keys=False)
    )

    loaded = load_config(config_path)

    assert (
        loaded.data["quality_rules"]
        ["backbone_distance"]
        ["minimum_distance_angstrom"]
        == 0.02
    )

    assert (
        loaded.data[
            "post_cleaning_geometric_validation"
        ]["minimum_triangle_angle_degrees"]
        == 5.0
    )



def test_duplicate_search_threshold_is_loaded() -> None:
    loaded = load_config(CONFIG_PATH)

    assert (
        loaded.data["duplicate_search"]
        ["near_duplicate_threshold_angstrom"]
        == 0.010
    )


@pytest.mark.parametrize(
    "value",
    [
        -0.001,
        0.0,
        True,
        "0.010",
        None,
        0.0105,
        float("inf"),
        float("nan"),
    ],
)
def test_invalid_duplicate_search_threshold_is_rejected(
    tmp_path: Path,
    value,
) -> None:
    data = yaml.safe_load(
        CONFIG_PATH.read_text()
    )

    data["duplicate_search"][
        "near_duplicate_threshold_angstrom"
    ] = value

    config_path = (
        tmp_path
        / "invalid_duplicate_search_threshold.yaml"
    )

    config_path.write_text(
        yaml.safe_dump(
            data,
            sort_keys=False,
        )
    )

    with pytest.raises(
        ConfigError,
        match="near_duplicate_threshold_angstrom",
    ):
        load_config(config_path)
