from pdbclean.acta_downstream_investigation import (
    _output_directory_name,
    _passes_resolution,
    _resolution_for_entry,
    _validate_downstream_config,
)


def _config():
    return {
        "stage": {
            "name": "acta_downstream_investigation",
            "version": "1.0",
        },
        "resolution_filter": {
            "enabled": True,
            "maximum_angstrom": 4.0,
            "operator": "less_than_or_equal",
            "missing_resolution_policy": "reject",
            "method_precedence": [
                "electron_microscopy",
                "diffraction",
            ],
            "electron_microscopy": {
                "experimental_method": (
                    "ELECTRON MICROSCOPY"
                ),
                "source_field": (
                    "em_3d_reconstruction_resolution"
                ),
                "aggregation": "minimum_numeric",
            },
            "diffraction": {
                "experimental_methods": [
                    "X-RAY DIFFRACTION",
                    "NEUTRON DIFFRACTION",
                    "FIBER DIFFRACTION",
                    "ELECTRON CRYSTALLOGRAPHY",
                ],
                "source_field": (
                    "refine_ls_d_res_high"
                ),
                "aggregation": "minimum_numeric",
            },
        },
        "pandda_group_deposition_filter": {
            "enabled": True,
            "metadata_field": (
                "deposit_group_mentions_pandda"
            ),
            "reject_when": True,
        },
        "same_deposition_filter": {
            "enabled": True,
            "query_field": "query_pdb_id",
            "subject_field": "subject_pdb_id",
            "reject_when_equal": True,
        },
        "manual_review": {
            "automatic_virus_filtering": False,
            "automatic_ribosome_filtering": False,
        },
    }


def test_xray_resolution():
    config = _config()

    row = {
        "experimental_methods": [
            "X-RAY DIFFRACTION"
        ],
        "refine_ls_d_res_high": [2.2],
        "em_3d_reconstruction_resolution": [],
    }

    value, basis = _resolution_for_entry(
        row,
        config["resolution_filter"],
    )

    assert value == 2.2
    assert basis == "diffraction"
    assert _passes_resolution(
        value,
        config["resolution_filter"],
    )


def test_em_precedes_refine():
    config = _config()

    row = {
        "experimental_methods": [
            "ELECTRON MICROSCOPY"
        ],
        "refine_ls_d_res_high": [1.5],
        "em_3d_reconstruction_resolution": [4.5],
    }

    value, basis = _resolution_for_entry(
        row,
        config["resolution_filter"],
    )

    assert value == 4.5
    assert basis == "electron_microscopy"
    assert not _passes_resolution(
        value,
        config["resolution_filter"],
    )


def test_em_uses_minimum_numeric():
    config = _config()

    row = {
        "experimental_methods": [
            "ELECTRON MICROSCOPY"
        ],
        "refine_ls_d_res_high": [],
        "em_3d_reconstruction_resolution": [
            3.9,
            3.3,
        ],
    }

    value, basis = _resolution_for_entry(
        row,
        config["resolution_filter"],
    )

    assert value == 3.3
    assert basis == "electron_microscopy"


def test_missing_resolution_rejected():
    config = _config()

    row = {
        "experimental_methods": [
            "SOLUTION NMR"
        ],
        "refine_ls_d_res_high": [],
        "em_3d_reconstruction_resolution": [],
    }

    value, basis = _resolution_for_entry(
        row,
        config["resolution_filter"],
    )

    assert value is None
    assert basis == "no_resolution_method"
    assert not _passes_resolution(
        value,
        config["resolution_filter"],
    )


def test_exact_four_angstrom_passes():
    config = _config()

    assert _passes_resolution(
        4.0,
        config["resolution_filter"],
    )


def test_downstream_config_validation():
    _validate_downstream_config(
        _config()
    )



def test_v1_publication_directory_is_preserved():
    config = _config()

    assert (
        _output_directory_name(config)
        == "acta_downstream_investigation"
    )


def test_v2_publication_directory_is_isolated():
    config = _config()
    config["stage"]["version"] = "2.0"

    _validate_downstream_config(config)

    assert (
        _output_directory_name(config)
        == "acta_downstream_investigation_v2"
    )



def test_v2_config_is_xray_crystallographic_only():
    from pathlib import Path

    import yaml

    path = Path(
        "config/pdbclean/"
        "acta_downstream_investigation_v2.yaml"
    )

    config = yaml.safe_load(
        path.read_text(encoding="utf-8")
    )

    _validate_downstream_config(config)

    assert config["stage"]["version"] == "2.0"

    assert (
        _output_directory_name(config)
        == "acta_downstream_investigation_v2"
    )

    resolution = config["resolution_filter"]

    assert resolution["method_precedence"] == [
        "diffraction"
    ]

    assert resolution["diffraction"]["experimental_methods"] == [
        "X-RAY DIFFRACTION"
    ]

    assert (
        resolution["diffraction"]["source_field"]
        == "refine_ls_d_res_high"
    )

    assert resolution["maximum_angstrom"] == 4.0
    assert resolution["operator"] == "less_than_or_equal"


def test_v2_accepts_xray_refine_resolution():
    from pathlib import Path

    import yaml

    config = yaml.safe_load(
        Path(
            "config/pdbclean/"
            "acta_downstream_investigation_v2.yaml"
        ).read_text(encoding="utf-8")
    )

    row = {
        "experimental_methods": [
            "X-RAY DIFFRACTION"
        ],
        "refine_ls_d_res_high": [2.2],
        "em_3d_reconstruction_resolution": [],
    }

    value, basis = _resolution_for_entry(
        row,
        config["resolution_filter"],
    )

    assert value == 2.2
    assert basis == "diffraction"
    assert _passes_resolution(
        value,
        config["resolution_filter"],
    )


def test_v2_rejects_electron_microscopy_scope():
    from pathlib import Path

    import yaml

    config = yaml.safe_load(
        Path(
            "config/pdbclean/"
            "acta_downstream_investigation_v2.yaml"
        ).read_text(encoding="utf-8")
    )

    row = {
        "experimental_methods": [
            "ELECTRON MICROSCOPY"
        ],
        "refine_ls_d_res_high": [2.0],
        "em_3d_reconstruction_resolution": [2.0],
    }

    value, basis = _resolution_for_entry(
        row,
        config["resolution_filter"],
    )

    assert value is None
    assert basis == "no_resolution_method"
    assert not _passes_resolution(
        value,
        config["resolution_filter"],
    )


def test_v2_preflight_oracle_is_frozen():
    from pathlib import Path

    import yaml

    config = yaml.safe_load(
        Path(
            "config/pdbclean/"
            "acta_downstream_investigation_v2.yaml"
        ).read_text(encoding="utf-8")
    )

    assert config["expected_preflight_counts"] == {
        "input_geometric_near_duplicate_pairs": 1072751,
        "after_resolution_filter": 95297,
        "after_pandda_group_deposition_filter": 95285,
        "after_same_deposition_filter": 11572,
    }
