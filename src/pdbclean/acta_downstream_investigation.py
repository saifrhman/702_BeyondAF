"""Acta paper-derived downstream investigation filtering.

This stage is deliberately downstream of geometric duplicate detection.

Scientific sequence
-------------------
Stage-8 complete-BRI near-duplicates
    -> both depositions have experimental resolution <= configured cutoff
    -> reject PanDDA Group-deposition entries
    -> reject hits involving two chains from the same PDB deposition
    -> detailed-inspection candidate pool

Virus/ribosome exclusions are NOT automated here.  In the Acta paper
those belong to the subsequent manual/detailed-evaluation stage.

This module never modifies Stage-8 or downstream-metadata publications.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path
from typing import Any

import pyarrow as pa
import pyarrow.parquet as pq
import yaml

from pdbclean.config import load_config


class ActaDownstreamInvestigationError(RuntimeError):
    """Raised when downstream filtering cannot be performed safely."""


SUMMARY_SCHEMA_NAME = (
    "pdbclean_acta_downstream_investigation_global_summary"
)
SUMMARY_SCHEMA_VERSION = "1.0"

SUCCESS_SCHEMA_NAME = (
    "pdbclean_acta_downstream_investigation_success"
)
SUCCESS_SCHEMA_VERSION = "1.0"


STAGE_OUTPUT_DIRECTORIES = {
    "1.0": "acta_downstream_investigation",
    "2.0": "acta_downstream_investigation_v2",
}


OUTPUT_SCHEMA = pa.schema(
    [
        pa.field("query_snapshot", pa.string(), nullable=False),
        pa.field("query_pdb_id", pa.string(), nullable=False),
        pa.field("query_model_id", pa.int64(), nullable=False),
        pa.field("query_label_chain_id", pa.string(), nullable=False),
        pa.field("subject_snapshot", pa.string(), nullable=False),
        pa.field("subject_pdb_id", pa.string(), nullable=False),
        pa.field("subject_model_id", pa.int64(), nullable=False),
        pa.field(
            "subject_label_chain_id",
            pa.string(),
            nullable=False,
        ),
        pa.field(
            "retained_residue_count",
            pa.int64(),
            nullable=False,
        ),
        pa.field(
            "d_brain_numerator_max",
            pa.int64(),
            nullable=True,
        ),
        pa.field(
            "d_brain",
            pa.float64(),
            nullable=True,
        ),
        pa.field("d_bri_mA", pa.int64(), nullable=False),
        pa.field("d_bri", pa.float64(), nullable=False),
        pa.field(
            "stage8_source",
            pa.string(),
            nullable=False,
        ),
        pa.field(
            "query_resolution_angstrom",
            pa.float64(),
            nullable=False,
        ),
        pa.field(
            "query_resolution_basis",
            pa.string(),
            nullable=False,
        ),
        pa.field(
            "subject_resolution_angstrom",
            pa.float64(),
            nullable=False,
        ),
        pa.field(
            "subject_resolution_basis",
            pa.string(),
            nullable=False,
        ),
    ],
    metadata={
        b"schema_name": (
            b"pdbclean_acta_detailed_inspection_candidates"
        ),
        b"schema_version": b"1.0",
        b"scientific_sequence": (
            b"resolution_then_pandda_then_same_deposition"
        ),
        b"geometric_input": (
            b"Stage-8 complete-BRI near-duplicates"
        ),
    },
)


def _validate_commit(value: str) -> str:
    value = value.strip().lower()

    if re.fullmatch(r"[0-9a-f]{40}", value) is None:
        raise ActaDownstreamInvestigationError(
            f"Expected full 40-character Git commit, found {value!r}"
        )

    return value


def _read_json(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise ActaDownstreamInvestigationError(
            f"Required JSON file is missing: {path}"
        )

    with path.open(
        "r",
        encoding="utf-8",
    ) as handle:
        value = json.load(handle)

    if not isinstance(value, dict):
        raise ActaDownstreamInvestigationError(
            f"Expected JSON object in {path}"
        )

    return value


def _write_json_atomic(
    data: dict[str, Any],
    path: Path,
) -> None:
    path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    temporary = path.with_suffix(
        path.suffix + ".tmp"
    )

    temporary.write_text(
        json.dumps(
            data,
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    temporary.replace(path)


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        while True:
            chunk = handle.read(1024 * 1024)

            if not chunk:
                break

            digest.update(chunk)

    return digest.hexdigest()


def _load_downstream_config(
    path: Path,
) -> tuple[dict[str, Any], str]:
    if not path.is_file():
        raise ActaDownstreamInvestigationError(
            f"Downstream config is missing: {path}"
        )

    raw = path.read_text(
        encoding="utf-8"
    )

    value = yaml.safe_load(raw)

    if not isinstance(value, dict):
        raise ActaDownstreamInvestigationError(
            "Downstream config must contain a YAML mapping"
        )

    return value, _sha256_file(path)


def _resolution_for_entry(
    row: dict[str, Any],
    resolution_config: dict[str, Any],
) -> tuple[float | None, str]:
    """Return the operational experimental resolution and basis.

    Method precedence is explicit in the downstream configuration.
    Missing or unsupported resolution returns ``(None, basis)``.
    """

    methods = {
        str(value)
        for value in row[
            "experimental_methods"
        ]
    }

    precedence = list(
        resolution_config[
            "method_precedence"
        ]
    )

    for basis in precedence:
        if basis == "electron_microscopy":
            rule = resolution_config[
                "electron_microscopy"
            ]

            method = str(
                rule[
                    "experimental_method"
                ]
            )

            if method not in methods:
                continue

            values = [
                float(value)
                for value in row[
                    str(
                        rule[
                            "source_field"
                        ]
                    )
                ]
            ]

            if not values:
                return (
                    None,
                    "electron_microscopy_missing",
                )

            if (
                rule.get("aggregation")
                != "minimum_numeric"
            ):
                raise ActaDownstreamInvestigationError(
                    "Unsupported EM resolution aggregation"
                )

            return (
                min(values),
                "electron_microscopy",
            )

        if basis == "diffraction":
            rule = resolution_config[
                "diffraction"
            ]

            applicable = {
                str(value)
                for value in rule[
                    "experimental_methods"
                ]
            }

            if not (
                methods & applicable
            ):
                continue

            values = [
                float(value)
                for value in row[
                    str(
                        rule[
                            "source_field"
                        ]
                    )
                ]
            ]

            if not values:
                return (
                    None,
                    "diffraction_missing",
                )

            if (
                rule.get("aggregation")
                != "minimum_numeric"
            ):
                raise ActaDownstreamInvestigationError(
                    "Unsupported diffraction resolution aggregation"
                )

            return (
                min(values),
                "diffraction",
            )

        raise ActaDownstreamInvestigationError(
            f"Unsupported resolution method basis: {basis!r}"
        )

    return (
        None,
        "no_resolution_method",
    )


def _passes_resolution(
    resolution: float | None,
    resolution_config: dict[str, Any],
) -> bool:
    if resolution is None:
        if (
            resolution_config[
                "missing_resolution_policy"
            ]
            != "reject"
        ):
            raise ActaDownstreamInvestigationError(
                "Unsupported missing-resolution policy"
            )

        return False

    if (
        resolution_config[
            "operator"
        ]
        != "less_than_or_equal"
    ):
        raise ActaDownstreamInvestigationError(
            "Unsupported resolution comparison operator"
        )

    threshold = float(
        resolution_config[
            "maximum_angstrom"
        ]
    )

    return resolution <= threshold


def _stage_version(
    config: dict[str, Any],
) -> str:
    stage = config.get("stage", {})

    if (
        stage.get("name")
        != "acta_downstream_investigation"
    ):
        raise ActaDownstreamInvestigationError(
            "Unexpected downstream stage identity"
        )

    version = str(
        stage.get("version")
    )

    if version not in STAGE_OUTPUT_DIRECTORIES:
        raise ActaDownstreamInvestigationError(
            f"Unsupported downstream stage version: {version!r}"
        )

    return version


def _output_directory_name(
    config: dict[str, Any],
) -> str:
    return STAGE_OUTPUT_DIRECTORIES[
        _stage_version(config)
    ]


def _validate_downstream_config(
    config: dict[str, Any],
) -> None:
    _stage_version(config)

    resolution = config[
        "resolution_filter"
    ]

    if not resolution.get(
        "enabled"
    ):
        raise ActaDownstreamInvestigationError(
            "Resolution filter must be enabled"
        )

    if float(
        resolution[
            "maximum_angstrom"
        ]
    ) <= 0.0:
        raise ActaDownstreamInvestigationError(
            "Resolution threshold must be positive"
        )

    pandda = config[
        "pandda_group_deposition_filter"
    ]

    if (
        not pandda.get("enabled")
        or pandda.get("metadata_field")
        != "deposit_group_mentions_pandda"
        or pandda.get("reject_when") is not True
    ):
        raise ActaDownstreamInvestigationError(
            "Unexpected PanDDA filter definition"
        )

    same = config[
        "same_deposition_filter"
    ]

    if (
        not same.get("enabled")
        or same.get("query_field")
        != "query_pdb_id"
        or same.get("subject_field")
        != "subject_pdb_id"
        or same.get("reject_when_equal") is not True
    ):
        raise ActaDownstreamInvestigationError(
            "Unexpected same-deposition filter definition"
        )

    manual = config[
        "manual_review"
    ]

    if (
        manual.get(
            "automatic_virus_filtering"
        )
        is not False
        or manual.get(
            "automatic_ribosome_filtering"
        )
        is not False
    ):
        raise ActaDownstreamInvestigationError(
            "Virus/ribosome exclusions must remain manual"
        )


def _expected_preflight_counts(
    config: dict[str, Any],
) -> dict[str, int]:
    raw = config.get(
        "expected_preflight_counts",
        {},
    )

    required = [
        "input_geometric_near_duplicate_pairs",
        "after_resolution_filter",
        "after_pandda_group_deposition_filter",
        "after_same_deposition_filter",
    ]

    result: dict[str, int] = {}

    for key in required:
        if key not in raw:
            raise ActaDownstreamInvestigationError(
                f"Missing preflight oracle count: {key}"
            )

        result[key] = int(
            raw[key]
        )

    return result


def _validate_upstream(
    *,
    canonical_config_hash: str,
    stage8_success: dict[str, Any],
    stage8_summary: dict[str, Any],
    metadata_success: dict[str, Any],
    metadata_summary: dict[str, Any],
) -> None:
    if (
        stage8_success.get(
            "success_schema_name"
        )
        != "pdbclean_stage8_full_bri_nn_success"
    ):
        raise ActaDownstreamInvestigationError(
            "Unexpected Stage-8 success schema"
        )

    if (
        stage8_summary.get(
            "summary_schema_name"
        )
        != "pdbclean_stage8_full_bri_nn_global_summary"
    ):
        raise ActaDownstreamInvestigationError(
            "Unexpected Stage-8 summary schema"
        )

    if not (
        stage8_summary.get(
            "exhaustive_oracle_pair_set_equal"
        )
        and stage8_summary.get(
            "exhaustive_oracle_distances_equal"
        )
        and stage8_summary.get(
            "pair_accounting_valid"
        )
        and stage8_summary.get(
            "processing_error_count"
        )
        == 0
    ):
        raise ActaDownstreamInvestigationError(
            "Stage-8 validation gate failed"
        )

    if (
        stage8_summary.get(
            "classification_performed"
        )
        is not False
    ):
        raise ActaDownstreamInvestigationError(
            "Stage-8 input unexpectedly performed classification"
        )

    if (
        metadata_success.get(
            "success_schema_name"
        )
        != "pdbclean_downstream_metadata_success"
    ):
        raise ActaDownstreamInvestigationError(
            "Unexpected downstream-metadata success schema"
        )

    if (
        metadata_summary.get(
            "summary_schema_name"
        )
        != "pdbclean_downstream_metadata_global_summary"
    ):
        raise ActaDownstreamInvestigationError(
            "Unexpected downstream-metadata summary schema"
        )

    if (
        metadata_summary.get(
            "scientific_filtering_performed"
        )
        is not False
        or metadata_summary.get(
            "processing_error_count"
        )
        != 0
    ):
        raise ActaDownstreamInvestigationError(
            "Metadata validation gate failed"
        )

    for source in (
        stage8_success,
        stage8_summary,
        metadata_success,
        metadata_summary,
    ):
        if (
            source.get(
                "config_sha256"
            )
            != canonical_config_hash
        ):
            raise ActaDownstreamInvestigationError(
                "Canonical configuration provenance mismatch"
            )

    stage8_commit = stage8_success[
        "full_bri_nn_pipeline_git_commit"
    ]

    if (
        stage8_summary.get(
            "full_bri_nn_pipeline_git_commit"
        )
        != stage8_commit
        or metadata_success.get(
            "full_bri_nn_pipeline_git_commit"
        )
        != stage8_commit
        or metadata_summary.get(
            "full_bri_nn_pipeline_git_commit"
        )
        != stage8_commit
    ):
        raise ActaDownstreamInvestigationError(
            "Stage-8 provenance linkage mismatch"
        )

    snapshots = {
        str(
            stage8_success[
                "snapshot"
            ]
        ),
        str(
            stage8_summary[
                "snapshot"
            ]
        ),
        str(
            metadata_success[
                "snapshot"
            ]
        ),
        str(
            metadata_summary[
                "snapshot"
            ]
        ),
    }

    if len(snapshots) != 1:
        raise ActaDownstreamInvestigationError(
            "Upstream snapshot mismatch"
        )


def _metadata_map(
    path: Path,
    *,
    snapshot: str,
    resolution_config: dict[str, Any],
) -> dict[str, dict[str, Any]]:
    table = pq.read_table(path)

    metadata: dict[
        str,
        dict[str, Any],
    ] = {}

    for row in table.to_pylist():
        if str(
            row["snapshot"]
        ) != snapshot:
            raise ActaDownstreamInvestigationError(
                "Metadata snapshot mismatch"
            )

        pdb_id = str(
            row["pdb_id"]
        ).lower()

        if pdb_id in metadata:
            raise ActaDownstreamInvestigationError(
                f"Duplicate metadata row for {pdb_id}"
            )

        resolution, basis = _resolution_for_entry(
            row,
            resolution_config,
        )

        metadata[pdb_id] = {
            "resolution": resolution,
            "resolution_basis": basis,
            "passes_resolution": (
                _passes_resolution(
                    resolution,
                    resolution_config,
                )
            ),
            "pandda": bool(
                row[
                    "deposit_group_mentions_pandda"
                ]
            ),
        }

    return metadata


def _process_pair_file(
    path: Path,
    *,
    source_name: str,
    metadata: dict[str, dict[str, Any]],
    writer: pq.ParquetWriter,
    snapshot: str,
    threshold_mA: int,
    counters: dict[str, int],
    final_pdb_ids: set[str],
) -> None:
    if not path.is_file():
        raise ActaDownstreamInvestigationError(
            f"Stage-8 pair file is missing: {path}"
        )

    pf = pq.ParquetFile(path)

    available = set(
        pf.schema_arrow.names
    )

    required = {
        "query_snapshot",
        "query_pdb_id",
        "query_model_id",
        "query_label_chain_id",
        "subject_snapshot",
        "subject_pdb_id",
        "subject_model_id",
        "subject_label_chain_id",
        "retained_residue_count",
        "d_bri_mA",
        "d_bri",
    }

    if not required.issubset(
        available
    ):
        raise ActaDownstreamInvestigationError(
            f"Unexpected Stage-8 schema in {path}"
        )

    has_brain = {
        "d_brain_numerator_max",
        "d_brain",
    }.issubset(available)

    columns = [
        "query_snapshot",
        "query_pdb_id",
        "query_model_id",
        "query_label_chain_id",
        "subject_snapshot",
        "subject_pdb_id",
        "subject_model_id",
        "subject_label_chain_id",
        "retained_residue_count",
        "d_bri_mA",
        "d_bri",
    ]

    if has_brain:
        columns.extend(
            [
                "d_brain_numerator_max",
                "d_brain",
            ]
        )

    for batch in pf.iter_batches(
        columns=columns,
        batch_size=65_536,
    ):
        retained_rows = []

        for row in pa.Table.from_batches(
            [batch]
        ).to_pylist():
            counters["input"] += 1

            q_snapshot = str(
                row["query_snapshot"]
            )
            s_snapshot = str(
                row["subject_snapshot"]
            )

            if (
                q_snapshot != snapshot
                or s_snapshot != snapshot
            ):
                raise ActaDownstreamInvestigationError(
                    "Stage-8 pair snapshot mismatch"
                )

            distance_mA = int(
                row["d_bri_mA"]
            )

            if (
                distance_mA < 0
                or distance_mA > threshold_mA
            ):
                raise ActaDownstreamInvestigationError(
                    "Stage-8 thresholded publication contains "
                    "an out-of-range complete-BRI distance"
                )

            q = str(
                row["query_pdb_id"]
            ).lower()
            s = str(
                row["subject_pdb_id"]
            ).lower()

            try:
                q_meta = metadata[q]
                s_meta = metadata[s]
            except KeyError as exc:
                raise ActaDownstreamInvestigationError(
                    f"Missing metadata for Stage-8 pair {q}, {s}"
                ) from exc

            # Paper-derived downstream step 1:
            # both structures have resolution cutoff or better.
            if not (
                q_meta[
                    "passes_resolution"
                ]
                and s_meta[
                    "passes_resolution"
                ]
            ):
                counters[
                    "resolution_rejected"
                ] += 1
                continue

            counters[
                "after_resolution"
            ] += 1

            # Paper-derived downstream step 2:
            # reject PanDDA Group-deposition entries.
            if (
                q_meta["pandda"]
                or s_meta["pandda"]
            ):
                counters[
                    "pandda_rejected"
                ] += 1
                continue

            counters[
                "after_pandda"
            ] += 1

            # Paper-derived downstream step 3:
            # reject two parts of the same deposition.
            if q == s:
                counters[
                    "same_deposition_rejected"
                ] += 1
                continue

            counters[
                "after_same_deposition"
            ] += 1

            final_pdb_ids.update(
                (q, s)
            )

            retained_rows.append(
                {
                    "query_snapshot": q_snapshot,
                    "query_pdb_id": row[
                        "query_pdb_id"
                    ],
                    "query_model_id": int(
                        row[
                            "query_model_id"
                        ]
                    ),
                    "query_label_chain_id": row[
                        "query_label_chain_id"
                    ],
                    "subject_snapshot": s_snapshot,
                    "subject_pdb_id": row[
                        "subject_pdb_id"
                    ],
                    "subject_model_id": int(
                        row[
                            "subject_model_id"
                        ]
                    ),
                    "subject_label_chain_id": row[
                        "subject_label_chain_id"
                    ],
                    "retained_residue_count": int(
                        row[
                            "retained_residue_count"
                        ]
                    ),
                    "d_brain_numerator_max": (
                        int(
                            row[
                                "d_brain_numerator_max"
                            ]
                        )
                        if has_brain
                        else None
                    ),
                    "d_brain": (
                        float(
                            row[
                                "d_brain"
                            ]
                        )
                        if has_brain
                        else None
                    ),
                    "d_bri_mA": distance_mA,
                    "d_bri": float(
                        row[
                            "d_bri"
                        ]
                    ),
                    "stage8_source": (
                        source_name
                    ),
                    "query_resolution_angstrom": float(
                        q_meta[
                            "resolution"
                        ]
                    ),
                    "query_resolution_basis": str(
                        q_meta[
                            "resolution_basis"
                        ]
                    ),
                    "subject_resolution_angstrom": float(
                        s_meta[
                            "resolution"
                        ]
                    ),
                    "subject_resolution_basis": str(
                        s_meta[
                            "resolution_basis"
                        ]
                    ),
                }
            )

        if retained_rows:
            writer.write_table(
                pa.Table.from_pylist(
                    retained_rows,
                    schema=OUTPUT_SCHEMA,
                )
            )


def main() -> int:
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--config",
        required=True,
        type=Path,
        help="Frozen canonical Protocol 3.2 configuration",
    )
    parser.add_argument(
        "--downstream-config",
        required=True,
        type=Path,
        help="Acta downstream investigation configuration",
    )
    parser.add_argument(
        "--pipeline-git-commit",
        required=True,
    )

    args = parser.parse_args()

    pipeline_commit = _validate_commit(
        args.pipeline_git_commit
    )

    canonical_loaded = load_config(
        args.config
    )
    canonical = canonical_loaded.data

    downstream, downstream_hash = (
        _load_downstream_config(
            args.downstream_config
        )
    )

    _validate_downstream_config(
        downstream
    )

    expected = (
        _expected_preflight_counts(
            downstream
        )
    )

    repo = Path.cwd()

    protocol = canonical[
        "release"
    ]["protocol_version"]

    storage_root = Path(
        canonical[
            "storage"
        ]["output_root"]
    )

    if not storage_root.is_absolute():
        storage_root = (
            repo / storage_root
        )

    stage8_roots = sorted(
        path.parent
        for path in storage_root.glob(
            f"*/{protocol}/full_bri_nn/_SUCCESS"
        )
    )

    if len(stage8_roots) != 1:
        raise ActaDownstreamInvestigationError(
            "Expected exactly one completed Stage-8 publication"
        )

    stage8_root = stage8_roots[0]
    stage_root = stage8_root.parent

    metadata_root = (
        stage_root
        / "downstream_metadata"
    )

    if not (
        metadata_root
        / "_SUCCESS"
    ).is_file():
        raise ActaDownstreamInvestigationError(
            "Finalized downstream metadata is missing"
        )

    stage8_success = _read_json(
        stage8_root
        / "_SUCCESS"
    )
    stage8_summary = _read_json(
        stage8_root
        / "global_summary.json"
    )
    metadata_success = _read_json(
        metadata_root
        / "_SUCCESS"
    )
    metadata_summary = _read_json(
        metadata_root
        / "global_summary.json"
    )

    _validate_upstream(
        canonical_config_hash=(
            canonical_loaded.sha256
        ),
        stage8_success=stage8_success,
        stage8_summary=stage8_summary,
        metadata_success=metadata_success,
        metadata_summary=metadata_summary,
    )

    snapshot = str(
        stage8_success[
            "snapshot"
        ]
    )

    threshold_mA = int(
        stage8_summary[
            "near_duplicate_threshold_mA"
        ]
    )

    if (
        int(
            stage8_summary[
                "total_near_duplicate_count"
            ]
        )
        != expected[
            "input_geometric_near_duplicate_pairs"
        ]
    ):
        raise ActaDownstreamInvestigationError(
            "Stage-8 population does not match the declared "
            "preflight oracle"
        )

    metadata_path = (
        metadata_root
        / metadata_success[
            "entry_metadata"
        ]
    )

    metadata = _metadata_map(
        metadata_path,
        snapshot=snapshot,
        resolution_config=downstream[
            "resolution_filter"
        ],
    )

    if len(metadata) != int(
        metadata_summary[
            "participating_deposition_count"
        ]
    ):
        raise ActaDownstreamInvestigationError(
            "Metadata population accounting mismatch"
        )

    output_root = (
        stage_root
        / _output_directory_name(
            downstream
        )
    )

    finalized = (
        output_root
        / "finalized"
    )

    finalized.mkdir(
        parents=True,
        exist_ok=True,
    )

    success_path = (
        output_root
        / "_SUCCESS"
    )

    if success_path.exists():
        success_path.unlink()

    output_path = (
        finalized
        / "detailed_inspection_candidates.parquet"
    )

    temporary = output_path.with_suffix(
        output_path.suffix + ".tmp"
    )

    if temporary.exists():
        temporary.unlink()

    counters = {
        "input": 0,
        "resolution_rejected": 0,
        "after_resolution": 0,
        "pandda_rejected": 0,
        "after_pandda": 0,
        "same_deposition_rejected": 0,
        "after_same_deposition": 0,
    }

    final_pdb_ids: set[str] = set()

    writer = pq.ParquetWriter(
        temporary,
        OUTPUT_SCHEMA,
        compression="zstd",
        version="2.6",
        use_dictionary=True,
    )

    try:
        _process_pair_file(
            stage8_root
            / stage8_success[
                "candidate_near_duplicates"
            ],
            source_name=(
                "candidate_near_duplicates"
            ),
            metadata=metadata,
            writer=writer,
            snapshot=snapshot,
            threshold_mA=threshold_mA,
            counters=counters,
            final_pdb_ids=final_pdb_ids,
        )

        _process_pair_file(
            stage8_root
            / stage8_success[
                "m1_near_duplicates"
            ],
            source_name=(
                "m1_near_duplicates"
            ),
            metadata=metadata,
            writer=writer,
            snapshot=snapshot,
            threshold_mA=threshold_mA,
            counters=counters,
            final_pdb_ids=final_pdb_ids,
        )
    except Exception:
        writer.close()

        if temporary.exists():
            temporary.unlink()

        raise

    writer.close()

    # Exact sequential accounting.
    if (
        counters["input"]
        != counters[
            "resolution_rejected"
        ]
        + counters[
            "after_resolution"
        ]
    ):
        raise ActaDownstreamInvestigationError(
            "Resolution-filter accounting failed"
        )

    if (
        counters[
            "after_resolution"
        ]
        != counters[
            "pandda_rejected"
        ]
        + counters[
            "after_pandda"
        ]
    ):
        raise ActaDownstreamInvestigationError(
            "PanDDA-filter accounting failed"
        )

    if (
        counters[
            "after_pandda"
        ]
        != counters[
            "same_deposition_rejected"
        ]
        + counters[
            "after_same_deposition"
        ]
    ):
        raise ActaDownstreamInvestigationError(
            "Same-deposition accounting failed"
        )

    observed = {
        "input_geometric_near_duplicate_pairs": (
            counters["input"]
        ),
        "after_resolution_filter": (
            counters[
                "after_resolution"
            ]
        ),
        "after_pandda_group_deposition_filter": (
            counters[
                "after_pandda"
            ]
        ),
        "after_same_deposition_filter": (
            counters[
                "after_same_deposition"
            ]
        ),
    }

    if observed != expected:
        if temporary.exists():
            temporary.unlink()

        raise ActaDownstreamInvestigationError(
            "Production counts differ from the independent "
            f"preflight oracle: expected={expected!r}, "
            f"observed={observed!r}"
        )

    output_rows = pq.read_metadata(
        temporary
    ).num_rows

    if (
        output_rows
        != counters[
            "after_same_deposition"
        ]
    ):
        temporary.unlink()

        raise ActaDownstreamInvestigationError(
            "Final candidate Parquet row-count mismatch"
        )

    temporary.replace(
        output_path
    )

    provenance = {
        "snapshot": snapshot,
        "cleaning_protocol": (
            stage8_success[
                "cleaning_protocol"
            ]
        ),
        "canonical_config_sha256": (
            canonical_loaded.sha256
        ),
        "downstream_config_sha256": (
            downstream_hash
        ),
        "full_bri_nn_pipeline_git_commit": (
            stage8_success[
                "full_bri_nn_pipeline_git_commit"
            ]
        ),
        "full_bri_nn_metadata_correction_git_commit": (
            stage8_success.get(
                "metadata_correction_git_commit"
            )
        ),
        "metadata_task_pipeline_git_commits": (
            metadata_success[
                "metadata_task_pipeline_git_commits"
            ]
        ),
        "metadata_finalizer_git_commit": (
            metadata_success[
                "metadata_finalizer_git_commit"
            ]
        ),
        "acta_downstream_pipeline_git_commit": (
            pipeline_commit
        ),
    }

    resolution_config = downstream[
        "resolution_filter"
    ]

    summary = {
        "summary_schema_name": (
            SUMMARY_SCHEMA_NAME
        ),
        "summary_schema_version": (
            SUMMARY_SCHEMA_VERSION
        ),
        **provenance,
        "scientific_filtering_performed": True,
        "scientific_sequence": (
            "resolution_then_pandda_then_same_deposition"
        ),
        "geometric_input_definition": (
            "complete_BRI_L_infinity_less_than_or_equal_threshold"
        ),
        "near_duplicate_threshold_mA": (
            threshold_mA
        ),
        "near_duplicate_threshold_angstrom": (
            float(
                stage8_summary[
                    "near_duplicate_threshold_angstrom"
                ]
            )
        ),
        "resolution_maximum_angstrom": (
            float(
                resolution_config[
                    "maximum_angstrom"
                ]
            )
        ),
        "resolution_operator": (
            resolution_config[
                "operator"
            ]
        ),
        "missing_resolution_policy": (
            resolution_config[
                "missing_resolution_policy"
            ]
        ),
        "input_geometric_near_duplicate_pair_count": (
            counters["input"]
        ),
        "resolution_rejected_pair_count": (
            counters[
                "resolution_rejected"
            ]
        ),
        "after_resolution_pair_count": (
            counters[
                "after_resolution"
            ]
        ),
        "pandda_group_rejected_pair_count": (
            counters[
                "pandda_rejected"
            ]
        ),
        "after_pandda_pair_count": (
            counters[
                "after_pandda"
            ]
        ),
        "same_deposition_rejected_pair_count": (
            counters[
                "same_deposition_rejected"
            ]
        ),
        "detailed_inspection_candidate_pair_count": (
            counters[
                "after_same_deposition"
            ]
        ),
        "detailed_inspection_participating_deposition_count": (
            len(final_pdb_ids)
        ),
        "automatic_virus_filtering": False,
        "automatic_ribosome_filtering": False,
        "manual_evaluation_performed": False,
        "preflight_oracle_counts_equal": True,
        "pair_accounting_valid": True,
        "processing_error_count": 0,
    }

    _write_json_atomic(
        summary,
        output_root
        / "global_summary.json",
    )

    success = {
        "success_schema_name": (
            SUCCESS_SCHEMA_NAME
        ),
        "success_schema_version": (
            SUCCESS_SCHEMA_VERSION
        ),
        **provenance,
        "detailed_inspection_candidates": (
            "finalized/"
            "detailed_inspection_candidates.parquet"
        ),
        "global_summary": (
            "global_summary.json"
        ),
    }

    # _SUCCESS is deliberately written last.
    _write_json_atomic(
        success,
        success_path,
    )

    print(
        "===== ACTA DOWNSTREAM INVESTIGATION ====="
    )
    print(
        "Input geometric near-duplicates:",
        f"{counters['input']:,}",
    )
    print(
        "After resolution <= "
        f"{float(resolution_config['maximum_angstrom']):g} A:",
        f"{counters['after_resolution']:,}",
    )
    print(
        "Resolution rejected:",
        f"{counters['resolution_rejected']:,}",
    )
    print(
        "After PanDDA Group-deposition rejection:",
        f"{counters['after_pandda']:,}",
    )
    print(
        "PanDDA additionally rejected:",
        f"{counters['pandda_rejected']:,}",
    )
    print(
        "After same-deposition removal:",
        f"{counters['after_same_deposition']:,}",
    )
    print(
        "Same-deposition additionally rejected:",
        f"{counters['same_deposition_rejected']:,}",
    )
    print(
        "Participating PDB depositions:",
        f"{len(final_pdb_ids):,}",
    )
    print(
        "Manual virus/ribosome filtering: NO"
    )
    print(
        "Preflight oracle equality: PASS"
    )
    print(
        "Pair accounting: PASS"
    )
    print(
        "ACTA DOWNSTREAM INVESTIGATION PUBLICATION: PASS"
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
