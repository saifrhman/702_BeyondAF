"""Validate and finalize distributed downstream metadata extraction."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

from pdbclean.config import load_config
from pdbclean.downstream_metadata_production import (
    ENTRY_METADATA_SCHEMA,
    _participating_pdb_ids,
    _read_json,
    _validate_commit,
)
from pdbclean.downstream_metadata_task import (
    TASK_SUMMARY_SCHEMA_NAME,
    TASK_SUMMARY_SCHEMA_VERSION,
    _task_count,
)


class DownstreamMetadataFinalizeError(RuntimeError):
    """Raised when metadata finalization cannot proceed safely."""


SUCCESS_SCHEMA_NAME = "pdbclean_downstream_metadata_success"
SUCCESS_SCHEMA_VERSION = "1.1"

GLOBAL_SUMMARY_SCHEMA_NAME = (
    "pdbclean_downstream_metadata_global_summary"
)
GLOBAL_SUMMARY_SCHEMA_VERSION = "1.1"


def _write_json_atomic(
    data: dict,
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


def main() -> int:
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--config",
        required=True,
        type=Path,
    )
    parser.add_argument(
        "--producer-git-commit",
        required=True,
    )
    parser.add_argument(
        "--finalizer-git-commit",
        required=True,
    )

    args = parser.parse_args()

    producer_commit = _validate_commit(
        args.producer_git_commit
    )
    finalizer_commit = _validate_commit(
        args.finalizer_git_commit
    )

    loaded = load_config(
        args.config
    )
    config = loaded.data

    repo = Path.cwd()

    protocol = config[
        "release"
    ]["protocol_version"]

    storage_root = Path(
        config["storage"]["output_root"]
    )

    if not storage_root.is_absolute():
        storage_root = (
            repo / storage_root
        )

    completed_nn = sorted(
        path.parent
        for path in storage_root.glob(
            f"*/{protocol}/full_bri_nn/_SUCCESS"
        )
    )

    if len(completed_nn) != 1:
        raise DownstreamMetadataFinalizeError(
            "Expected exactly one completed Stage-8 publication"
        )

    stage8_root = completed_nn[0]
    stage_root = stage8_root.parent

    stage8_success = _read_json(
        stage8_root / "_SUCCESS"
    )
    stage8_summary = _read_json(
        stage8_root / "global_summary.json"
    )

    if not (
        stage8_summary.get(
            "exhaustive_oracle_pair_set_equal"
        )
        and stage8_summary.get(
            "exhaustive_oracle_distances_equal"
        )
        and stage8_summary.get(
            "processing_error_count"
        ) == 0
    ):
        raise DownstreamMetadataFinalizeError(
            "Stage-8 validation gate failed"
        )

    snapshot = str(
        stage8_success["snapshot"]
    )

    participating = _participating_pdb_ids(
        [
            stage8_root
            / stage8_success[
                "candidate_near_duplicates"
            ],
            stage8_root
            / stage8_success[
                "m1_near_duplicates"
            ],
        ]
    )

    ordered_ids = sorted(
        participating
    )

    batch_size = int(
        config["execution"]["batch_size"]
    )

    expected_task_count = _task_count(
        len(ordered_ids),
        batch_size,
    )

    output_root = (
        stage_root
        / "downstream_metadata"
    )

    # Once finalization starts, no old success marker may survive.
    success_path = (
        output_root
        / "_SUCCESS"
    )

    if success_path.exists():
        success_path.unlink()

    tables: list[pa.Table] = []
    total_verified_bytes = 0

    for task_id in range(
        expected_task_count
    ):
        shard_path = (
            output_root
            / "tasks/metadata"
            / f"task_{task_id}.parquet"
        )
        summary_path = (
            output_root
            / "tasks/summaries"
            / f"task_{task_id}.json"
        )

        if not shard_path.is_file():
            raise DownstreamMetadataFinalizeError(
                f"Missing metadata shard for task {task_id}"
            )

        if not summary_path.is_file():
            raise DownstreamMetadataFinalizeError(
                f"Missing summary for task {task_id}"
            )

        summary = _read_json(
            summary_path
        )

        if (
            summary.get(
                "summary_schema_name"
            )
            != TASK_SUMMARY_SCHEMA_NAME
            or summary.get(
                "summary_schema_version"
            )
            != TASK_SUMMARY_SCHEMA_VERSION
        ):
            raise DownstreamMetadataFinalizeError(
                f"Unexpected summary schema for task {task_id}"
            )

        checks = {
            "task_id": task_id,
            "task_count": (
                expected_task_count
            ),
            "batch_size": batch_size,
            "snapshot": snapshot,
            "cleaning_protocol": (
                stage8_success[
                    "cleaning_protocol"
                ]
            ),
            "full_bri_nn_pipeline_git_commit": (
                stage8_success[
                    "full_bri_nn_pipeline_git_commit"
                ]
            ),
            "metadata_task_pipeline_git_commit": (
                producer_commit
            ),
            "config_sha256": loaded.sha256,
            "processing_error_count": 0,
            "scientific_filtering_performed": False,
        }

        for field, expected in checks.items():
            if summary.get(field) != expected:
                raise DownstreamMetadataFinalizeError(
                    f"Task {task_id} summary mismatch "
                    f"for {field}: expected {expected!r}, "
                    f"found {summary.get(field)!r}"
                )

        start = (
            task_id
            * batch_size
        )
        stop = min(
            start + batch_size,
            len(ordered_ids),
        )

        expected_ids = ordered_ids[
            start:stop
        ]

        if (
            summary.get(
                "input_deposition_count"
            )
            != len(expected_ids)
            or summary.get(
                "successful_deposition_count"
            )
            != len(expected_ids)
        ):
            raise DownstreamMetadataFinalizeError(
                f"Task {task_id} source accounting mismatch"
            )

        table = pq.read_table(
            shard_path
        )

        if table.num_rows != len(
            expected_ids
        ):
            raise DownstreamMetadataFinalizeError(
                f"Task {task_id} Parquet row-count mismatch"
            )

        observed_ids = [
            str(value).lower()
            for value in table[
                "pdb_id"
            ].to_pylist()
        ]

        if observed_ids != expected_ids:
            raise DownstreamMetadataFinalizeError(
                f"Task {task_id} PDB identity/order mismatch"
            )

        snapshots = set(
            table[
                "snapshot"
            ].to_pylist()
        )

        if snapshots != {snapshot}:
            raise DownstreamMetadataFinalizeError(
                f"Task {task_id} snapshot mismatch"
            )

        total_verified_bytes += int(
            summary[
                "verified_source_bytes"
            ]
        )

        tables.append(
            table.cast(
                ENTRY_METADATA_SCHEMA
            )
        )

    merged = pa.concat_tables(
        tables,
        promote_options="none",
    )

    if merged.num_rows != len(
        ordered_ids
    ):
        raise DownstreamMetadataFinalizeError(
            "Merged metadata row-count mismatch"
        )

    merged_ids = [
        str(value).lower()
        for value in merged[
            "pdb_id"
        ].to_pylist()
    ]

    if merged_ids != ordered_ids:
        raise DownstreamMetadataFinalizeError(
            "Merged metadata identities/order mismatch"
        )

    if len(
        set(merged_ids)
    ) != len(
        merged_ids
    ):
        raise DownstreamMetadataFinalizeError(
            "Merged metadata contains duplicate PDB IDs"
        )

    # Re-check exact immutable source identities against bronze.
    manifest_path = (
        stage_root.parent
        / "bronze/source_manifest.parquet"
    )

    manifest = pq.read_table(
        manifest_path,
        columns=[
            "pdb_id",
            "s3_key",
            "size_bytes",
            "etag",
        ],
    ).to_pylist()

    expected_source = {}

    for row in manifest:
        pdb_id = str(
            row["pdb_id"]
        ).lower()

        if pdb_id not in participating:
            continue

        expected_source[pdb_id] = (
            str(row["s3_key"]),
            int(row["size_bytes"]),
            str(row["etag"]),
        )

    if set(
        expected_source
    ) != participating:
        raise DownstreamMetadataFinalizeError(
            "Bronze source identity population mismatch"
        )

    for row in merged.to_pylist():
        pdb_id = row[
            "pdb_id"
        ]

        observed = (
            row["s3_key"],
            int(row["size_bytes"]),
            row["etag"],
        )

        if (
            expected_source[pdb_id]
            != observed
        ):
            raise DownstreamMetadataFinalizeError(
                f"Final source identity mismatch for {pdb_id}"
            )

    finalized = (
        output_root
        / "finalized"
    )

    finalized.mkdir(
        parents=True,
        exist_ok=True,
    )

    output_path = (
        finalized
        / "entry_metadata.parquet"
    )

    temporary = output_path.with_suffix(
        output_path.suffix + ".tmp"
    )

    if temporary.exists():
        temporary.unlink()

    pq.write_table(
        merged,
        temporary,
        compression="zstd",
        version="2.6",
        use_dictionary=True,
    )

    if (
        pq.read_metadata(
            temporary
        ).num_rows
        != len(ordered_ids)
    ):
        temporary.unlink()

        raise DownstreamMetadataFinalizeError(
            "Finalized Parquet validation failed"
        )

    temporary.replace(
        output_path
    )

    rows = merged.to_pylist()

    methods_present = sum(
        bool(row["experimental_methods"])
        for row in rows
    )
    refine_present = sum(
        bool(row["refine_ls_d_res_high"])
        for row in rows
    )
    em_present = sum(
        bool(
            row[
                "em_3d_reconstruction_resolution"
            ]
        )
        for row in rows
    )
    group_present = sum(
        bool(row["has_deposit_group"])
        for row in rows
    )
    group_pandda = sum(
        bool(
            row[
                "deposit_group_mentions_pandda"
            ]
        )
        for row in rows
    )
    entry_pandda = sum(
        bool(
            row[
                "entry_mentions_pandda"
            ]
        )
        for row in rows
    )

    provenance = {
        "snapshot": snapshot,
        "cleaning_protocol": (
            stage8_success[
                "cleaning_protocol"
            ]
        ),
        "full_bri_nn_pipeline_git_commit": (
            stage8_success[
                "full_bri_nn_pipeline_git_commit"
            ]
        ),
        "metadata_task_pipeline_git_commit": (
            producer_commit
        ),
        "metadata_finalizer_git_commit": (
            finalizer_commit
        ),
        "config_sha256": loaded.sha256,
    }

    summary = {
        "summary_schema_name": (
            GLOBAL_SUMMARY_SCHEMA_NAME
        ),
        "summary_schema_version": (
            GLOBAL_SUMMARY_SCHEMA_VERSION
        ),
        **provenance,
        "scientific_filtering_performed": False,
        "participating_deposition_count": (
            len(ordered_ids)
        ),
        "batch_size": batch_size,
        "task_count": expected_task_count,
        "verified_source_object_count": (
            len(ordered_ids)
        ),
        "verified_source_bytes": (
            total_verified_bytes
        ),
        "experimental_method_present_count": (
            methods_present
        ),
        "refine_resolution_present_count": (
            refine_present
        ),
        "em_resolution_present_count": (
            em_present
        ),
        "deposit_group_present_count": (
            group_present
        ),
        "deposit_group_mentions_pandda_count": (
            group_pandda
        ),
        "entry_mentions_pandda_count": (
            entry_pandda
        ),
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
        "batch_size": batch_size,
        "task_count": expected_task_count,
        "entry_metadata": (
            "finalized/entry_metadata.parquet"
        ),
        "global_summary": (
            "global_summary.json"
        ),
    }

    # Completion marker strictly last.
    _write_json_atomic(
        success,
        success_path,
    )

    print(
        "===== DOWNSTREAM METADATA FINALIZATION ====="
    )
    print(
        "Depositions:",
        f"{len(ordered_ids):,}",
    )
    print(
        "Logical tasks:",
        expected_task_count,
    )
    print(
        "Verified compressed GiB:",
        f"{total_verified_bytes / 1024**3:.3f}",
    )
    print(
        "Experimental method present:",
        f"{methods_present:,}",
    )
    print(
        "Refine resolution present:",
        f"{refine_present:,}",
    )
    print(
        "EM resolution present:",
        f"{em_present:,}",
    )
    print(
        "Deposit-group metadata present:",
        f"{group_present:,}",
    )
    print(
        "Deposit group mentions PanDDA:",
        f"{group_pandda:,}",
    )
    print(
        "Entry metadata mentions PanDDA:",
        f"{entry_pandda:,}",
    )
    print(
        "Scientific filtering performed: NO"
    )
    print(
        "DOWNSTREAM METADATA FINALIZATION: PASS"
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
