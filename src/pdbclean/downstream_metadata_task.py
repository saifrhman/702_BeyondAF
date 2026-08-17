"""Distributed producer for snapshot-consistent downstream metadata.

Each logical task processes one deterministic partition of the PDB
depositions participating in the frozen Stage-8 near-duplicate output.

This stage extracts metadata only.  It performs no downstream
scientific filtering.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import json
from pathlib import Path
import time

import pyarrow as pa
import pyarrow.parquet as pq

from pdbclean.config import load_config
from pdbclean.downstream_metadata import parse_entry_metadata_bytes
from pdbclean.downstream_metadata_production import (
    ENTRY_METADATA_SCHEMA,
    _participating_pdb_ids,
    _read_json,
    _validate_commit,
)
from pdbclean.snapshot import (
    SnapshotTransportError,
    download_verified_s3_object_bytes,
)


class DownstreamMetadataTaskError(RuntimeError):
    """Raised when a distributed metadata task cannot complete safely."""


TASK_SUMMARY_SCHEMA_NAME = "pdbclean_downstream_metadata_task_summary"
TASK_SUMMARY_SCHEMA_VERSION = "1.0"


def _write_json_atomic(data: dict, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)

    temporary = path.with_suffix(path.suffix + ".tmp")

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


def _task_count(
    row_count: int,
    batch_size: int,
) -> int:
    if row_count <= 0:
        raise DownstreamMetadataTaskError(
            "Participating deposition population must be non-empty"
        )

    if batch_size <= 0:
        raise DownstreamMetadataTaskError(
            "Metadata batch size must be positive"
        )

    return (
        row_count
        + batch_size
        - 1
    ) // batch_size


def _resolve_inputs(
    *,
    config: dict,
    repo: Path,
) -> tuple[
    Path,
    Path,
    dict,
    str,
    list[dict],
]:
    protocol = config["release"]["protocol_version"]

    storage_root = Path(
        config["storage"]["output_root"]
    )

    if not storage_root.is_absolute():
        storage_root = repo / storage_root

    completed_nn = sorted(
        path.parent
        for path in storage_root.glob(
            f"*/{protocol}/full_bri_nn/_SUCCESS"
        )
    )

    if len(completed_nn) != 1:
        raise DownstreamMetadataTaskError(
            "Expected exactly one completed Stage-8 "
            "paper-faithful NN publication"
        )

    stage8_root = completed_nn[0]
    stage_root = stage8_root.parent

    stage8_success = _read_json(
        stage8_root / "_SUCCESS"
    )
    stage8_summary = _read_json(
        stage8_root / "global_summary.json"
    )

    if (
        stage8_success.get("success_schema_name")
        != "pdbclean_stage8_full_bri_nn_success"
        or stage8_success.get("success_schema_version")
        != "1.0"
    ):
        raise DownstreamMetadataTaskError(
            "Unexpected Stage-8 success schema"
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
        raise DownstreamMetadataTaskError(
            "Stage-8 publication is not fully validated"
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

    manifest_path = (
        stage_root.parent
        / "bronze/source_manifest.parquet"
    )

    if not manifest_path.is_file():
        raise DownstreamMetadataTaskError(
            "Frozen bronze source manifest is missing"
        )

    manifest_rows = pq.read_table(
        manifest_path,
        columns=[
            "snapshot",
            "pdb_id",
            "s3_key",
            "size_bytes",
            "etag",
        ],
    ).to_pylist()

    selected: dict[str, dict] = {}

    for row in manifest_rows:
        pdb_id = str(
            row["pdb_id"]
        ).lower()

        if pdb_id not in participating:
            continue

        if pdb_id in selected:
            raise DownstreamMetadataTaskError(
                f"Duplicate manifest entry for {pdb_id}"
            )

        if str(row["snapshot"]) != snapshot:
            raise DownstreamMetadataTaskError(
                f"Snapshot mismatch for {pdb_id}"
            )

        selected[pdb_id] = {
            "snapshot": snapshot,
            "pdb_id": pdb_id,
            "s3_key": str(
                row["s3_key"]
            ),
            "size_bytes": int(
                row["size_bytes"]
            ),
            "etag": str(
                row["etag"]
            ),
        }

    if set(selected) != participating:
        missing = participating - set(selected)

        raise DownstreamMetadataTaskError(
            "Bronze manifest does not exactly cover "
            "Stage-8 participating depositions; "
            f"missing={len(missing)}"
        )

    # Cross-check source object identity against the canonical BRI.
    source_by_pdb: dict[
        str,
        tuple[str, str],
    ] = {}

    bri_path = (
        stage_root
        / "bri/finalized/bri.parquet"
    )

    pf = pq.ParquetFile(bri_path)

    for batch in pf.iter_batches(
        columns=[
            "pdb_id",
            "source_mmcif_key",
            "source_etag",
        ],
        batch_size=65_536,
    ):
        values = batch.to_pydict()

        for pdb_id, key, etag in zip(
            values["pdb_id"],
            values["source_mmcif_key"],
            values["source_etag"],
            strict=True,
        ):
            pdb_id = str(pdb_id).lower()

            if pdb_id not in participating:
                continue

            identity = (
                str(key),
                str(etag),
            )

            previous = source_by_pdb.get(
                pdb_id
            )

            if (
                previous is not None
                and previous != identity
            ):
                raise DownstreamMetadataTaskError(
                    "Conflicting canonical BRI source "
                    f"identity for {pdb_id}"
                )

            source_by_pdb[pdb_id] = identity

    if set(source_by_pdb) != participating:
        raise DownstreamMetadataTaskError(
            "Canonical BRI provenance does not cover "
            "the participating population"
        )

    for pdb_id, row in selected.items():
        expected = (
            row["s3_key"],
            row["etag"],
        )

        if source_by_pdb[pdb_id] != expected:
            raise DownstreamMetadataTaskError(
                "Bronze/BRI source identity mismatch "
                f"for {pdb_id}"
            )

    ordered = [
        selected[pdb_id]
        for pdb_id in sorted(selected)
    ]

    return (
        stage_root,
        stage8_root,
        stage8_success,
        snapshot,
        ordered,
    )


def main() -> int:
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--config",
        required=True,
        type=Path,
    )
    parser.add_argument(
        "--task-id",
        required=True,
        type=int,
    )
    parser.add_argument(
        "--task-count",
        required=True,
        type=int,
    )
    parser.add_argument(
        "--pipeline-git-commit",
        required=True,
    )

    args = parser.parse_args()

    producer_commit = _validate_commit(
        args.pipeline_git_commit
    )

    loaded = load_config(
        args.config
    )
    config = loaded.data

    repo = Path.cwd()

    (
        stage_root,
        stage8_root,
        stage8_success,
        snapshot,
        ordered_manifest,
    ) = _resolve_inputs(
        config=config,
        repo=repo,
    )

    batch_size = int(
        config["execution"]["batch_size"]
    )

    expected_task_count = _task_count(
        len(ordered_manifest),
        batch_size,
    )

    if args.task_count != expected_task_count:
        raise DownstreamMetadataTaskError(
            "Task-count mismatch: "
            f"expected {expected_task_count}, "
            f"received {args.task_count}"
        )

    if (
        args.task_id < 0
        or args.task_id >= expected_task_count
    ):
        raise DownstreamMetadataTaskError(
            f"Task ID {args.task_id} outside "
            f"0..{expected_task_count - 1}"
        )

    start = (
        args.task_id
        * batch_size
    )
    stop = min(
        start + batch_size,
        len(ordered_manifest),
    )

    partition = ordered_manifest[
        start:stop
    ]

    if not partition:
        raise DownstreamMetadataTaskError(
            "Resolved task partition is empty"
        )

    bucket_url = config[
        "snapshot"
    ]["bucket_url"]

    timeout_seconds = int(
        config["execution"][
            "connection_timeout_seconds"
        ]
    )
    max_retries = int(
        config["execution"][
            "max_retries"
        ]
    )
    download_concurrency = int(
        config["execution"][
            "download_concurrency"
        ]
    )

    if download_concurrency <= 0:
        raise DownstreamMetadataTaskError(
            "download_concurrency must be positive"
        )

    def extract_one(
        source: dict,
    ) -> dict:
        compressed_bytes = None

        for attempt in range(
            max_retries + 1
        ):
            try:
                compressed_bytes = (
                    download_verified_s3_object_bytes(
                        bucket_url=bucket_url,
                        s3_key=source["s3_key"],
                        expected_size_bytes=(
                            source["size_bytes"]
                        ),
                        expected_etag=source["etag"],
                        timeout_seconds=(
                            timeout_seconds
                        ),
                    )
                )
                break

            except SnapshotTransportError:
                if attempt >= max_retries:
                    raise

                time.sleep(
                    min(
                        30,
                        2 ** attempt,
                    )
                )

        if compressed_bytes is None:
            raise DownstreamMetadataTaskError(
                "Verified download unexpectedly absent "
                f"for {source['pdb_id']}"
            )

        metadata = parse_entry_metadata_bytes(
            compressed_bytes,
            pdb_id=source["pdb_id"],
        )

        return {
            "snapshot": snapshot,
            "pdb_id": source["pdb_id"],
            "s3_key": source["s3_key"],
            "size_bytes": source[
                "size_bytes"
            ],
            "etag": source["etag"],
            "experimental_methods": list(
                metadata.experimental_methods
            ),
            "refine_ls_d_res_high": list(
                metadata.refine_ls_d_res_high
            ),
            "em_3d_reconstruction_resolution": list(
                metadata.em_3d_reconstruction_resolution
            ),
            "initial_deposition_date": (
                metadata.initial_deposition_date
            ),
            "struct_title": (
                metadata.struct_title
            ),
            "struct_keywords_text": list(
                metadata.struct_keywords_text
            ),
            "deposit_group_ids": list(
                metadata.deposit_group_ids
            ),
            "deposit_group_titles": list(
                metadata.deposit_group_titles
            ),
            "deposit_group_descriptions": list(
                metadata.deposit_group_descriptions
            ),
            "deposit_group_types": list(
                metadata.deposit_group_types
            ),
            "has_deposit_group": (
                metadata.has_deposit_group
            ),
            "deposit_group_mentions_pandda": (
                metadata.deposit_group_mentions_pandda
            ),
            "entry_mentions_pandda": (
                metadata.entry_mentions_pandda
            ),
        }

    print(
        "===== DOWNSTREAM METADATA TASK ====="
    )
    print(
        "Task:",
        f"{args.task_id}/{expected_task_count - 1}",
    )
    print(
        "Input depositions:",
        f"{len(partition):,}",
    )
    print(
        "PDB range:",
        partition[0]["pdb_id"],
        "through",
        partition[-1]["pdb_id"],
    )
    print(
        "Download concurrency:",
        download_concurrency,
    )

    extracted: list[dict] = []

    with ThreadPoolExecutor(
        max_workers=download_concurrency
    ) as executor:
        futures = {
            executor.submit(
                extract_one,
                source,
            ): source["pdb_id"]
            for source in partition
        }

        completed = 0

        for future in as_completed(
            futures
        ):
            pdb_id = futures[future]

            try:
                result = future.result()
            except Exception as exc:
                raise DownstreamMetadataTaskError(
                    f"Task {args.task_id} failed on "
                    f"{pdb_id}: {type(exc).__name__}: {exc}"
                ) from exc

            extracted.append(result)
            completed += 1

            if (
                completed % 100 == 0
                or completed == len(partition)
            ):
                print(
                    "Verified/extracted:",
                    f"{completed:,}/{len(partition):,}",
                    flush=True,
                )

    extracted.sort(
        key=lambda row: row["pdb_id"]
    )

    expected_ids = [
        row["pdb_id"]
        for row in partition
    ]
    observed_ids = [
        row["pdb_id"]
        for row in extracted
    ]

    if observed_ids != expected_ids:
        raise DownstreamMetadataTaskError(
            "Task output identities differ from "
            "deterministic input partition"
        )

    output_root = (
        stage_root
        / "downstream_metadata"
    )

    shard_dir = (
        output_root
        / "tasks/metadata"
    )
    summary_dir = (
        output_root
        / "tasks/summaries"
    )

    shard_dir.mkdir(
        parents=True,
        exist_ok=True,
    )
    summary_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    shard_path = (
        shard_dir
        / f"task_{args.task_id}.parquet"
    )
    temporary = shard_path.with_suffix(
        shard_path.suffix + ".tmp"
    )

    if temporary.exists():
        temporary.unlink()

    table = pa.Table.from_pylist(
        extracted,
        schema=ENTRY_METADATA_SCHEMA,
    )

    pq.write_table(
        table,
        temporary,
        compression="zstd",
        version="2.6",
        use_dictionary=True,
    )

    if (
        pq.read_metadata(
            temporary
        ).num_rows
        != len(partition)
    ):
        temporary.unlink()

        raise DownstreamMetadataTaskError(
            "Task Parquet row-count validation failed"
        )

    temporary.replace(shard_path)

    verified_bytes = sum(
        row["size_bytes"]
        for row in extracted
    )

    summary = {
        "summary_schema_name": (
            TASK_SUMMARY_SCHEMA_NAME
        ),
        "summary_schema_version": (
            TASK_SUMMARY_SCHEMA_VERSION
        ),
        "task_id": args.task_id,
        "task_count": expected_task_count,
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
        "input_deposition_count": len(
            partition
        ),
        "successful_deposition_count": len(
            extracted
        ),
        "verified_source_bytes": (
            verified_bytes
        ),
        "first_pdb_id": (
            extracted[0]["pdb_id"]
        ),
        "last_pdb_id": (
            extracted[-1]["pdb_id"]
        ),
        "processing_error_count": 0,
        "scientific_filtering_performed": False,
        "metadata_shard": (
            f"tasks/metadata/task_{args.task_id}.parquet"
        ),
    }

    _write_json_atomic(
        summary,
        summary_dir
        / f"task_{args.task_id}.json",
    )

    print()
    print(
        "Task output:",
        shard_path,
    )
    print(
        "Verified bytes:",
        f"{verified_bytes:,}",
    )
    print(
        "Scientific filtering performed: NO"
    )
    print(
        "DOWNSTREAM METADATA TASK: PASS"
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
