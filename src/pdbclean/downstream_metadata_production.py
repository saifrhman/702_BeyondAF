"""Extract snapshot-consistent metadata for Stage-8 near duplicates.

No downstream scientific filters are applied by this stage.
"""

from __future__ import annotations

import argparse
from concurrent.futures import (
    ThreadPoolExecutor,
    as_completed,
)
import json
from pathlib import Path
import time

import pyarrow as pa
import pyarrow.parquet as pq

from pdbclean.config import load_config
from pdbclean.downstream_metadata import (
    parse_entry_metadata_bytes,
)
from pdbclean.snapshot import (
    SnapshotTransportError,
    download_verified_s3_object_bytes,
)


class DownstreamMetadataProductionError(RuntimeError):
    """Raised when metadata extraction cannot be published safely."""


ENTRY_METADATA_SCHEMA = pa.schema(
    [
        pa.field("snapshot", pa.string(), nullable=False),
        pa.field("pdb_id", pa.string(), nullable=False),
        pa.field("s3_key", pa.string(), nullable=False),
        pa.field("size_bytes", pa.int64(), nullable=False),
        pa.field("etag", pa.string(), nullable=False),

        pa.field(
            "experimental_methods",
            pa.list_(pa.string()),
            nullable=False,
        ),
        pa.field(
            "refine_ls_d_res_high",
            pa.list_(pa.float64()),
            nullable=False,
        ),
        pa.field(
            "em_3d_reconstruction_resolution",
            pa.list_(pa.float64()),
            nullable=False,
        ),

        pa.field(
            "initial_deposition_date",
            pa.string(),
            nullable=True,
        ),
        pa.field(
            "struct_title",
            pa.string(),
            nullable=True,
        ),
        pa.field(
            "struct_keywords_text",
            pa.list_(pa.string()),
            nullable=False,
        ),

        pa.field(
            "deposit_group_ids",
            pa.list_(pa.string()),
            nullable=False,
        ),
        pa.field(
            "deposit_group_titles",
            pa.list_(pa.string()),
            nullable=False,
        ),
        pa.field(
            "deposit_group_descriptions",
            pa.list_(pa.string()),
            nullable=False,
        ),
        pa.field(
            "deposit_group_types",
            pa.list_(pa.string()),
            nullable=False,
        ),

        pa.field(
            "has_deposit_group",
            pa.bool_(),
            nullable=False,
        ),
        pa.field(
            "deposit_group_mentions_pandda",
            pa.bool_(),
            nullable=False,
        ),
        pa.field(
            "entry_mentions_pandda",
            pa.bool_(),
            nullable=False,
        ),
    ],
    metadata={
        b"schema_name": (
            b"pdbclean_downstream_entry_metadata"
        ),
        b"schema_version": b"1.0",
        b"scientific_filtering": b"none",
    },
)


def _read_json(path: Path) -> dict:
    if not path.is_file():
        raise DownstreamMetadataProductionError(
            f"Required JSON file missing: {path}"
        )

    return json.loads(
        path.read_text(
            encoding="utf-8"
        )
    )


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


def _validate_commit(value: str) -> str:
    value = value.strip().lower()

    if (
        len(value) != 40
        or any(
            character not in "0123456789abcdef"
            for character in value
        )
    ):
        raise DownstreamMetadataProductionError(
            "Metadata extractor commit must be "
            "a full 40-character Git SHA"
        )

    return value


def _participating_pdb_ids(
    paths: list[Path],
) -> set[str]:
    result: set[str] = set()

    for path in paths:
        if not path.is_file():
            raise DownstreamMetadataProductionError(
                f"Stage-8 result missing: {path}"
            )

        pf = pq.ParquetFile(path)

        for batch in pf.iter_batches(
            columns=[
                "query_pdb_id",
                "subject_pdb_id",
            ],
            batch_size=131_072,
        ):
            values = batch.to_pydict()

            result.update(
                str(value).lower()
                for value
                in values["query_pdb_id"]
            )
            result.update(
                str(value).lower()
                for value
                in values["subject_pdb_id"]
            )

    return result


def main() -> int:
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--config",
        required=True,
        type=Path,
    )
    parser.add_argument(
        "--pipeline-git-commit",
        required=True,
    )

    args = parser.parse_args()

    metadata_commit = _validate_commit(
        args.pipeline_git_commit
    )

    repo = Path.cwd()

    loaded = load_config(
        args.config
    )
    config = loaded.data

    protocol = (
        config["release"]["protocol_version"]
    )

    bucket_url = (
        config["snapshot"]["bucket_url"]
    )

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
        raise DownstreamMetadataProductionError(
            "Expected exactly one completed "
            "paper-faithful Stage-8 publication"
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
        stage8_success.get(
            "success_schema_name"
        )
        != "pdbclean_stage8_full_bri_nn_success"
        or stage8_success.get(
            "success_schema_version"
        )
        != "1.0"
    ):
        raise DownstreamMetadataProductionError(
            "Unexpected Stage-8 success schema"
        )

    if (
        stage8_summary.get(
            "summary_schema_name"
        )
        != "pdbclean_stage8_full_bri_nn_global_summary"
        or stage8_summary.get(
            "summary_schema_version"
        )
        != "1.0"
    ):
        raise DownstreamMetadataProductionError(
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
            "processing_error_count"
        ) == 0
    ):
        raise DownstreamMetadataProductionError(
            "Stage-8 publication is not fully validated"
        )

    snapshot = str(
        stage8_success["snapshot"]
    )

    candidate_path = (
        stage8_root
        / stage8_success[
            "candidate_near_duplicates"
        ]
    )
    m1_path = (
        stage8_root
        / stage8_success[
            "m1_near_duplicates"
        ]
    )

    participating = (
        _participating_pdb_ids(
            [
                candidate_path,
                m1_path,
            ]
        )
    )

    if not participating:
        raise DownstreamMetadataProductionError(
            "Stage-8 near-duplicate population is empty"
        )

    manifest_path = (
        stage_root.parent
        / "bronze/source_manifest.parquet"
    )

    if not manifest_path.is_file():
        raise DownstreamMetadataProductionError(
            "Frozen snapshot source manifest is missing"
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

    manifest_by_pdb: dict[str, dict] = {}

    for row in manifest_rows:
        pdb_id = str(
            row["pdb_id"]
        ).lower()

        if pdb_id not in participating:
            continue

        if pdb_id in manifest_by_pdb:
            raise DownstreamMetadataProductionError(
                f"Duplicate manifest row for {pdb_id}"
            )

        if str(
            row["snapshot"]
        ) != snapshot:
            raise DownstreamMetadataProductionError(
                f"Snapshot mismatch for manifest row {pdb_id}"
            )

        manifest_by_pdb[pdb_id] = row

    if set(
        manifest_by_pdb
    ) != participating:
        missing = (
            participating
            - set(manifest_by_pdb)
        )

        raise DownstreamMetadataProductionError(
            "Manifest does not exactly cover Stage-8 "
            f"participating depositions; missing={len(missing)}"
        )

    # --------------------------------------------------------
    # Cross-check preserved BRI source provenance.
    # --------------------------------------------------------

    bri_path = (
        stage_root
        / "bri/finalized/bri.parquet"
    )

    source_by_pdb: dict[
        str,
        tuple[str, str],
    ] = {}

    bri_pf = pq.ParquetFile(
        bri_path
    )

    for batch in bri_pf.iter_batches(
        columns=[
            "pdb_id",
            "source_mmcif_key",
            "source_etag",
        ],
        batch_size=65_536,
    ):
        data = batch.to_pydict()

        for pdb_id, key, etag in zip(
            data["pdb_id"],
            data["source_mmcif_key"],
            data["source_etag"],
        ):
            pdb_id = str(
                pdb_id
            ).lower()

            if pdb_id not in participating:
                continue

            value = (
                str(key),
                str(etag),
            )

            previous = source_by_pdb.get(
                pdb_id
            )

            if (
                previous is not None
                and previous != value
            ):
                raise DownstreamMetadataProductionError(
                    f"Conflicting canonical BRI source "
                    f"provenance for {pdb_id}"
                )

            source_by_pdb[
                pdb_id
            ] = value

    if set(
        source_by_pdb
    ) != participating:
        raise DownstreamMetadataProductionError(
            "Canonical BRI source provenance does not "
            "cover every participating deposition"
        )

    for pdb_id in sorted(
        participating
    ):
        manifest_row = (
            manifest_by_pdb[pdb_id]
        )

        expected = (
            str(
                manifest_row["s3_key"]
            ),
            str(
                manifest_row["etag"]
            ),
        )

        if source_by_pdb[
            pdb_id
        ] != expected:
            raise DownstreamMetadataProductionError(
                f"Manifest/BRI source provenance mismatch "
                f"for {pdb_id}"
            )

    print(
        "Participating depositions:",
        f"{len(participating):,}",
    )
    print(
        "Manifest + BRI source provenance:",
        "PASS",
    )

    # --------------------------------------------------------
    # Verified download + metadata extraction.
    # --------------------------------------------------------

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

    worker_count = int(
        config["execution"][
            "download_concurrency"
        ]
    )

    if worker_count < 1:
        raise DownstreamMetadataProductionError(
            "download_concurrency must be positive"
        )

    ordered_manifest = [
        manifest_by_pdb[pdb_id]
        for pdb_id in sorted(
            participating
        )
    ]

    def extract_one(
        row: dict,
    ) -> dict:
        pdb_id = str(
            row["pdb_id"]
        ).lower()

        compressed_bytes = None

        for attempt in range(
            max_retries + 1
        ):
            try:
                compressed_bytes = (
                    download_verified_s3_object_bytes(
                        bucket_url=bucket_url,
                        s3_key=str(
                            row["s3_key"]
                        ),
                        expected_size_bytes=int(
                            row["size_bytes"]
                        ),
                        expected_etag=str(
                            row["etag"]
                        ),
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
            raise DownstreamMetadataProductionError(
                f"Download unexpectedly absent for {pdb_id}"
            )

        metadata = (
            parse_entry_metadata_bytes(
                compressed_bytes,
                pdb_id=pdb_id,
            )
        )

        return {
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

    extracted: list[dict] = []

    print(
        "Downloading and verifying exact snapshot mmCIFs..."
    )
    print(
        "Download concurrency:",
        worker_count,
    )

    with ThreadPoolExecutor(
        max_workers=worker_count
    ) as executor:
        futures = {
            executor.submit(
                extract_one,
                row,
            ): str(
                row["pdb_id"]
            ).lower()
            for row in ordered_manifest
        }

        completed = 0

        for future in as_completed(
            futures
        ):
            pdb_id = futures[
                future
            ]

            try:
                row = future.result()
            except Exception as exc:
                raise DownstreamMetadataProductionError(
                    f"Metadata extraction failed for {pdb_id}: "
                    f"{type(exc).__name__}: {exc}"
                ) from exc

            extracted.append(
                row
            )

            completed += 1

            if (
                completed % 250 == 0
                or completed
                == len(
                    ordered_manifest
                )
            ):
                print(
                    "Verified/extracted:",
                    f"{completed:,}/"
                    f"{len(ordered_manifest):,}",
                )

    extracted.sort(
        key=lambda row: row[
            "pdb_id"
        ]
    )

    if len(
        extracted
    ) != len(
        participating
    ):
        raise DownstreamMetadataProductionError(
            "Metadata extraction count mismatch"
        )

    if len(
        {
            row["pdb_id"]
            for row in extracted
        }
    ) != len(
        extracted
    ):
        raise DownstreamMetadataProductionError(
            "Duplicate extracted metadata entry"
        )

    # --------------------------------------------------------
    # Publish metadata only.
    # --------------------------------------------------------

    output_root = (
        stage_root
        / "downstream_metadata"
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
        output_root / "_SUCCESS"
    )

    if success_path.exists():
        success_path.unlink()

    output_path = (
        finalized
        / "entry_metadata.parquet"
    )

    temporary = (
        output_path.with_suffix(
            output_path.suffix
            + ".tmp"
        )
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
        != len(
            participating
        )
    ):
        temporary.unlink()

        raise DownstreamMetadataProductionError(
            "Published metadata row count mismatch"
        )

    temporary.replace(
        output_path
    )

    methods_present = sum(
        bool(
            row[
                "experimental_methods"
            ]
        )
        for row in extracted
    )

    refine_present = sum(
        bool(
            row[
                "refine_ls_d_res_high"
            ]
        )
        for row in extracted
    )

    em_present = sum(
        bool(
            row[
                "em_3d_reconstruction_resolution"
            ]
        )
        for row in extracted
    )

    group_present = sum(
        bool(
            row[
                "has_deposit_group"
            ]
        )
        for row in extracted
    )

    group_pandda = sum(
        bool(
            row[
                "deposit_group_mentions_pandda"
            ]
        )
        for row in extracted
    )

    entry_pandda = sum(
        bool(
            row[
                "entry_mentions_pandda"
            ]
        )
        for row in extracted
    )

    total_verified_bytes = sum(
        int(
            row[
                "size_bytes"
            ]
        )
        for row in extracted
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
        "metadata_extraction_pipeline_git_commit": (
            metadata_commit
        ),
        "config_sha256": (
            loaded.sha256
        ),
    }

    summary = {
        "summary_schema_name": (
            "pdbclean_downstream_metadata_global_summary"
        ),
        "summary_schema_version": "1.0",
        **provenance,
        "scientific_filtering_performed": False,
        "participating_deposition_count": (
            len(
                participating
            )
        ),
        "verified_source_object_count": (
            len(
                extracted
            )
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
            "pdbclean_downstream_metadata_success"
        ),
        "success_schema_version": "1.0",
        **provenance,
        "entry_metadata": (
            "finalized/entry_metadata.parquet"
        ),
        "global_summary": (
            "global_summary.json"
        ),
    }

    _write_json_atomic(
        success,
        success_path,
    )

    print()
    print(
        "===== DOWNSTREAM METADATA EXTRACTION ====="
    )
    print(
        "Depositions:",
        f"{len(extracted):,}",
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
        "Scientific filtering performed:",
        "NO",
    )
    print(
        "DOWNSTREAM METADATA EXTRACTION: PASS"
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
