"""Build the non-filtering Acta manual-review manifest.

The input is the frozen Acta downstream-investigation publication.

This stage performs NO scientific filtering.  It groups chain-pair
evidence by unordered PDB deposition pair so that the subsequent
manual evaluation is tractable.

The complete Stage-11 chain-pair publication remains the authoritative
underlying evidence.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from collections import defaultdict
from pathlib import Path
from typing import Any

import pyarrow as pa
import pyarrow.parquet as pq

from pdbclean.config import load_config


class ActaManualReviewManifestError(RuntimeError):
    """Raised when the manual-review manifest cannot be built safely."""


SUMMARY_SCHEMA_NAME = (
    "pdbclean_acta_manual_review_manifest_global_summary"
)
SUMMARY_SCHEMA_VERSION = "1.0"

SUCCESS_SCHEMA_NAME = (
    "pdbclean_acta_manual_review_manifest_success"
)
SUCCESS_SCHEMA_VERSION = "1.0"


MANIFEST_SCHEMA = pa.schema(
    [
        pa.field("deposition_pair_id", pa.string(), nullable=False),
        pa.field("left_pdb_id", pa.string(), nullable=False),
        pa.field("right_pdb_id", pa.string(), nullable=False),
        pa.field("chain_pair_count", pa.int64(), nullable=False),
        pa.field(
            "exact_bri_chain_pair_count",
            pa.int64(),
            nullable=False,
        ),
        pa.field(
            "nonzero_near_duplicate_chain_pair_count",
            pa.int64(),
            nullable=False,
        ),
        pa.field(
            "minimum_d_bri_mA",
            pa.int64(),
            nullable=False,
        ),
        pa.field(
            "maximum_d_bri_mA",
            pa.int64(),
            nullable=False,
        ),
        pa.field(
            "minimum_d_bri",
            pa.float64(),
            nullable=False,
        ),
        pa.field(
            "maximum_d_bri",
            pa.float64(),
            nullable=False,
        ),
        pa.field(
            "retained_residue_counts",
            pa.list_(pa.int64()),
            nullable=False,
        ),
        pa.field(
            "left_chain_count",
            pa.int64(),
            nullable=False,
        ),
        pa.field(
            "right_chain_count",
            pa.int64(),
            nullable=False,
        ),
        pa.field(
            "left_chain_ids",
            pa.list_(pa.string()),
            nullable=False,
        ),
        pa.field(
            "right_chain_ids",
            pa.list_(pa.string()),
            nullable=False,
        ),
        pa.field(
            "left_experimental_methods",
            pa.list_(pa.string()),
            nullable=False,
        ),
        pa.field(
            "right_experimental_methods",
            pa.list_(pa.string()),
            nullable=False,
        ),
        pa.field(
            "left_resolution_angstrom",
            pa.float64(),
            nullable=False,
        ),
        pa.field(
            "left_resolution_basis",
            pa.string(),
            nullable=False,
        ),
        pa.field(
            "right_resolution_angstrom",
            pa.float64(),
            nullable=False,
        ),
        pa.field(
            "right_resolution_basis",
            pa.string(),
            nullable=False,
        ),
        pa.field(
            "left_initial_deposition_date",
            pa.string(),
            nullable=True,
        ),
        pa.field(
            "right_initial_deposition_date",
            pa.string(),
            nullable=True,
        ),
        pa.field(
            "left_struct_title",
            pa.string(),
            nullable=False,
        ),
        pa.field(
            "right_struct_title",
            pa.string(),
            nullable=False,
        ),
        pa.field(
            "left_struct_keywords_text",
            pa.list_(pa.string()),
            nullable=False,
        ),
        pa.field(
            "right_struct_keywords_text",
            pa.list_(pa.string()),
            nullable=False,
        ),
        pa.field(
            "review_status",
            pa.string(),
            nullable=False,
        ),
        pa.field(
            "manual_exclusion_reason",
            pa.string(),
            nullable=True,
        ),
        pa.field(
            "manual_notes",
            pa.string(),
            nullable=True,
        ),
    ],
    metadata={
        b"schema_name": (
            b"pdbclean_acta_manual_review_manifest"
        ),
        b"schema_version": b"1.0",
        b"scientific_filtering": b"none",
        b"grouping_unit": (
            b"unordered_PDB_deposition_pair"
        ),
        b"underlying_evidence": (
            b"Stage-11 chain-pair publication"
        ),
    },
)


def _validate_commit(value: str) -> str:
    value = value.strip().lower()

    if re.fullmatch(r"[0-9a-f]{40}", value) is None:
        raise ActaManualReviewManifestError(
            "Expected full 40-character Git commit"
        )

    return value


def _read_json(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise ActaManualReviewManifestError(
            f"Missing JSON file: {path}"
        )

    value = json.loads(
        path.read_text(
            encoding="utf-8"
        )
    )

    if not isinstance(value, dict):
        raise ActaManualReviewManifestError(
            f"Expected JSON object: {path}"
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


def _canonical_deposition_pair(
    query_pdb_id: str,
    query_chain_id: str,
    subject_pdb_id: str,
    subject_chain_id: str,
) -> tuple[str, str, str, str]:
    q = query_pdb_id.lower()
    s = subject_pdb_id.lower()

    if q == s:
        raise ActaManualReviewManifestError(
            f"Same-deposition pair reached manual manifest: {q}"
        )

    if q < s:
        return (
            q,
            s,
            query_chain_id,
            subject_chain_id,
        )

    return (
        s,
        q,
        subject_chain_id,
        query_chain_id,
    )


def _metadata_map(
    path: Path,
) -> dict[str, dict[str, Any]]:
    table = pq.read_table(
        path,
        columns=[
            "pdb_id",
            "experimental_methods",
            "initial_deposition_date",
            "struct_title",
            "struct_keywords_text",
        ],
    )

    result: dict[
        str,
        dict[str, Any],
    ] = {}

    for row in table.to_pylist():
        pdb_id = str(
            row["pdb_id"]
        ).lower()

        if pdb_id in result:
            raise ActaManualReviewManifestError(
                f"Duplicate metadata row: {pdb_id}"
            )

        if not row["struct_title"]:
            raise ActaManualReviewManifestError(
                f"Missing title for {pdb_id}"
            )

        if not row[
            "struct_keywords_text"
        ]:
            raise ActaManualReviewManifestError(
                f"Missing keywords for {pdb_id}"
            )

        result[pdb_id] = row

    return result


def _write_review_csv(
    rows: list[dict[str, Any]],
    path: Path,
) -> None:
    temporary = path.with_suffix(
        path.suffix + ".tmp"
    )

    if temporary.exists():
        temporary.unlink()

    fields = [
        "deposition_pair_id",
        "left_pdb_id",
        "right_pdb_id",
        "chain_pair_count",
        "exact_bri_chain_pair_count",
        "nonzero_near_duplicate_chain_pair_count",
        "minimum_d_bri_mA",
        "maximum_d_bri_mA",
        "retained_residue_counts",
        "left_experimental_methods",
        "right_experimental_methods",
        "left_resolution_angstrom",
        "right_resolution_angstrom",
        "left_initial_deposition_date",
        "right_initial_deposition_date",
        "left_struct_title",
        "right_struct_title",
        "left_struct_keywords_text",
        "right_struct_keywords_text",
        "review_status",
        "manual_exclusion_reason",
        "manual_notes",
    ]

    with temporary.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
        )

        writer.writeheader()

        for row in rows:
            output = {
                field: row.get(field)
                for field in fields
            }

            for field in (
                "retained_residue_counts",
                "left_experimental_methods",
                "right_experimental_methods",
                "left_struct_keywords_text",
                "right_struct_keywords_text",
            ):
                output[field] = " | ".join(
                    str(value)
                    for value in output[field]
                )

            writer.writerow(output)

    temporary.replace(path)



def _publication_layout(
    stage11_version: str,
) -> tuple[str, str]:
    """Return Stage-11 input and manual-manifest output directories.

    v1 paths are preserved exactly for reproducibility.
    v2 is physically isolated from the frozen v1 publications.
    """

    layouts = {
        "1.0": (
            "acta_downstream_investigation",
            "acta_manual_review_manifest",
        ),
        "2.0": (
            "acta_downstream_investigation_v2",
            "acta_manual_review_manifest_v2",
        ),
    }

    try:
        return layouts[str(stage11_version)]
    except KeyError as exc:
        raise ActaManualReviewManifestError(
            f"Unsupported Stage-11 version: {stage11_version!r}"
        ) from exc


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

    parser.add_argument(
        "--stage11-version",
        choices=("1.0", "2.0"),
        default="1.0",
        help=(
            "Stage-11 publication to consume; "
            "defaults to frozen v1 for backward compatibility"
        ),
    )

    args = parser.parse_args()

    (
        stage11_directory,
        manifest_output_directory,
    ) = _publication_layout(
        args.stage11_version
    )

    pipeline_commit = _validate_commit(
        args.pipeline_git_commit
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
        config[
            "storage"
        ]["output_root"]
    )

    if not storage_root.is_absolute():
        storage_root = repo / storage_root

    stage11_roots = sorted(
        path.parent
        for path in storage_root.glob(
            f"*/{protocol}/"
            f"{stage11_directory}/_SUCCESS"
        )
    )

    if len(stage11_roots) != 1:
        raise ActaManualReviewManifestError(
            "Expected exactly one completed Stage-11 publication"
        )

    stage11_root = stage11_roots[0]
    stage_root = stage11_root.parent

    metadata_root = (
        stage_root
        / "downstream_metadata"
    )

    stage11_success = _read_json(
        stage11_root
        / "_SUCCESS"
    )
    stage11_summary = _read_json(
        stage11_root
        / "global_summary.json"
    )
    metadata_success = _read_json(
        metadata_root
        / "_SUCCESS"
    )

    if (
        stage11_success.get(
            "success_schema_name"
        )
        != "pdbclean_acta_downstream_investigation_success"
    ):
        raise ActaManualReviewManifestError(
            "Unexpected Stage-11 success schema"
        )

    if (
        stage11_summary.get(
            "summary_schema_name"
        )
        != "pdbclean_acta_downstream_investigation_global_summary"
    ):
        raise ActaManualReviewManifestError(
            "Unexpected Stage-11 summary schema"
        )

    if not (
        stage11_summary.get(
            "pair_accounting_valid"
        )
        and stage11_summary.get(
            "preflight_oracle_counts_equal"
        )
        and stage11_summary.get(
            "processing_error_count"
        )
        == 0
    ):
        raise ActaManualReviewManifestError(
            "Stage-11 validation gate failed"
        )

    if (
        stage11_summary.get(
            "manual_evaluation_performed"
        )
        is not False
        or stage11_summary.get(
            "automatic_virus_filtering"
        )
        is not False
        or stage11_summary.get(
            "automatic_ribosome_filtering"
        )
        is not False
    ):
        raise ActaManualReviewManifestError(
            "Unexpected manual-review state in Stage 11"
        )

    if (
        stage11_success.get(
            "canonical_config_sha256"
        )
        != loaded.sha256
    ):
        raise ActaManualReviewManifestError(
            "Canonical config provenance mismatch"
        )

    pair_path = (
        stage11_root
        / stage11_success[
            "detailed_inspection_candidates"
        ]
    )

    metadata_path = (
        metadata_root
        / metadata_success[
            "entry_metadata"
        ]
    )

    metadata = _metadata_map(
        metadata_path
    )

    groups: dict[
        tuple[str, str],
        dict[str, Any],
    ] = defaultdict(
        lambda: {
            "chain_pair_count": 0,
            "exact_count": 0,
            "minimum_d_bri_mA": None,
            "maximum_d_bri_mA": None,
            "retained_residue_counts": set(),
            "left_chain_ids": set(),
            "right_chain_ids": set(),
            "left_resolution_angstrom": None,
            "left_resolution_basis": None,
            "right_resolution_angstrom": None,
            "right_resolution_basis": None,
        }
    )

    pf = pq.ParquetFile(
        pair_path
    )

    input_chain_pairs = 0

    for batch in pf.iter_batches(
        columns=[
            "query_pdb_id",
            "query_label_chain_id",
            "subject_pdb_id",
            "subject_label_chain_id",
            "retained_residue_count",
            "d_bri_mA",
            "query_resolution_angstrom",
            "query_resolution_basis",
            "subject_resolution_angstrom",
            "subject_resolution_basis",
        ],
        batch_size=65_536,
    ):
        for row in pa.Table.from_batches(
            [batch]
        ).to_pylist():
            input_chain_pairs += 1

            q = str(
                row["query_pdb_id"]
            ).lower()
            s = str(
                row["subject_pdb_id"]
            ).lower()

            (
                left,
                right,
                left_chain,
                right_chain,
            ) = _canonical_deposition_pair(
                q,
                str(
                    row[
                        "query_label_chain_id"
                    ]
                ),
                s,
                str(
                    row[
                        "subject_label_chain_id"
                    ]
                ),
            )

            if left not in metadata:
                raise ActaManualReviewManifestError(
                    f"Missing metadata for {left}"
                )

            if right not in metadata:
                raise ActaManualReviewManifestError(
                    f"Missing metadata for {right}"
                )

            if q == left:
                left_resolution = float(
                    row[
                        "query_resolution_angstrom"
                    ]
                )
                left_basis = str(
                    row[
                        "query_resolution_basis"
                    ]
                )
                right_resolution = float(
                    row[
                        "subject_resolution_angstrom"
                    ]
                )
                right_basis = str(
                    row[
                        "subject_resolution_basis"
                    ]
                )
            else:
                left_resolution = float(
                    row[
                        "subject_resolution_angstrom"
                    ]
                )
                left_basis = str(
                    row[
                        "subject_resolution_basis"
                    ]
                )
                right_resolution = float(
                    row[
                        "query_resolution_angstrom"
                    ]
                )
                right_basis = str(
                    row[
                        "query_resolution_basis"
                    ]
                )

            g = groups[
                (left, right)
            ]

            g[
                "chain_pair_count"
            ] += 1

            distance = int(
                row["d_bri_mA"]
            )

            if distance == 0:
                g["exact_count"] += 1

            minimum = g[
                "minimum_d_bri_mA"
            ]
            maximum = g[
                "maximum_d_bri_mA"
            ]

            if (
                minimum is None
                or distance < minimum
            ):
                g[
                    "minimum_d_bri_mA"
                ] = distance

            if (
                maximum is None
                or distance > maximum
            ):
                g[
                    "maximum_d_bri_mA"
                ] = distance

            g[
                "retained_residue_counts"
            ].add(
                int(
                    row[
                        "retained_residue_count"
                    ]
                )
            )

            g[
                "left_chain_ids"
            ].add(
                left_chain
            )
            g[
                "right_chain_ids"
            ].add(
                right_chain
            )

            for field, value in (
                (
                    "left_resolution_angstrom",
                    left_resolution,
                ),
                (
                    "left_resolution_basis",
                    left_basis,
                ),
                (
                    "right_resolution_angstrom",
                    right_resolution,
                ),
                (
                    "right_resolution_basis",
                    right_basis,
                ),
            ):
                previous = g[field]

                if previous is None:
                    g[field] = value
                elif previous != value:
                    raise ActaManualReviewManifestError(
                        "Inconsistent Stage-11 resolution metadata "
                        f"for deposition pair {left}/{right}"
                    )

    expected_chain_pairs = int(
        stage11_summary[
            "detailed_inspection_candidate_pair_count"
        ]
    )

    if input_chain_pairs != expected_chain_pairs:
        raise ActaManualReviewManifestError(
            "Stage-11 chain-pair accounting mismatch"
        )

    manifest_rows: list[
        dict[str, Any]
    ] = []

    for (
        left,
        right,
    ), g in sorted(
        groups.items()
    ):
        left_meta = metadata[left]
        right_meta = metadata[right]

        count = int(
            g[
                "chain_pair_count"
            ]
        )
        exact = int(
            g["exact_count"]
        )

        minimum_mA = int(
            g[
                "minimum_d_bri_mA"
            ]
        )
        maximum_mA = int(
            g[
                "maximum_d_bri_mA"
            ]
        )

        manifest_rows.append(
            {
                "deposition_pair_id": (
                    f"{left}__{right}"
                ),
                "left_pdb_id": left,
                "right_pdb_id": right,
                "chain_pair_count": count,
                "exact_bri_chain_pair_count": exact,
                "nonzero_near_duplicate_chain_pair_count": (
                    count - exact
                ),
                "minimum_d_bri_mA": minimum_mA,
                "maximum_d_bri_mA": maximum_mA,
                "minimum_d_bri": (
                    minimum_mA
                    / 1000.0
                ),
                "maximum_d_bri": (
                    maximum_mA
                    / 1000.0
                ),
                "retained_residue_counts": sorted(
                    g[
                        "retained_residue_counts"
                    ]
                ),
                "left_chain_count": len(
                    g[
                        "left_chain_ids"
                    ]
                ),
                "right_chain_count": len(
                    g[
                        "right_chain_ids"
                    ]
                ),
                "left_chain_ids": sorted(
                    g[
                        "left_chain_ids"
                    ]
                ),
                "right_chain_ids": sorted(
                    g[
                        "right_chain_ids"
                    ]
                ),
                "left_experimental_methods": list(
                    left_meta[
                        "experimental_methods"
                    ]
                ),
                "right_experimental_methods": list(
                    right_meta[
                        "experimental_methods"
                    ]
                ),
                "left_resolution_angstrom": float(
                    g[
                        "left_resolution_angstrom"
                    ]
                ),
                "left_resolution_basis": str(
                    g[
                        "left_resolution_basis"
                    ]
                ),
                "right_resolution_angstrom": float(
                    g[
                        "right_resolution_angstrom"
                    ]
                ),
                "right_resolution_basis": str(
                    g[
                        "right_resolution_basis"
                    ]
                ),
                "left_initial_deposition_date": (
                    left_meta[
                        "initial_deposition_date"
                    ]
                ),
                "right_initial_deposition_date": (
                    right_meta[
                        "initial_deposition_date"
                    ]
                ),
                "left_struct_title": str(
                    left_meta[
                        "struct_title"
                    ]
                ),
                "right_struct_title": str(
                    right_meta[
                        "struct_title"
                    ]
                ),
                "left_struct_keywords_text": list(
                    left_meta[
                        "struct_keywords_text"
                    ]
                ),
                "right_struct_keywords_text": list(
                    right_meta[
                        "struct_keywords_text"
                    ]
                ),
                "review_status": "unreviewed",
                "manual_exclusion_reason": None,
                "manual_notes": None,
            }
        )

    if (
        sum(
            row[
                "chain_pair_count"
            ]
            for row in manifest_rows
        )
        != input_chain_pairs
    ):
        raise ActaManualReviewManifestError(
            "Grouped chain-pair accounting failed"
        )

    participating = {
        pdb_id
        for row in manifest_rows
        for pdb_id in (
            row["left_pdb_id"],
            row["right_pdb_id"],
        )
    }

    if len(
        participating
    ) != int(
        stage11_summary[
            "detailed_inspection_participating_deposition_count"
        ]
    ):
        raise ActaManualReviewManifestError(
            "Participating deposition accounting mismatch"
        )

    output_root = (
        stage_root
        / manifest_output_directory
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

    parquet_path = (
        finalized
        / "deposition_pair_manifest.parquet"
    )

    parquet_tmp = parquet_path.with_suffix(
        parquet_path.suffix + ".tmp"
    )

    if parquet_tmp.exists():
        parquet_tmp.unlink()

    table = pa.Table.from_pylist(
        manifest_rows,
        schema=MANIFEST_SCHEMA,
    )

    pq.write_table(
        table,
        parquet_tmp,
        compression="zstd",
        version="2.6",
        use_dictionary=True,
    )

    if (
        pq.read_metadata(
            parquet_tmp
        ).num_rows
        != len(manifest_rows)
    ):
        parquet_tmp.unlink()

        raise ActaManualReviewManifestError(
            "Manifest Parquet validation failed"
        )

    parquet_tmp.replace(
        parquet_path
    )

    csv_path = (
        finalized
        / "manual_review_template.csv"
    )

    _write_review_csv(
        manifest_rows,
        csv_path,
    )

    exact_groups = sum(
        row[
            "exact_bri_chain_pair_count"
        ]
        > 0
        for row in manifest_rows
    )

    provenance = {
        "source_stage11_version": str(
            args.stage11_version
        ),
        "source_stage11_publication_directory": (
            stage11_directory
        ),
        "manual_review_manifest_publication_directory": (
            manifest_output_directory
        ),
        "snapshot": (
            stage11_success[
                "snapshot"
            ]
        ),
        "cleaning_protocol": (
            stage11_success[
                "cleaning_protocol"
            ]
        ),
        "canonical_config_sha256": (
            loaded.sha256
        ),
        "downstream_config_sha256": (
            stage11_success[
                "downstream_config_sha256"
            ]
        ),
        "full_bri_nn_pipeline_git_commit": (
            stage11_success[
                "full_bri_nn_pipeline_git_commit"
            ]
        ),
        "metadata_finalizer_git_commit": (
            stage11_success[
                "metadata_finalizer_git_commit"
            ]
        ),
        "acta_downstream_pipeline_git_commit": (
            stage11_success[
                "acta_downstream_pipeline_git_commit"
            ]
        ),
        "manual_review_manifest_pipeline_git_commit": (
            pipeline_commit
        ),
    }

    summary = {
        "summary_schema_name": (
            SUMMARY_SCHEMA_NAME
        ),
        "summary_schema_version": (
            SUMMARY_SCHEMA_VERSION
        ),
        **provenance,
        "scientific_filtering_performed": False,
        "manual_evaluation_performed": False,
        "grouping_unit": (
            "unordered_PDB_deposition_pair"
        ),
        "input_chain_pair_count": (
            input_chain_pairs
        ),
        "deposition_pair_count": (
            len(manifest_rows)
        ),
        "participating_deposition_count": (
            len(participating)
        ),
        "deposition_pairs_with_exact_bri_evidence": (
            exact_groups
        ),
        "deposition_pairs_with_only_nonzero_near_duplicate_evidence": (
            len(manifest_rows)
            - exact_groups
        ),
        "all_review_status_unreviewed": True,
        "underlying_chain_pair_accounting_valid": True,
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
        "deposition_pair_manifest": (
            "finalized/"
            "deposition_pair_manifest.parquet"
        ),
        "manual_review_template": (
            "finalized/"
            "manual_review_template.csv"
        ),
        "global_summary": (
            "global_summary.json"
        ),
    }

    # Success marker is written last.
    _write_json_atomic(
        success,
        success_path,
    )

    print(
        "===== ACTA MANUAL-REVIEW MANIFEST ====="
    )
    print(
        "Underlying Stage-11 chain pairs:",
        f"{input_chain_pairs:,}",
    )
    print(
        "Unique unordered deposition pairs:",
        f"{len(manifest_rows):,}",
    )
    print(
        "Participating depositions:",
        f"{len(participating):,}",
    )
    print(
        "Pairs containing exact BRI evidence:",
        f"{exact_groups:,}",
    )
    print(
        "Pairs with only nonzero <=10 mA evidence:",
        f"{len(manifest_rows) - exact_groups:,}",
    )
    print(
        "Scientific filtering performed: NO"
    )
    print(
        "Manual evaluation performed: NO"
    )
    print(
        "Underlying chain-pair accounting: PASS"
    )
    print(
        "ACTA MANUAL-REVIEW MANIFEST PUBLICATION: PASS"
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
