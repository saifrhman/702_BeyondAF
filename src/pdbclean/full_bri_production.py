"""Stage 8 production: exact full-BRI L-infinity comparisons."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

from pdbclean.config import load_config
from pdbclean.full_bri_compare import (
    bri_to_integer_mA,
)
from pdbclean.length_buckets import (
    validate_stage5_brain_publication,
)


class FullBRIProductionError(RuntimeError):
    """Raised when Stage-8 production cannot proceed safely."""


CANDIDATE_COMPARISON_SCHEMA = pa.schema(
    [
        pa.field("query_snapshot", pa.string(), nullable=False),
        pa.field("query_pdb_id", pa.string(), nullable=False),
        pa.field("query_model_id", pa.int64(), nullable=False),
        pa.field("query_label_chain_id", pa.string(), nullable=False),

        pa.field("subject_snapshot", pa.string(), nullable=False),
        pa.field("subject_pdb_id", pa.string(), nullable=False),
        pa.field("subject_model_id", pa.int64(), nullable=False),
        pa.field("subject_label_chain_id", pa.string(), nullable=False),

        pa.field(
            "retained_residue_count",
            pa.int64(),
            nullable=False,
        ),
        pa.field(
            "d_brain_numerator_max",
            pa.int64(),
            nullable=False,
        ),
        pa.field(
            "d_brain",
            pa.float64(),
            nullable=False,
        ),
        pa.field(
            "d_bri_mA",
            pa.int64(),
            nullable=False,
        ),
        pa.field(
            "d_bri",
            pa.float64(),
            nullable=False,
        ),
    ],
    metadata={
        b"schema_name": (
            b"pdbclean_stage8_full_bri_candidate_comparisons"
        ),
        b"schema_version": b"1.0",
        b"metric": b"full BRI L-infinity",
        b"distance_representation": (
            b"exact integer milliangstrom"
        ),
    },
)


M1_COMPARISON_SCHEMA = pa.schema(
    [
        pa.field("query_snapshot", pa.string(), nullable=False),
        pa.field("query_pdb_id", pa.string(), nullable=False),
        pa.field("query_model_id", pa.int64(), nullable=False),
        pa.field("query_label_chain_id", pa.string(), nullable=False),

        pa.field("subject_snapshot", pa.string(), nullable=False),
        pa.field("subject_pdb_id", pa.string(), nullable=False),
        pa.field("subject_model_id", pa.int64(), nullable=False),
        pa.field("subject_label_chain_id", pa.string(), nullable=False),

        pa.field(
            "retained_residue_count",
            pa.int64(),
            nullable=False,
        ),
        pa.field(
            "d_bri_mA",
            pa.int64(),
            nullable=False,
        ),
        pa.field(
            "d_bri",
            pa.float64(),
            nullable=False,
        ),
    ],
    metadata={
        b"schema_name": (
            b"pdbclean_stage8_full_bri_m1_comparisons"
        ),
        b"schema_version": b"1.0",
        b"metric": b"full BRI L-infinity",
        b"search_mode": b"m1_direct_all_pairs",
        b"distance_representation": (
            b"exact integer milliangstrom"
        ),
    },
)


def _read_json(path: Path) -> dict:
    try:
        value = json.loads(
            path.read_text(encoding="utf-8")
        )
    except Exception as exc:
        raise FullBRIProductionError(
            f"Cannot read JSON {path}: {exc}"
        ) from exc

    if not isinstance(value, dict):
        raise FullBRIProductionError(
            f"JSON artifact is not an object: {path}"
        )

    return value


def _validate_commit(value: str) -> str:
    if (
        not isinstance(value, str)
        or len(value) != 40
        or any(
            c not in "0123456789abcdef"
            for c in value
        )
    ):
        raise FullBRIProductionError(
            "Stage-8 pipeline commit must be a "
            "40-character lowercase hexadecimal SHA"
        )

    return value


def _write_json_atomic(
    payload: dict,
    path: Path,
) -> Path:
    path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    temporary = path.with_suffix(
        path.suffix + ".tmp"
    )

    temporary.write_text(
        json.dumps(
            payload,
            sort_keys=True,
            indent=2,
            ensure_ascii=True,
        )
        + "\n",
        encoding="utf-8",
    )

    temporary.replace(path)

    return path


def _identity(
    row: dict,
) -> tuple[str, str, int, str]:
    return (
        row["snapshot"],
        row["pdb_id"],
        row["model_id"],
        row["label_chain_id"],
    )


def _candidate_identity_indices(
    table: pa.Table,
    *,
    prefix: str,
    identity_to_index: dict[
        tuple[str, str, int, str],
        int,
    ],
) -> np.ndarray:
    snapshots = table[
        f"{prefix}_snapshot"
    ].to_pylist()
    pdb_ids = table[
        f"{prefix}_pdb_id"
    ].to_pylist()
    model_ids = table[
        f"{prefix}_model_id"
    ].to_pylist()
    chains = table[
        f"{prefix}_label_chain_id"
    ].to_pylist()

    result = np.empty(
        table.num_rows,
        dtype=np.int64,
    )

    for i in range(table.num_rows):
        identity = (
            snapshots[i],
            pdb_ids[i],
            model_ids[i],
            chains[i],
        )

        try:
            result[i] = identity_to_index[
                identity
            ]
        except KeyError as exc:
            raise FullBRIProductionError(
                f"Stage-7 {prefix} identity absent "
                f"from canonical BRI: {identity!r}"
            ) from exc

    return result


def _batched_full_bri_distances(
    bri: np.ndarray,
    query_indices: np.ndarray,
    subject_indices: np.ndarray,
    *,
    m: int,
) -> np.ndarray:
    """Vectorized exact integer-mA L-infinity comparisons."""

    if (
        query_indices.shape
        != subject_indices.shape
    ):
        raise FullBRIProductionError(
            "Query/subject index shape mismatch"
        )

    count = query_indices.shape[0]

    result = np.empty(
        count,
        dtype=np.int64,
    )

    # Bound temporary advanced-index arrays. This does not
    # affect the mathematics, only memory usage.
    target_coordinate_elements = 4_000_000

    chunk_size = max(
        1,
        min(
            65_536,
            target_coordinate_elements
            // (m * 9),
        ),
    )

    for start in range(
        0,
        count,
        chunk_size,
    ):
        stop = min(
            count,
            start + chunk_size,
        )

        q = bri[
            query_indices[start:stop]
        ]
        s = bri[
            subject_indices[start:stop]
        ]

        result[start:stop] = np.max(
            np.abs(q - s),
            axis=(1, 2),
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

    stage8_commit = _validate_commit(
        args.pipeline_git_commit
    )

    repo = Path.cwd()

    config = load_config(
        args.config
    ).data

    protocol = (
        config["release"]["protocol_version"]
    )

    storage_root = Path(
        config["storage"]["output_root"]
    )

    if not storage_root.is_absolute():
        storage_root = (
            repo / storage_root
        )

    completed_brain = sorted(
        path.parent
        for path in storage_root.glob(
            f"*/{protocol}/brain/_SUCCESS"
        )
    )

    if len(completed_brain) != 1:
        raise FullBRIProductionError(
            "Expected exactly one completed "
            "Stage-5 Brain publication"
        )

    stage5 = validate_stage5_brain_publication(
        completed_brain[0]
    )

    stage_root = completed_brain[0].parent

    # --------------------------------------------------------
    # Validate canonical Stage-7 input.
    # --------------------------------------------------------

    stage7_root = (
        stage_root / "brain_prefilter"
    )

    stage7_success = _read_json(
        stage7_root / "_SUCCESS"
    )
    stage7_summary = _read_json(
        stage7_root / "global_summary.json"
    )

    if (
        stage7_success.get(
            "success_schema_name"
        )
        != "pdbclean_stage7_brain_prefilter_success"
        or stage7_success.get(
            "success_schema_version"
        )
        != "1.0"
    ):
        raise FullBRIProductionError(
            "Unexpected Stage-7 _SUCCESS schema"
        )

    if (
        stage7_summary.get(
            "summary_schema_name"
        )
        != "pdbclean_stage7_brain_prefilter_global_summary"
        or stage7_summary.get(
            "summary_schema_version"
        )
        != "1.0"
    ):
        raise FullBRIProductionError(
            "Unexpected Stage-7 global-summary schema"
        )

    if (
        stage7_success["snapshot"]
        != stage5.snapshot
        or stage7_summary["snapshot"]
        != stage5.snapshot
    ):
        raise FullBRIProductionError(
            "Stage-7 snapshot mismatch"
        )

    if (
        stage7_success["cleaning_protocol"]
        != stage5.cleaning_protocol
        or stage7_summary["cleaning_protocol"]
        != stage5.cleaning_protocol
    ):
        raise FullBRIProductionError(
            "Stage-7 cleaning-protocol mismatch"
        )

    for field, expected in (
        stage5.provenance.items()
    ):
        if (
            stage7_success.get(field)
            != expected
            or stage7_summary.get(field)
            != expected
        ):
            raise FullBRIProductionError(
                f"Stage-7 provenance mismatch: {field}"
            )

    stage7_commit = stage7_success.get(
        "brain_prefilter_pipeline_git_commit"
    )

    if (
        stage7_commit
        != stage7_summary.get(
            "brain_prefilter_pipeline_git_commit"
        )
    ):
        raise FullBRIProductionError(
            "Stage-7 producer provenance mismatch"
        )

    candidate_path = (
        stage7_root
        / stage7_success["candidate_pairs"]
    )
    bypass_path = (
        stage7_root
        / stage7_success["m1_bypass"]
    )
    bucket_path = (
        stage7_root
        / stage7_success["bucket_summary"]
    )

    candidate_count = pq.read_metadata(
        candidate_path
    ).num_rows
    bypass_count = pq.read_metadata(
        bypass_path
    ).num_rows

    if (
        candidate_count
        != stage7_summary[
            "candidate_pair_count"
        ]
    ):
        raise FullBRIProductionError(
            "Stage-7 candidate-count mismatch"
        )

    if (
        bypass_count
        != stage7_summary[
            "m1_bypass_chain_count"
        ]
    ):
        raise FullBRIProductionError(
            "Stage-7 m1-bypass count mismatch"
        )

    if (
        stage7_summary[
            "processing_error_count"
        ]
        != 0
    ):
        raise FullBRIProductionError(
            "Stage 8 requires zero Stage-7 errors"
        )

    buckets = pq.read_table(
        bucket_path
    ).to_pylist()

    required_counts: dict[int, int] = {}

    for row in buckets:
        m = row[
            "retained_residue_count"
        ]

        if m == 1:
            required_counts[1] = row[
                "chain_count"
            ]
            continue

        if row[
            "candidate_pair_count"
        ] > 0:
            required_counts[m] = row[
                "chain_count"
            ]

    if required_counts.get(1) != bypass_count:
        raise FullBRIProductionError(
            "Stage-7 m1 bucket accounting mismatch"
        )

    print(
        "Required BRI length buckets:",
        f"{len(required_counts):,}",
    )

    # --------------------------------------------------------
    # Load canonical Stage-3 BRI exactly once.
    # --------------------------------------------------------

    bri_path = (
        stage_root
        / "bri/finalized/bri.parquet"
    )

    if not bri_path.is_file():
        raise FullBRIProductionError(
            "Canonical Stage-3 BRI population is missing"
        )

    if (
        pq.read_metadata(
            bri_path
        ).num_rows
        != stage5.input_bri_chain_count
    ):
        raise FullBRIProductionError(
            "Canonical Stage-3 BRI population count mismatch"
        )

    bri_by_m = {
        m: np.empty(
            (count, m, 9),
            dtype=np.int64,
        )
        for m, count
        in required_counts.items()
    }

    identity_to_index = {
        m: {}
        for m in required_counts
    }

    fill_count = {
        m: 0
        for m in required_counts
    }

    m1_identity_by_index = [
        None
    ] * required_counts[1]

    bri_pf = pq.ParquetFile(
        bri_path
    )

    print(
        "Loading required canonical BRI chains..."
    )

    for batch in bri_pf.iter_batches(
        columns=[
            "snapshot",
            "pdb_id",
            "model_id",
            "label_chain_id",
            "retained_residue_count",
            "bri",
        ],
        batch_size=256,
    ):
        for row in batch.to_pylist():
            m = row[
                "retained_residue_count"
            ]

            if m not in required_counts:
                continue

            index = fill_count[m]

            if index >= required_counts[m]:
                raise FullBRIProductionError(
                    f"Too many canonical BRI chains for m={m}"
                )

            identity = _identity(
                row
            )

            if (
                identity
                in identity_to_index[m]
            ):
                raise FullBRIProductionError(
                    "Duplicate canonical BRI identity: "
                    f"{identity!r}"
                )

            integer_bri = (
                bri_to_integer_mA(
                    np.asarray(
                        row["bri"],
                        dtype=np.float64,
                    ),
                    name=(
                        "stage8_canonical_bri"
                    ),
                )
            )

            if integer_bri.shape != (
                m,
                9,
            ):
                raise FullBRIProductionError(
                    f"Canonical BRI shape mismatch for "
                    f"{identity!r}"
                )

            bri_by_m[m][
                index
            ] = integer_bri

            identity_to_index[m][
                identity
            ] = index

            if m == 1:
                m1_identity_by_index[
                    index
                ] = identity

            fill_count[m] += 1

    for m, expected in (
        required_counts.items()
    ):
        if fill_count[m] != expected:
            raise FullBRIProductionError(
                f"Canonical BRI load mismatch for m={m}: "
                f"{fill_count[m]} != {expected}"
            )

    print(
        "Canonical BRI load: PASS"
    )

    # Verify Stage-7 m=1 bypass identities exactly.
    bypass = pq.read_table(
        bypass_path,
        columns=[
            "snapshot",
            "pdb_id",
            "model_id",
            "label_chain_id",
        ],
    ).to_pylist()

    bypass_identities = {
        (
            row["snapshot"],
            row["pdb_id"],
            row["model_id"],
            row["label_chain_id"],
        )
        for row in bypass
    }

    if (
        bypass_identities
        != set(
            identity_to_index[1]
        )
    ):
        raise FullBRIProductionError(
            "Stage-7 m1 bypass identities do not "
            "exactly match canonical m1 BRI identities"
        )

    # --------------------------------------------------------
    # Publish Stage-7-candidate full-BRI comparisons.
    # --------------------------------------------------------

    output_root = (
        stage_root / "full_bri_compare"
    )
    finalized = (
        output_root / "finalized"
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

    candidate_output = (
        finalized
        / "candidate_comparisons.parquet"
    )
    candidate_tmp = candidate_output.with_suffix(
        candidate_output.suffix
        + ".tmp"
    )

    if candidate_tmp.exists():
        candidate_tmp.unlink()

    writer = pq.ParquetWriter(
        candidate_tmp,
        CANDIDATE_COMPARISON_SCHEMA,
        compression="zstd",
        version="2.6",
        use_dictionary=True,
    )

    candidate_pf = pq.ParquetFile(
        candidate_path
    )

    candidate_rows_written = 0
    seen_candidate_m = set()

    try:
        for rg in range(
            candidate_pf.metadata.num_row_groups
        ):
            table = candidate_pf.read_row_group(
                rg
            )

            if table.num_rows == 0:
                continue

            m_values = set(
                table[
                    "retained_residue_count"
                ].to_pylist()
            )

            if len(m_values) != 1:
                raise FullBRIProductionError(
                    "Stage-7 candidate row group "
                    "contains mixed m values"
                )

            m = next(
                iter(m_values)
            )

            if m < 2:
                raise FullBRIProductionError(
                    "Stage-7 Brain candidate contains m < 2"
                )

            if m in seen_candidate_m:
                raise FullBRIProductionError(
                    f"Stage-7 candidate bucket m={m} "
                    "appears in multiple row groups"
                )

            seen_candidate_m.add(
                m
            )

            if m not in bri_by_m:
                raise FullBRIProductionError(
                    f"No canonical BRI cache for candidate m={m}"
                )

            q_indices = (
                _candidate_identity_indices(
                    table,
                    prefix="query",
                    identity_to_index=(
                        identity_to_index[m]
                    ),
                )
            )

            s_indices = (
                _candidate_identity_indices(
                    table,
                    prefix="subject",
                    identity_to_index=(
                        identity_to_index[m]
                    ),
                )
            )

            if np.any(
                q_indices == s_indices
            ):
                raise FullBRIProductionError(
                    f"Self-pair found in Stage-7 candidates for m={m}"
                )

            d_bri_mA = (
                _batched_full_bri_distances(
                    bri_by_m[m],
                    q_indices,
                    s_indices,
                    m=m,
                )
            )

            d_num = np.asarray(
                table[
                    "d_brain_numerator_max"
                ]
                .combine_chunks()
                .to_numpy(
                    zero_copy_only=False
                ),
                dtype=np.int64,
            )

            # Exact integer form of Lemma 6.1:
            #
            # dBrain = d_num / (1000*(m-1))
            # dBRI   = d_bri_mA / 1000
            #
            # therefore dBrain <= dBRI iff
            # d_num <= d_bri_mA * (m-1).
            if np.any(
                d_num
                > d_bri_mA
                * (m - 1)
            ):
                raise FullBRIProductionError(
                    f"Lemma 6.1 violation in real "
                    f"Stage-8 bucket m={m}"
                )

            d_brain = np.asarray(
                table[
                    "d_brain"
                ]
                .combine_chunks()
                .to_numpy(
                    zero_copy_only=False
                ),
                dtype=np.float64,
            )

            expected_d_brain = (
                d_num.astype(
                    np.float64
                )
                / (
                    1000.0
                    * float(m - 1)
                )
            )

            if not np.array_equal(
                d_brain,
                expected_d_brain,
            ):
                raise FullBRIProductionError(
                    f"Stage-7 exact dBrain representation "
                    f"mismatch for m={m}"
                )

            d_bri = (
                d_bri_mA.astype(
                    np.float64
                )
                / 1000.0
            )

            output = pa.Table.from_arrays(
                [
                    table["query_snapshot"],
                    table["query_pdb_id"],
                    table["query_model_id"],
                    table["query_label_chain_id"],

                    table["subject_snapshot"],
                    table["subject_pdb_id"],
                    table["subject_model_id"],
                    table["subject_label_chain_id"],

                    table["retained_residue_count"],
                    table["d_brain_numerator_max"],
                    table["d_brain"],

                    pa.array(
                        d_bri_mA,
                        type=pa.int64(),
                    ),
                    pa.array(
                        d_bri,
                        type=pa.float64(),
                    ),
                ],
                schema=(
                    CANDIDATE_COMPARISON_SCHEMA
                ),
            )

            writer.write_table(
                output
            )

            candidate_rows_written += (
                output.num_rows
            )

            print(
                f"m={m:5d} "
                f"pairs={output.num_rows:10,d}"
            )

    finally:
        writer.close()

    if (
        candidate_rows_written
        != candidate_count
    ):
        if candidate_tmp.exists():
            candidate_tmp.unlink()

        raise FullBRIProductionError(
            "Stage-8 candidate comparison "
            "population accounting failed"
        )

    if (
        pq.read_metadata(
            candidate_tmp
        ).num_rows
        != candidate_count
    ):
        candidate_tmp.unlink()

        raise FullBRIProductionError(
            "Stage-8 candidate Parquet count mismatch"
        )

    candidate_tmp.replace(
        candidate_output
    )

    # --------------------------------------------------------
    # m=1 bypass: direct all-pairs full-BRI comparison.
    # --------------------------------------------------------

    m1_count = required_counts[1]

    query_indices, subject_indices = (
        np.triu_indices(
            m1_count,
            k=1,
        )
    )

    query_indices = query_indices.astype(
        np.int64,
        copy=False,
    )
    subject_indices = (
        subject_indices.astype(
            np.int64,
            copy=False,
        )
    )

    expected_m1_pairs = (
        m1_count
        * (m1_count - 1)
        // 2
    )

    if (
        query_indices.shape[0]
        != expected_m1_pairs
    ):
        raise FullBRIProductionError(
            "m1 direct-pair generation failed"
        )

    m1_d_bri_mA = (
        _batched_full_bri_distances(
            bri_by_m[1],
            query_indices,
            subject_indices,
            m=1,
        )
    )

    m1_d_bri = (
        m1_d_bri_mA.astype(
            np.float64
        )
        / 1000.0
    )

    q_identity = [
        m1_identity_by_index[
            int(i)
        ]
        for i in query_indices
    ]
    s_identity = [
        m1_identity_by_index[
            int(i)
        ]
        for i in subject_indices
    ]

    m1_table = pa.Table.from_arrays(
        [
            pa.array(
                [x[0] for x in q_identity],
                type=pa.string(),
            ),
            pa.array(
                [x[1] for x in q_identity],
                type=pa.string(),
            ),
            pa.array(
                [x[2] for x in q_identity],
                type=pa.int64(),
            ),
            pa.array(
                [x[3] for x in q_identity],
                type=pa.string(),
            ),

            pa.array(
                [x[0] for x in s_identity],
                type=pa.string(),
            ),
            pa.array(
                [x[1] for x in s_identity],
                type=pa.string(),
            ),
            pa.array(
                [x[2] for x in s_identity],
                type=pa.int64(),
            ),
            pa.array(
                [x[3] for x in s_identity],
                type=pa.string(),
            ),

            pa.array(
                np.ones(
                    expected_m1_pairs,
                    dtype=np.int64,
                )
            ),
            pa.array(
                m1_d_bri_mA,
                type=pa.int64(),
            ),
            pa.array(
                m1_d_bri,
                type=pa.float64(),
            ),
        ],
        schema=M1_COMPARISON_SCHEMA,
    )

    m1_output = (
        finalized
        / "m1_comparisons.parquet"
    )
    m1_tmp = m1_output.with_suffix(
        m1_output.suffix
        + ".tmp"
    )

    if m1_tmp.exists():
        m1_tmp.unlink()

    pq.write_table(
        m1_table,
        m1_tmp,
        compression="zstd",
        version="2.6",
        use_dictionary=True,
    )

    if (
        pq.read_metadata(
            m1_tmp
        ).num_rows
        != expected_m1_pairs
    ):
        m1_tmp.unlink()

        raise FullBRIProductionError(
            "Stage-8 m1 comparison count mismatch"
        )

    m1_tmp.replace(
        m1_output
    )

    total_comparisons = (
        candidate_rows_written
        + expected_m1_pairs
    )

    # --------------------------------------------------------
    # Global publication metadata. No duplicate classification
    # is performed here; that remains Stage 10.
    # --------------------------------------------------------

    provenance = {
        **stage5.provenance,
        "brain_prefilter_pipeline_git_commit": (
            stage7_commit
        ),
        "full_bri_pipeline_git_commit": (
            stage8_commit
        ),
    }

    summary = {
        "summary_schema_name": (
            "pdbclean_stage8_full_bri_global_summary"
        ),
        "summary_schema_version": "1.0",
        "snapshot": stage5.snapshot,
        "cleaning_protocol": (
            stage5.cleaning_protocol
        ),
        **provenance,
        "metric": "L_infinity",
        "distance_representation": (
            "exact_integer_milliangstrom"
        ),
        "brain_candidate_comparison_count": (
            candidate_rows_written
        ),
        "m1_bypass_chain_count": (
            m1_count
        ),
        "m1_direct_comparison_count": (
            expected_m1_pairs
        ),
        "total_full_bri_comparison_count": (
            total_comparisons
        ),
        "processing_error_count": 0,
        "pair_accounting_valid": True,
        "lemma_6_1_valid_for_all_brain_candidates": (
            True
        ),
        "classification_performed": False,
    }

    _write_json_atomic(
        summary,
        output_root / "global_summary.json",
    )

    success = {
        "success_schema_name": (
            "pdbclean_stage8_full_bri_success"
        ),
        "success_schema_version": "1.0",
        "snapshot": stage5.snapshot,
        "cleaning_protocol": (
            stage5.cleaning_protocol
        ),
        **provenance,
        "global_summary": (
            "global_summary.json"
        ),
        "candidate_comparisons": (
            "finalized/candidate_comparisons.parquet"
        ),
        "m1_comparisons": (
            "finalized/m1_comparisons.parquet"
        ),
    }

    # Completion marker strictly last.
    _write_json_atomic(
        success,
        success_path,
    )

    print()
    print("===== STAGE-8 FINAL ACCOUNTING =====")
    print(
        "Brain candidate comparisons:",
        f"{candidate_rows_written:,}",
    )
    print(
        "m=1 direct comparisons:",
        f"{expected_m1_pairs:,}",
    )
    print(
        "Total full-BRI comparisons:",
        f"{total_comparisons:,}",
    )
    print(
        "Lemma 6.1 violations: 0"
    )
    print()
    print(
        "STAGE-8 FULL-BRI PUBLICATION: PASS"
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
