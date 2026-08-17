"""Paper-faithful Stage 8: fast complete-BRI nearest-neighbour search.

Scientific cascade
-------------------
exact cleaned chain length
    -> Stage-7 Brain filtering
    -> compressed-cover-tree search on complete BRI
    -> exact complete-BRI L-infinity <= configured threshold

The configured threshold is the paper-derived near-duplicate threshold
(0.010 A in the canonical COMP702 protocol).

The existing exhaustive ``full_bri_compare`` publication is never
modified.  It is used only as an independent oracle for exact
pair-set and distance validation.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

from pdbclean.config import load_config
from pdbclean.full_bri_compare import bri_to_integer_mA
from pdbclean.full_bri_nn import (
    search_all_pairs_radius,
    search_brain_candidate_components,
)
from pdbclean.full_bri_production import (
    CANDIDATE_COMPARISON_SCHEMA,
    M1_COMPARISON_SCHEMA,
    _candidate_identity_indices,
    _identity,
    _read_json,
    _validate_commit,
    _write_json_atomic,
)
from pdbclean.length_buckets import (
    validate_stage5_brain_publication,
)


class FullBRINNProductionError(RuntimeError):
    """Raised when paper-faithful Stage-8 NN production cannot proceed."""


def _validate_threshold(config: dict) -> tuple[float, int]:
    value = float(
        config["duplicate_search"][
            "near_duplicate_threshold_angstrom"
        ]
    )

    scaled = value * 1000.0
    rounded = int(round(scaled))

    if abs(scaled - rounded) > 1.0e-9:
        raise FullBRINNProductionError(
            "Configured near-duplicate threshold is not exactly "
            "representable in integer milliangstroms"
        )

    if rounded <= 0:
        raise FullBRINNProductionError(
            "Configured near-duplicate threshold must be positive"
        )

    return value, rounded


def _canonical_pair_key(
    query_snapshot: str,
    query_pdb_id: str,
    query_model_id: int,
    query_chain: str,
    subject_snapshot: str,
    subject_pdb_id: str,
    subject_model_id: int,
    subject_chain: str,
) -> tuple[
    tuple[str, str, int, str],
    tuple[str, str, int, str],
]:
    left = (
        query_snapshot,
        query_pdb_id,
        int(query_model_id),
        query_chain,
    )
    right = (
        subject_snapshot,
        subject_pdb_id,
        int(subject_model_id),
        subject_chain,
    )

    if right < left:
        left, right = right, left

    return left, right


def _comparison_map(
    path: Path,
    *,
    threshold_mA: int,
) -> dict[
    tuple[
        tuple[str, str, int, str],
        tuple[str, str, int, str],
    ],
    int,
]:
    """Read <=threshold pairs from a Stage-8 comparison Parquet."""

    if not path.is_file():
        raise FullBRINNProductionError(
            f"Oracle comparison file is missing: {path}"
        )

    result: dict[
        tuple[
            tuple[str, str, int, str],
            tuple[str, str, int, str],
        ],
        int,
    ] = {}

    pf = pq.ParquetFile(path)

    columns = [
        "query_snapshot",
        "query_pdb_id",
        "query_model_id",
        "query_label_chain_id",
        "subject_snapshot",
        "subject_pdb_id",
        "subject_model_id",
        "subject_label_chain_id",
        "d_bri_mA",
    ]

    for batch in pf.iter_batches(
        columns=columns,
        batch_size=65_536,
    ):
        table = pa.Table.from_batches([batch])

        mask = pc.less_equal(
            table["d_bri_mA"],
            pa.scalar(
                threshold_mA,
                type=pa.int64(),
            ),
        )

        table = table.filter(mask)

        if table.num_rows == 0:
            continue

        rows = table.to_pylist()

        for row in rows:
            key = _canonical_pair_key(
                row["query_snapshot"],
                row["query_pdb_id"],
                row["query_model_id"],
                row["query_label_chain_id"],
                row["subject_snapshot"],
                row["subject_pdb_id"],
                row["subject_model_id"],
                row["subject_label_chain_id"],
            )

            distance = int(
                row["d_bri_mA"]
            )

            previous = result.get(key)

            if previous is not None:
                if previous != distance:
                    raise FullBRINNProductionError(
                        "Oracle contains one pair with conflicting "
                        "complete-BRI distances"
                    )

                raise FullBRINNProductionError(
                    "Oracle contains a duplicate unordered pair"
                )

            result[key] = distance

    return result


def _output_comparison_map(
    candidate_path: Path,
    m1_path: Path,
) -> dict[
    tuple[
        tuple[str, str, int, str],
        tuple[str, str, int, str],
    ],
    int,
]:
    """Read every pair from the new thresholded publication."""

    result: dict[
        tuple[
            tuple[str, str, int, str],
            tuple[str, str, int, str],
        ],
        int,
    ] = {}

    for path in (
        candidate_path,
        m1_path,
    ):
        if not path.is_file():
            raise FullBRINNProductionError(
                f"NN output file is missing: {path}"
            )

        pf = pq.ParquetFile(path)

        columns = [
            "query_snapshot",
            "query_pdb_id",
            "query_model_id",
            "query_label_chain_id",
            "subject_snapshot",
            "subject_pdb_id",
            "subject_model_id",
            "subject_label_chain_id",
            "d_bri_mA",
        ]

        for batch in pf.iter_batches(
            columns=columns,
            batch_size=65_536,
        ):
            for row in pa.Table.from_batches(
                [batch]
            ).to_pylist():
                key = _canonical_pair_key(
                    row["query_snapshot"],
                    row["query_pdb_id"],
                    row["query_model_id"],
                    row["query_label_chain_id"],
                    row["subject_snapshot"],
                    row["subject_pdb_id"],
                    row["subject_model_id"],
                    row["subject_label_chain_id"],
                )

                distance = int(
                    row["d_bri_mA"]
                )

                if key in result:
                    raise FullBRINNProductionError(
                        "NN publication contains a duplicate "
                        "unordered pair"
                    )

                result[key] = distance

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

    stage8_nn_commit = _validate_commit(
        args.pipeline_git_commit
    )

    repo = Path.cwd()

    loaded_config = load_config(
        args.config
    )
    config = loaded_config.data

    threshold_A, threshold_mA = (
        _validate_threshold(
            config
        )
    )

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
        raise FullBRINNProductionError(
            "Expected exactly one completed Stage-5 Brain publication"
        )

    stage5 = validate_stage5_brain_publication(
        completed_brain[0]
    )

    stage_root = (
        completed_brain[0].parent
    )

    # --------------------------------------------------------
    # Validate Stage 7.
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
        raise FullBRINNProductionError(
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
        raise FullBRINNProductionError(
            "Unexpected Stage-7 global-summary schema"
        )

    if (
        stage7_success["snapshot"]
        != stage5.snapshot
        or stage7_summary["snapshot"]
        != stage5.snapshot
    ):
        raise FullBRINNProductionError(
            "Stage-7 snapshot mismatch"
        )

    if (
        stage7_success["cleaning_protocol"]
        != stage5.cleaning_protocol
        or stage7_summary["cleaning_protocol"]
        != stage5.cleaning_protocol
    ):
        raise FullBRINNProductionError(
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
            raise FullBRINNProductionError(
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
        raise FullBRINNProductionError(
            "Stage-7 producer provenance mismatch"
        )

    if (
        stage7_summary[
            "processing_error_count"
        ]
        != 0
    ):
        raise FullBRINNProductionError(
            "Fast Stage 8 requires zero Stage-7 errors"
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

    candidate_count = (
        pq.read_metadata(
            candidate_path
        ).num_rows
    )
    bypass_count = (
        pq.read_metadata(
            bypass_path
        ).num_rows
    )

    if (
        candidate_count
        != stage7_summary[
            "candidate_pair_count"
        ]
    ):
        raise FullBRINNProductionError(
            "Stage-7 candidate-count mismatch"
        )

    if (
        bypass_count
        != stage7_summary[
            "m1_bypass_chain_count"
        ]
    ):
        raise FullBRINNProductionError(
            "Stage-7 m1-bypass count mismatch"
        )

    # --------------------------------------------------------
    # Determine exactly which length buckets need BRI loaded.
    # --------------------------------------------------------

    buckets = pq.read_table(
        bucket_path
    ).to_pylist()

    required_counts: dict[int, int] = {}

    for row in buckets:
        m = int(
            row["retained_residue_count"]
        )

        if m == 1:
            required_counts[1] = int(
                row["chain_count"]
            )
            continue

        if int(
            row["candidate_pair_count"]
        ) > 0:
            required_counts[m] = int(
                row["chain_count"]
            )

    if (
        required_counts.get(1, 0)
        != bypass_count
    ):
        raise FullBRINNProductionError(
            "Stage-7 m1 bucket accounting mismatch"
        )

    print(
        "Near-duplicate threshold:",
        f"{threshold_A:.3f} A "
        f"({threshold_mA} mA)",
    )
    print(
        "Required BRI length buckets:",
        f"{len(required_counts):,}",
    )

    # --------------------------------------------------------
    # Load canonical Stage-3 BRI.
    # --------------------------------------------------------

    bri_path = (
        stage_root
        / "bri/finalized/bri.parquet"
    )

    if not bri_path.is_file():
        raise FullBRINNProductionError(
            "Canonical Stage-3 BRI population is missing"
        )

    if (
        pq.read_metadata(
            bri_path
        ).num_rows
        != stage5.input_bri_chain_count
    ):
        raise FullBRINNProductionError(
            "Canonical Stage-3 BRI population count mismatch"
        )

    bri_by_m = {
        m: np.empty(
            (count, m, 9),
            dtype=np.int64,
        )
        for m, count in required_counts.items()
    }

    identity_to_index = {
        m: {}
        for m in required_counts
    }

    identity_by_index = {
        m: [None] * count
        for m, count in required_counts.items()
    }

    fill_count = {
        m: 0
        for m in required_counts
    }

    print(
        "Loading required canonical BRI chains..."
    )

    bri_pf = pq.ParquetFile(
        bri_path
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
            m = int(
                row[
                    "retained_residue_count"
                ]
            )

            if m not in required_counts:
                continue

            index = fill_count[m]

            if index >= required_counts[m]:
                raise FullBRINNProductionError(
                    f"Too many canonical BRI chains for m={m}"
                )

            identity = _identity(
                row
            )

            if (
                identity
                in identity_to_index[m]
            ):
                raise FullBRINNProductionError(
                    "Duplicate canonical BRI identity: "
                    f"{identity!r}"
                )

            integer_bri = (
                bri_to_integer_mA(
                    np.asarray(
                        row["bri"],
                        dtype=np.float64,
                    ),
                    name="stage8_nn_canonical_bri",
                )
            )

            if integer_bri.shape != (
                m,
                9,
            ):
                raise FullBRINNProductionError(
                    f"Canonical BRI shape mismatch for {identity!r}"
                )

            bri_by_m[m][
                index
            ] = integer_bri

            identity_to_index[m][
                identity
            ] = index

            identity_by_index[m][
                index
            ] = identity

            fill_count[m] += 1

    for m, expected in (
        required_counts.items()
    ):
        if fill_count[m] != expected:
            raise FullBRINNProductionError(
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
        raise FullBRINNProductionError(
            "Stage-7 m1 bypass identities do not exactly "
            "match canonical m1 BRI identities"
        )

    # --------------------------------------------------------
    # New paper-faithful publication namespace.
    # --------------------------------------------------------

    output_root = (
        stage_root / "full_bri_nn"
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
        / "candidate_near_duplicates.parquet"
    )
    candidate_tmp = (
        candidate_output.with_suffix(
            candidate_output.suffix
            + ".tmp"
        )
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

    seen_candidate_m: set[int] = set()

    candidate_hit_count = 0
    participating_chain_count = 0
    component_count = 0
    tree_count = 0
    non_candidate_radius_hit_count = 0

    try:
        for rg in range(
            candidate_pf.metadata.num_row_groups
        ):
            table = (
                candidate_pf.read_row_group(
                    rg
                )
            )

            if table.num_rows == 0:
                continue

            m_values = set(
                table[
                    "retained_residue_count"
                ].to_pylist()
            )

            if len(m_values) != 1:
                raise FullBRINNProductionError(
                    "Stage-7 candidate row group contains mixed m values"
                )

            m = int(
                next(iter(m_values))
            )

            if m < 2:
                raise FullBRINNProductionError(
                    "Stage-7 Brain candidate contains m < 2"
                )

            if m in seen_candidate_m:
                raise FullBRINNProductionError(
                    f"Stage-7 candidate bucket m={m} "
                    "appears in multiple row groups"
                )

            seen_candidate_m.add(m)

            if m not in bri_by_m:
                raise FullBRINNProductionError(
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
                raise FullBRINNProductionError(
                    f"Self-pair found in Stage-7 candidates for m={m}"
                )

            result = (
                search_brain_candidate_components(
                    bri_by_m[m],
                    q_indices,
                    s_indices,
                    radius_mA=threshold_mA,
                )
            )

            # A true complete-BRI hit <= threshold must also pass
            # the Brain threshold by dBrain <= dBRI.
            if (
                result.non_candidate_radius_hit_count
                != 0
            ):
                raise FullBRINNProductionError(
                    "Complete-BRI radius hit absent from the "
                    f"Stage-7 Brain edge relation for m={m}: "
                    f"{result.non_candidate_radius_hit_count}"
                )

            rows = result.row_indices

            selected = table.take(
                pa.array(
                    rows,
                    type=pa.int64(),
                )
            )

            d_bri_mA = (
                result.distances_mA
            )

            if (
                selected.num_rows
                != d_bri_mA.shape[0]
            ):
                raise FullBRINNProductionError(
                    "Fast-NN selected-row accounting mismatch"
                )

            d_num = np.asarray(
                selected[
                    "d_brain_numerator_max"
                ]
                .combine_chunks()
                .to_numpy(
                    zero_copy_only=False
                ),
                dtype=np.int64,
            )

            # Exact integer Lemma 6.1 validation on every emitted hit.
            if np.any(
                d_num
                > d_bri_mA
                * (m - 1)
            ):
                raise FullBRINNProductionError(
                    f"Lemma 6.1 violation in NN bucket m={m}"
                )

            d_brain = np.asarray(
                selected[
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
                raise FullBRINNProductionError(
                    "Stage-7 exact dBrain representation "
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
                    selected["query_snapshot"],
                    selected["query_pdb_id"],
                    selected["query_model_id"],
                    selected["query_label_chain_id"],

                    selected["subject_snapshot"],
                    selected["subject_pdb_id"],
                    selected["subject_model_id"],
                    selected["subject_label_chain_id"],

                    selected["retained_residue_count"],
                    selected["d_brain_numerator_max"],
                    selected["d_brain"],

                    pa.array(
                        d_bri_mA,
                        type=pa.int64(),
                    ),
                    pa.array(
                        d_bri,
                        type=pa.float64(),
                    ),
                ],
                schema=CANDIDATE_COMPARISON_SCHEMA,
            )

            if output.num_rows:
                writer.write_table(output)

            candidate_hit_count += (
                output.num_rows
            )
            participating_chain_count += (
                result.participating_chain_count
            )
            component_count += (
                result.component_count
            )
            tree_count += (
                result.compressed_tree_count
            )
            non_candidate_radius_hit_count += (
                result.non_candidate_radius_hit_count
            )

            print(
                f"m={m:5d} "
                f"brain_pairs={table.num_rows:10,d} "
                f"hits={output.num_rows:9,d} "
                f"components={result.component_count:6,d}"
            )

    finally:
        writer.close()

    candidate_tmp.replace(
        candidate_output
    )

    # --------------------------------------------------------
    # m=1: Brain undefined, so search the full m=1 BRI bucket.
    # --------------------------------------------------------

    m1_result = search_all_pairs_radius(
        bri_by_m[1],
        radius_mA=threshold_mA,
    )

    q_identity = [
        identity_by_index[1][
            int(index)
        ]
        for index in (
            m1_result.query_indices
        )
    ]
    s_identity = [
        identity_by_index[1][
            int(index)
        ]
        for index in (
            m1_result.subject_indices
        )
    ]

    m1_d_bri_mA = (
        m1_result.distances_mA
    )
    m1_d_bri = (
        m1_d_bri_mA.astype(
            np.float64
        )
        / 1000.0
    )

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
                    m1_result.query_indices.shape[0],
                    dtype=np.int64,
                ),
                type=pa.int64(),
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
        / "m1_near_duplicates.parquet"
    )
    m1_tmp = (
        m1_output.with_suffix(
            m1_output.suffix
            + ".tmp"
        )
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

    m1_tmp.replace(
        m1_output
    )

    # --------------------------------------------------------
    # Independent validation against the existing exhaustive oracle.
    # --------------------------------------------------------

    oracle_root = (
        stage_root / "full_bri_compare"
    )

    oracle_success = _read_json(
        oracle_root / "_SUCCESS"
    )

    if (
        oracle_success.get(
            "success_schema_name"
        )
        != "pdbclean_stage8_full_bri_success"
        or oracle_success.get(
            "success_schema_version"
        )
        != "1.0"
    ):
        raise FullBRINNProductionError(
            "Unexpected exhaustive Stage-8 oracle schema"
        )

    if (
        oracle_success.get("snapshot")
        != stage5.snapshot
        or oracle_success.get(
            "cleaning_protocol"
        )
        != stage5.cleaning_protocol
    ):
        raise FullBRINNProductionError(
            "Exhaustive Stage-8 oracle population mismatch"
        )

    if (
        oracle_success.get(
            "brain_prefilter_pipeline_git_commit"
        )
        != stage7_commit
    ):
        raise FullBRINNProductionError(
            "Exhaustive oracle Stage-7 provenance mismatch"
        )

    oracle_candidate_path = (
        oracle_root
        / oracle_success[
            "candidate_comparisons"
        ]
    )
    oracle_m1_path = (
        oracle_root
        / oracle_success[
            "m1_comparisons"
        ]
    )

    print()
    print(
        "Validating fast-NN pair set against exhaustive oracle..."
    )

    oracle_pairs = _comparison_map(
        oracle_candidate_path,
        threshold_mA=threshold_mA,
    )

    oracle_m1_pairs = _comparison_map(
        oracle_m1_path,
        threshold_mA=threshold_mA,
    )

    overlap = (
        set(oracle_pairs)
        & set(oracle_m1_pairs)
    )

    if overlap:
        raise FullBRINNProductionError(
            "Oracle candidate and m1 populations overlap"
        )

    oracle_pairs.update(
        oracle_m1_pairs
    )

    nn_pairs = _output_comparison_map(
        candidate_output,
        m1_output,
    )

    if (
        set(nn_pairs)
        != set(oracle_pairs)
    ):
        missing = (
            set(oracle_pairs)
            - set(nn_pairs)
        )
        extra = (
            set(nn_pairs)
            - set(oracle_pairs)
        )

        raise FullBRINNProductionError(
            "Fast-NN pair-set mismatch against exhaustive oracle: "
            f"missing={len(missing):,}, "
            f"extra={len(extra):,}"
        )

    mismatched_distances = sum(
        nn_pairs[key]
        != oracle_pairs[key]
        for key in nn_pairs
    )

    if mismatched_distances:
        raise FullBRINNProductionError(
            "Fast-NN complete-BRI distance mismatch against "
            f"exhaustive oracle: {mismatched_distances:,}"
        )

    total_hit_count = (
        candidate_hit_count
        + m1_table.num_rows
    )

    if (
        total_hit_count
        != len(oracle_pairs)
        or total_hit_count
        != len(nn_pairs)
    ):
        raise FullBRINNProductionError(
            "Final near-duplicate pair accounting mismatch"
        )

    print(
        "Exhaustive-oracle pair-set equality: PASS"
    )
    print(
        "Exhaustive-oracle distance equality: PASS"
    )

    # --------------------------------------------------------
    # Publication metadata.
    # --------------------------------------------------------

    provenance = {
        **stage5.provenance,
        "brain_prefilter_pipeline_git_commit": (
            stage7_commit
        ),
        "full_bri_nn_pipeline_git_commit": (
            stage8_nn_commit
        ),
        "config_sha256": (
            loaded_config.sha256
        ),
    }

    summary = {
        "summary_schema_name": (
            "pdbclean_stage8_full_bri_nn_global_summary"
        ),
        "summary_schema_version": "1.0",
        "snapshot": stage5.snapshot,
        "cleaning_protocol": (
            stage5.cleaning_protocol
        ),
        **provenance,

        "scientific_sequence": (
            "exact_length_then_brain_then_complete_bri_fast_nn"
        ),
        "search_engine": (
            "elkin_kurlin_compressed_cover_tree"
        ),
        "metric": "L_infinity",
        "distance_representation": (
            "exact_integer_milliangstrom"
        ),
        "near_duplicate_threshold_angstrom": (
            threshold_A
        ),
        "near_duplicate_threshold_mA": (
            threshold_mA
        ),
        "threshold_operator": (
            "less_than_or_equal"
        ),

        "brain_candidate_pair_count": (
            candidate_count
        ),
        "brain_participating_chain_count": (
            participating_chain_count
        ),
        "brain_connected_component_count": (
            component_count
        ),
        "compressed_tree_count": (
            tree_count
        ),
        "candidate_near_duplicate_count": (
            candidate_hit_count
        ),

        "m1_bypass_chain_count": (
            required_counts[1]
        ),
        "m1_unique_bri_point_count": (
            m1_result.unique_point_count
        ),
        "m1_near_duplicate_count": (
            m1_table.num_rows
        ),

        "total_near_duplicate_count": (
            total_hit_count
        ),
        "non_candidate_radius_hit_count": (
            non_candidate_radius_hit_count
        ),

        "exhaustive_oracle_pair_count": (
            len(oracle_pairs)
        ),
        "exhaustive_oracle_pair_set_equal": (
            True
        ),
        "exhaustive_oracle_distances_equal": (
            True
        ),
        "lemma_6_1_valid_for_all_emitted_brain_pairs": (
            True
        ),

        "processing_error_count": 0,
        "pair_accounting_valid": True,
        # Stage 8 applies the geometric near-duplicate threshold.
        # Zero/nonzero duplicate classification remains a separate
        # downstream analysis stage.
        "classification_performed": False,
    }

    _write_json_atomic(
        summary,
        output_root / "global_summary.json",
    )

    success = {
        "success_schema_name": (
            "pdbclean_stage8_full_bri_nn_success"
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
        "candidate_near_duplicates": (
            "finalized/candidate_near_duplicates.parquet"
        ),
        "m1_near_duplicates": (
            "finalized/m1_near_duplicates.parquet"
        ),
    }

    # Completion marker strictly last.
    _write_json_atomic(
        success,
        success_path,
    )

    print()
    print(
        "===== PAPER-FAITHFUL STAGE-8 FINAL ACCOUNTING ====="
    )
    print(
        "Stage-7 Brain candidates:",
        f"{candidate_count:,}",
    )
    print(
        "Brain connected components:",
        f"{component_count:,}",
    )
    print(
        "Compressed cover trees:",
        f"{tree_count:,}",
    )
    print(
        "Candidate near-duplicates:",
        f"{candidate_hit_count:,}",
    )
    print(
        "m=1 near-duplicates:",
        f"{m1_table.num_rows:,}",
    )
    print(
        "Total dBRI <= threshold:",
        f"{total_hit_count:,}",
    )
    print(
        "Non-candidate radius hits:",
        f"{non_candidate_radius_hit_count:,}",
    )
    print(
        "STAGE-8 FAST COMPLETE-BRI NN PUBLICATION: PASS"
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
