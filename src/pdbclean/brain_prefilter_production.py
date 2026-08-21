"""Stage 7 production publication for the lossless Brain prefilter."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

from pdbclean.brain_prefilter import (
    BRAIN_PREFILTER_TAU_ANGSTROM,
    BRAIN_PREFILTER_TAU_MA,
    brain_candidate_pairs,
    brain_integer_numerators,
)
from pdbclean.config import load_config
from pdbclean.defaults import (
    brain_filter_threshold_milliangstrom,
    require_implemented_precision,
)
from pdbclean.length_buckets import (
    validate_stage5_brain_publication,
)


CANDIDATE_SCHEMA = pa.schema(
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
    ],
    metadata={
        b"schema_name": b"pdbclean_stage7_brain_candidates",
        b"schema_version": b"1.0",
        b"metric": b"9D L-infinity Brain distance",
        b"threshold_angstrom": b"0.010",
        b"search": b"scipy cKDTree p=inf eps=0 plus exact integer post-filter",
    },
)


BYPASS_SCHEMA = pa.schema(
    [
        pa.field("snapshot", pa.string(), nullable=False),
        pa.field("pdb_id", pa.string(), nullable=False),
        pa.field("model_id", pa.int64(), nullable=False),
        pa.field("label_chain_id", pa.string(), nullable=False),
        pa.field(
            "retained_residue_count",
            pa.int64(),
            nullable=False,
        ),
        pa.field(
            "undefined_reason",
            pa.string(),
            nullable=False,
        ),
    ],
    metadata={
        b"schema_name": b"pdbclean_stage7_m1_bypass",
        b"schema_version": b"1.0",
    },
)


BUCKET_SCHEMA = pa.schema(
    [
        pa.field(
            "retained_residue_count",
            pa.int64(),
            nullable=False,
        ),
        pa.field(
            "chain_count",
            pa.int64(),
            nullable=False,
        ),
        pa.field(
            "possible_pair_count",
            pa.int64(),
            nullable=False,
        ),
        pa.field(
            "candidate_pair_count",
            pa.int64(),
            nullable=False,
        ),
        pa.field(
            "search_mode",
            pa.string(),
            nullable=False,
        ),
    ],
    metadata={
        b"schema_name": b"pdbclean_stage7_brain_bucket_summary",
        b"schema_version": b"1.0",
    },
)


def _format_threshold_angstrom(threshold_mA: int) -> str:
    """Render an integer-milliangstrom threshold as an angstrom string.

    ``10`` renders as ``"0.010"``, exactly reproducing the literal that the
    frozen Stage-7 publication carries in its Parquet schema metadata.
    """

    return f"{threshold_mA / 1000.0:.3f}"


def _candidate_schema(threshold_mA: int) -> pa.Schema:
    """Return the Stage-7 candidate schema for a configured threshold.

    The schema's declared threshold must describe the threshold the data was
    actually produced with, so it is derived from the resolved configuration
    rather than from a literal.  At the validated default (10 mA) this
    reproduces the frozen metadata byte for byte.
    """

    metadata = {
        key.decode("utf-8"): value.decode("utf-8")
        for key, value in (CANDIDATE_SCHEMA.metadata or {}).items()
    }

    metadata["threshold_angstrom"] = _format_threshold_angstrom(threshold_mA)

    return CANDIDATE_SCHEMA.with_metadata(metadata)


def _write_json_atomic(payload, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)

    tmp = path.with_suffix(path.suffix + ".tmp")

    tmp.write_text(
        json.dumps(
            payload,
            sort_keys=True,
            indent=2,
            ensure_ascii=True,
        )
        + "\n",
        encoding="utf-8",
    )

    tmp.replace(path)


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

    if (
        len(args.pipeline_git_commit) != 40
        or any(
            c not in "0123456789abcdef"
            for c in args.pipeline_git_commit
        )
    ):
        raise RuntimeError(
            "pipeline commit must be lowercase 40-character SHA"
        )

    repo = Path.cwd()

    config = load_config(args.config).data
    protocol = config["release"]["protocol_version"]

    # The Brain filtering threshold comes from the resolved configuration.
    # A legacy protocol configuration has no `brain_filter` section, so the
    # validated default is used, which is exactly what the frozen Stage-7
    # publication was produced with.
    # The Brain prefilter works in exact integer milliangstrom sums, so it
    # is only valid on the implemented precision grid.
    require_implemented_precision(config, stage="stage7_brain_prefilter")

    brain_tau_mA = brain_filter_threshold_milliangstrom(
        {
            "brain_filter": config.get(
                "brain_filter",
                {"threshold_angstrom": BRAIN_PREFILTER_TAU_ANGSTROM},
            )
        }
    )

    brain_tau_angstrom = brain_tau_mA / 1000.0

    storage_root = Path(
        config["storage"]["output_root"]
    )

    if not storage_root.is_absolute():
        storage_root = repo / storage_root

    completed = sorted(
        path.parent
        for path in storage_root.glob(
            f"*/{protocol}/brain/_SUCCESS"
        )
    )

    if len(completed) != 1:
        raise RuntimeError(
            "Expected exactly one completed Stage-5 publication"
        )

    upstream = validate_stage5_brain_publication(
        completed[0]
    )

    stage_root = completed[0].parent

    bucket_path = (
        stage_root
        / "length_buckets/finalized/bucket_index.parquet"
    )

    if not bucket_path.is_file():
        raise RuntimeError(
            "Canonical Stage-6 bucket index is missing"
        )

    stage6 = {
        row["retained_residue_count"]: row
        for row in pq.read_table(
            bucket_path
        ).to_pylist()
    }

    output_root = (
        stage_root / "brain_prefilter"
    )
    finalized = output_root / "finalized"

    finalized.mkdir(
        parents=True,
        exist_ok=True,
    )

    success_path = output_root / "_SUCCESS"

    if success_path.exists():
        success_path.unlink()

    brain = pq.read_table(
        upstream.brain_path,
        columns=[
            "snapshot",
            "pdb_id",
            "model_id",
            "label_chain_id",
            "retained_residue_count",
            "brain",
        ],
    )

    m_values = np.asarray(
        brain[
            "retained_residue_count"
        ].combine_chunks().to_numpy(
            zero_copy_only=False
        ),
        dtype=np.int64,
    )

    brain_column = (
        brain["brain"].combine_chunks()
    )

    brain_values = np.asarray(
        brain_column.values.to_numpy(
            zero_copy_only=False
        ),
        dtype=np.float64,
    ).reshape(-1, 9)

    if brain_values.shape != (
        upstream.brain_chain_count,
        9,
    ):
        raise RuntimeError(
            "Unexpected finalized Brain matrix shape"
        )

    # Stable exact-m grouping without changing canonical
    # within-bucket chain order.
    order = np.argsort(
        m_values,
        kind="stable",
    )

    sorted_m = m_values[order]

    unique_m, starts, counts = np.unique(
        sorted_m,
        return_index=True,
        return_counts=True,
    )

    candidate_tmp = (
        finalized / "candidates.parquet.tmp"
    )
    candidate_path = (
        finalized / "candidates.parquet"
    )

    if candidate_tmp.exists():
        candidate_tmp.unlink()

    candidate_schema = _candidate_schema(brain_tau_mA)

    writer = pq.ParquetWriter(
        candidate_tmp,
        candidate_schema,
        compression="zstd",
        version="2.6",
        use_dictionary=True,
    )

    bucket_rows = []

    total_candidates = 0
    total_defined = 0

    try:
        for m_raw, start_raw, count_raw in zip(
            unique_m,
            starts,
            counts,
            strict=True,
        ):
            m = int(m_raw)
            start = int(start_raw)
            count = int(count_raw)

            if m < 2:
                raise RuntimeError(
                    "Brain-defined population contains m < 2"
                )

            indices = order[
                start:start + count
            ]

            vectors = brain_values[
                indices
            ]

            result = brain_candidate_pairs(
                vectors,
                m=m,
                tau_mA=brain_tau_mA,
            )

            pairs = result.pairs
            pair_count = result.pair_count

            total_defined += count
            total_candidates += pair_count

            expected = stage6.get(m)

            if expected is None:
                raise RuntimeError(
                    f"m={m} absent from Stage-6 index"
                )

            if (
                expected["brain_defined_count"]
                != count
            ):
                raise RuntimeError(
                    f"Stage-6 population mismatch for m={m}"
                )

            bucket_rows.append(
                {
                    "retained_residue_count": m,
                    "chain_count": count,
                    "possible_pair_count": (
                        count * (count - 1) // 2
                    ),
                    "candidate_pair_count": pair_count,
                    "search_mode": "brain_ckdtree",
                }
            )

            if pair_count == 0:
                continue

            numerators = (
                brain_integer_numerators(
                    vectors,
                    m=m,
                )
            )

            d_num = np.max(
                np.abs(
                    numerators[
                        pairs[:, 0]
                    ]
                    - numerators[
                        pairs[:, 1]
                    ]
                ),
                axis=1,
            ).astype(
                np.int64,
                copy=False,
            )

            exact_radius = (
                brain_tau_mA
                * (m - 1)
            )

            if np.any(
                d_num > exact_radius
            ):
                raise RuntimeError(
                    f"Stage-7 exact threshold violation for m={m}"
                )

            d_brain = (
                d_num.astype(
                    np.float64
                )
                / (
                    1000.0
                    * float(m - 1)
                )
            )

            if np.any(
                d_brain
                > brain_tau_angstrom
            ):
                raise RuntimeError(
                    f"Stage-7 floating distance violation for m={m}"
                )

            query_indices = indices[
                pairs[:, 0]
            ]
            subject_indices = indices[
                pairs[:, 1]
            ]

            qtake = pa.array(
                query_indices,
                type=pa.int64(),
            )
            stake = pa.array(
                subject_indices,
                type=pa.int64(),
            )

            table = pa.Table.from_arrays(
                [
                    pc.take(brain["snapshot"], qtake),
                    pc.take(brain["pdb_id"], qtake),
                    pc.take(brain["model_id"], qtake),
                    pc.take(brain["label_chain_id"], qtake),

                    pc.take(brain["snapshot"], stake),
                    pc.take(brain["pdb_id"], stake),
                    pc.take(brain["model_id"], stake),
                    pc.take(brain["label_chain_id"], stake),

                    pa.array(
                        np.full(
                            pair_count,
                            m,
                            dtype=np.int64,
                        )
                    ),
                    pa.array(d_num),
                    pa.array(d_brain),
                ],
                schema=CANDIDATE_SCHEMA,
            )

            writer.write_table(
                table
            )

    finally:
        writer.close()

    if total_defined != upstream.brain_chain_count:
        raise RuntimeError(
            "Stage-7 Brain-defined accounting failed"
        )

    if (
        pq.read_metadata(
            candidate_tmp
        ).num_rows
        != total_candidates
    ):
        raise RuntimeError(
            "Stage-7 candidate Parquet count mismatch"
        )

    candidate_tmp.replace(
        candidate_path
    )

    # m=1 bypass is explicit and does not pass through Brain.
    undefined = pq.read_table(
        upstream.undefined_path,
        columns=[
            "snapshot",
            "pdb_id",
            "model_id",
            "label_chain_id",
            "retained_residue_count",
            "undefined_reason",
        ],
    )

    if undefined.num_rows != upstream.undefined_chain_count:
        raise RuntimeError(
            "Stage-7 m1 bypass count mismatch"
        )

    if any(
        m != 1
        for m in undefined[
            "retained_residue_count"
        ].to_pylist()
    ):
        raise RuntimeError(
            "Stage-7 bypass contains m != 1"
        )

    bypass = pa.Table.from_arrays(
        [
            undefined["snapshot"],
            undefined["pdb_id"],
            undefined["model_id"],
            undefined["label_chain_id"],
            undefined["retained_residue_count"],
            undefined["undefined_reason"],
        ],
        schema=BYPASS_SCHEMA,
    )

    bypass_tmp = (
        finalized / "m1_bypass.parquet.tmp"
    )
    bypass_path = (
        finalized / "m1_bypass.parquet"
    )

    pq.write_table(
        bypass,
        bypass_tmp,
        compression="zstd",
        version="2.6",
    )
    bypass_tmp.replace(
        bypass_path
    )

    bucket_rows.append(
        {
            "retained_residue_count": 1,
            "chain_count": bypass.num_rows,
            "possible_pair_count": (
                bypass.num_rows
                * (bypass.num_rows - 1)
                // 2
            ),
            "candidate_pair_count": 0,
            "search_mode": "full_bri_bypass",
        }
    )

    bucket_rows.sort(
        key=lambda row:
        row["retained_residue_count"]
    )

    bucket_table = pa.Table.from_pylist(
        bucket_rows,
        schema=BUCKET_SCHEMA,
    )

    bucket_tmp = (
        finalized
        / "bucket_summary.parquet.tmp"
    )
    bucket_final = (
        finalized
        / "bucket_summary.parquet"
    )

    pq.write_table(
        bucket_table,
        bucket_tmp,
        compression="zstd",
        version="2.6",
    )
    bucket_tmp.replace(
        bucket_final
    )

    if bucket_table.num_rows != len(stage6):
        raise RuntimeError(
            "Stage-7 bucket-count mismatch"
        )

    if (
        total_defined + bypass.num_rows
        != upstream.input_bri_chain_count
    ):
        raise RuntimeError(
            "Stage-7 population accounting failed"
        )

    summary = {
        "summary_schema_name": (
            "pdbclean_stage7_brain_prefilter_global_summary"
        ),
        "summary_schema_version": "1.0",
        "snapshot": upstream.snapshot,
        "cleaning_protocol": upstream.cleaning_protocol,
        **upstream.provenance,
        "brain_prefilter_pipeline_git_commit": (
            args.pipeline_git_commit
        ),
        "brain_threshold_angstrom": (
            brain_tau_angstrom
        ),
        "brain_threshold_mA": (
            brain_tau_mA
        ),
        "brain_threshold_operator": (
            "less_than_or_equal"
        ),
        "metric": "L_infinity",
        "search_implementation": (
            "scipy_cKDTree_p_inf_eps_0_exact_integer_postfilter"
        ),
        "input_chain_count": (
            upstream.input_bri_chain_count
        ),
        "brain_defined_chain_count": (
            total_defined
        ),
        "m1_bypass_chain_count": (
            bypass.num_rows
        ),
        "length_bucket_count": (
            bucket_table.num_rows
        ),
        "brain_searched_bucket_count": (
            len(unique_m)
        ),
        "candidate_pair_count": (
            total_candidates
        ),
        "processing_error_count": 0,
        "population_accounting_valid": True,
    }

    _write_json_atomic(
        summary,
        output_root / "global_summary.json",
    )

    success = {
        "success_schema_name": (
            "pdbclean_stage7_brain_prefilter_success"
        ),
        "success_schema_version": "1.0",
        "snapshot": upstream.snapshot,
        "cleaning_protocol": upstream.cleaning_protocol,
        **upstream.provenance,
        "brain_prefilter_pipeline_git_commit": (
            args.pipeline_git_commit
        ),
        "global_summary": "global_summary.json",
        "candidate_pairs": (
            "finalized/candidates.parquet"
        ),
        "m1_bypass": (
            "finalized/m1_bypass.parquet"
        ),
        "bucket_summary": (
            "finalized/bucket_summary.parquet"
        ),
    }

    # Completion marker strictly last.
    _write_json_atomic(
        success,
        success_path,
    )

    print("Input chains:", f"{upstream.input_bri_chain_count:,}")
    print("Brain-defined:", f"{total_defined:,}")
    print("m=1 bypass:", f"{bypass.num_rows:,}")
    print("Length buckets:", f"{bucket_table.num_rows:,}")
    print("Candidate pairs:", f"{total_candidates:,}")
    print()
    print("STAGE-7 BRAIN PREFILTER PUBLICATION: PASS")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
