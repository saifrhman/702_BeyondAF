#!/usr/bin/env python3

import argparse
import hashlib
import json
import shutil
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def canonical_key_from_source(d, i):
    return (
        str(d["snapshot"][i]),
        str(d["pdb_id"][i]).lower(),
        int(d["model_id"][i]),
        str(d["label_chain_id"][i]),
    )


def canonical_key_from_mapping(r):
    return (
        str(r["snapshot"]),
        str(r["pdb_id"]).lower(),
        int(r["model_id"]),
        str(r["label_chain_id"]),
    )


def representative_key(r):
    return (
        str(r["representative_snapshot"]),
        str(r["representative_pdb_id"]).lower(),
        int(r["representative_model_id"]),
        str(r["representative_label_chain_id"]),
    )


def parse_args():
    p = argparse.ArgumentParser()

    p.add_argument("--protocol-root", type=Path, required=True)
    p.add_argument("--policy-config", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)

    return p.parse_args()


def main():
    args = parse_args()

    root = args.protocol_root

    canonical = (
        root
        / "geometric_validation"
        / "finalized"
        / "eligible.parquet"
    )

    representative_root = (
        root
        / "stage14_representative_selection_v1"
    )

    mapping_path = (
        representative_root
        / "representative_mapping.parquet"
    )

    representatives_path = (
        representative_root
        / "representatives.parquet"
    )

    representative_summary = (
        representative_root
        / "global_summary.json"
    )

    representative_success = (
        representative_root
        / "_SUCCESS"
    )

    graph_root = (
        root
        / "stage14_geometric_graph"
    )

    component_summary = (
        graph_root
        / "component_summary.parquet"
    )

    graph_summary = (
        graph_root
        / "global_summary.json"
    )

    graph_success = (
        graph_root
        / "_SUCCESS"
    )

    full_bri_root = (
        root
        / "full_bri_nn"
    )

    edges_mge2 = (
        full_bri_root
        / "finalized"
        / "candidate_near_duplicates.parquet"
    )

    edges_m1 = (
        full_bri_root
        / "finalized"
        / "m1_near_duplicates.parquet"
    )

    stage8_summary = (
        full_bri_root
        / "global_summary.json"
    )

    stage8_success = (
        full_bri_root
        / "_SUCCESS"
    )

    source_manifest = (
        root.parent
        / "bronze"
        / "source_manifest.parquet"
    )

    stage13_success = (
        root
        / "acta_detailed_review_v2"
        / "_SUCCESS"
    )

    required = [
        canonical,
        mapping_path,
        representatives_path,
        representative_summary,
        representative_success,
        component_summary,
        graph_summary,
        graph_success,
        edges_mge2,
        edges_m1,
        stage8_summary,
        stage8_success,
        source_manifest,
        stage13_success,
        args.policy_config,
    ]

    for path in required:
        assert path.is_file(), path

    assert not args.output_dir.exists(), (
        f"Release already exists: {args.output_dir}"
    )

    # ========================================================
    # Validate frozen upstream publications.
    # ========================================================

    with representative_success.open(
        "r",
        encoding="utf-8",
    ) as f:
        rep_success = json.load(f)

    assert rep_success["status"] == "PASS"
    assert rep_success["input_chains"] == 578_524
    assert rep_success["retained_chains"] == 499_770
    assert rep_success["removed_chains"] == 78_754
    assert rep_success["direct_edge_safety"] is True
    assert rep_success["m1_deduplication_performed"] is False

    with representative_summary.open(
        "r",
        encoding="utf-8",
    ) as f:
        rep_summary = json.load(f)

    assert rep_summary["canonical_input_chain_count"] == 578_524
    assert rep_summary["final_retained_chain_count"] == 499_770
    assert rep_summary["final_removed_chain_count"] == 78_754

    assert (
        rep_summary[
            "every_removed_chain_has_direct_representative_edge"
        ]
        is True
    )

    assert (
        rep_summary["connectedness_treated_as_equivalence"]
        is False
    )

    assert (
        rep_summary["m1_deduplication_performed"]
        is False
    )

    assert (
        rep_summary["old_snapshot_comparison_used"]
        is False
    )

    # ========================================================
    # Load representative mapping.
    # ========================================================

    mapping = pq.read_table(
        mapping_path
    ).to_pylist()

    assert len(mapping) == 99_854

    mapping_by_key = {
        canonical_key_from_mapping(r): r
        for r in mapping
    }

    assert len(mapping_by_key) == 99_854

    removed_by_key = {
        canonical_key_from_mapping(r): r
        for r in mapping
        if r["action"] == "remove"
    }

    representative_keys = {
        canonical_key_from_mapping(r)
        for r in mapping
        if r["action"] == "retain_representative"
    }

    assert len(removed_by_key) == 78_754
    assert len(representative_keys) == 21_100

    assert not (
        set(removed_by_key)
        & representative_keys
    )

    # Every removal points to an actual retained representative.
    for key, row in removed_by_key.items():
        rep = representative_key(row)

        assert rep in representative_keys
        assert 0 <= int(row["direct_d_bri_mA"]) <= 10

    # ========================================================
    # Create release layout.
    # ========================================================

    data_dir = args.output_dir / "data"
    audit_dir = args.output_dir / "audit"
    provenance_dir = args.output_dir / "provenance"

    data_dir.mkdir(parents=True)
    audit_dir.mkdir(parents=True)
    provenance_dir.mkdir(parents=True)

    retained_output = (
        data_dir
        / "retained_chains.parquet"
    )

    removed_output = (
        audit_dir
        / "removed_chain_audit.parquet"
    )

    # ========================================================
    # Stream canonical 578,524-chain source.
    #
    # retained_chains.parquet preserves the exact canonical
    # post-geometry schema. No derived BRI column is injected.
    #
    # removed_chain_audit.parquet contains the corresponding
    # canonical chain evidence plus the direct representative
    # justification.
    # ========================================================

    source_pf = pq.ParquetFile(canonical)

    assert source_pf.metadata.num_rows == 578_524

    source_schema = source_pf.schema_arrow

    extra_removed_schema = pa.schema([
        pa.field(
            "dedup_component_id",
            pa.string(),
            nullable=False,
        ),
        pa.field(
            "dedup_component_is_clique",
            pa.bool_(),
            nullable=False,
        ),
        pa.field(
            "dedup_representative_snapshot",
            pa.string(),
            nullable=False,
        ),
        pa.field(
            "dedup_representative_pdb_id",
            pa.string(),
            nullable=False,
        ),
        pa.field(
            "dedup_representative_model_id",
            pa.int64(),
            nullable=False,
        ),
        pa.field(
            "dedup_representative_label_chain_id",
            pa.string(),
            nullable=False,
        ),
        pa.field(
            "dedup_direct_d_bri_mA",
            pa.int64(),
            nullable=False,
        ),
        pa.field(
            "dedup_policy_version",
            pa.string(),
            nullable=False,
        ),
    ])

    removed_schema = pa.schema(
        list(source_schema)
        + list(extra_removed_schema)
    )

    retained_writer = pq.ParquetWriter(
        retained_output,
        source_schema,
        compression="zstd",
    )

    removed_writer = pq.ParquetWriter(
        removed_output,
        removed_schema,
        compression="zstd",
    )

    source_seen = 0
    retained_count = 0
    removed_count = 0
    m1_count = 0

    removed_found = set()
    representative_found = set()

    key_columns = [
        "snapshot",
        "pdb_id",
        "model_id",
        "label_chain_id",
        "retained_residue_count",
    ]

    try:
        for batch in source_pf.iter_batches(
            batch_size=32768
        ):
            d = batch.select(
                [
                    batch.schema.get_field_index(c)
                    for c in key_columns
                ]
            ).to_pydict()

            retain_mask = []
            remove_mask = []

            removed_mapping_rows = []

            for i in range(batch.num_rows):
                key = canonical_key_from_source(
                    d,
                    i,
                )

                m = int(
                    d["retained_residue_count"][i]
                )

                if m == 1:
                    m1_count += 1

                if key in removed_by_key:
                    retain_mask.append(False)
                    remove_mask.append(True)

                    row = removed_by_key[key]
                    removed_mapping_rows.append(row)
                    removed_found.add(key)

                else:
                    retain_mask.append(True)
                    remove_mask.append(False)

                    if key in representative_keys:
                        representative_found.add(key)

            retained_batch = batch.filter(
                pa.array(retain_mask)
            )

            if retained_batch.num_rows:
                retained_writer.write_batch(
                    retained_batch
                )

                retained_count += (
                    retained_batch.num_rows
                )

            removed_batch = batch.filter(
                pa.array(remove_mask)
            )

            if removed_batch.num_rows:
                assert (
                    removed_batch.num_rows
                    == len(removed_mapping_rows)
                )

                table = pa.Table.from_batches(
                    [removed_batch]
                )

                table = table.append_column(
                    "dedup_component_id",
                    pa.array(
                        [
                            str(r["component_id"])
                            for r
                            in removed_mapping_rows
                        ],
                        type=pa.string(),
                    ),
                )

                table = table.append_column(
                    "dedup_component_is_clique",
                    pa.array(
                        [
                            bool(
                                r[
                                    "component_is_clique"
                                ]
                            )
                            for r
                            in removed_mapping_rows
                        ],
                        type=pa.bool_(),
                    ),
                )

                table = table.append_column(
                    "dedup_representative_snapshot",
                    pa.array(
                        [
                            str(
                                r[
                                    "representative_snapshot"
                                ]
                            )
                            for r
                            in removed_mapping_rows
                        ],
                        type=pa.string(),
                    ),
                )

                table = table.append_column(
                    "dedup_representative_pdb_id",
                    pa.array(
                        [
                            str(
                                r[
                                    "representative_pdb_id"
                                ]
                            ).lower()
                            for r
                            in removed_mapping_rows
                        ],
                        type=pa.string(),
                    ),
                )

                table = table.append_column(
                    "dedup_representative_model_id",
                    pa.array(
                        [
                            int(
                                r[
                                    "representative_model_id"
                                ]
                            )
                            for r
                            in removed_mapping_rows
                        ],
                        type=pa.int64(),
                    ),
                )

                table = table.append_column(
                    "dedup_representative_label_chain_id",
                    pa.array(
                        [
                            str(
                                r[
                                    "representative_label_chain_id"
                                ]
                            )
                            for r
                            in removed_mapping_rows
                        ],
                        type=pa.string(),
                    ),
                )

                table = table.append_column(
                    "dedup_direct_d_bri_mA",
                    pa.array(
                        [
                            int(
                                r[
                                    "direct_d_bri_mA"
                                ]
                            )
                            for r
                            in removed_mapping_rows
                        ],
                        type=pa.int64(),
                    ),
                )

                table = table.append_column(
                    "dedup_policy_version",
                    pa.array(
                        [
                            str(
                                r[
                                    "policy_version"
                                ]
                            )
                            for r
                            in removed_mapping_rows
                        ],
                        type=pa.string(),
                    ),
                )

                # append_column() infers appended fields as nullable even
                # though the release audit schema deliberately declares these
                # provenance fields non-nullable.  Validate names/types first,
                # then reconstruct the table with the canonical target schema.
                assert table.schema.names == removed_schema.names

                for actual, expected in zip(
                    table.schema,
                    removed_schema,
                ):
                    assert actual.type == expected.type, (
                        actual,
                        expected,
                    )

                table = pa.Table.from_arrays(
                    list(table.columns),
                    schema=removed_schema,
                )

                assert table.schema == removed_schema

                removed_writer.write_table(table)

                removed_count += (
                    removed_batch.num_rows
                )

            source_seen += batch.num_rows

            if source_seen % 150000 < batch.num_rows:
                print(
                    "Canonical chains processed:",
                    f"{source_seen:,}/578,524",
                    flush=True,
                )

    finally:
        retained_writer.close()
        removed_writer.close()

    # ========================================================
    # Population and membership audit.
    # ========================================================

    assert source_seen == 578_524
    assert retained_count == 499_770
    assert removed_count == 78_754
    assert m1_count == 764

    assert (
        retained_count + removed_count
        == source_seen
    )

    assert removed_found == set(
        removed_by_key
    )

    assert representative_found == (
        representative_keys
    )

    assert (
        pq.ParquetFile(
            retained_output
        ).metadata.num_rows
        == 499_770
    )

    assert (
        pq.ParquetFile(
            removed_output
        ).metadata.num_rows
        == 78_754
    )

    # ========================================================
    # Copy immutable audit/provenance artifacts.
    # ========================================================

    copies = {
        audit_dir / "representative_mapping.parquet":
            mapping_path,

        audit_dir / "representatives.parquet":
            representatives_path,

        audit_dir / "component_summary.parquet":
            component_summary,

        audit_dir / "near_duplicate_edges_mge2.parquet":
            edges_mge2,

        audit_dir / "near_duplicate_edges_m1_retained.parquet":
            edges_m1,

        provenance_dir / "source_manifest.parquet":
            source_manifest,

        provenance_dir / "representative_policy_v1.yaml":
            args.policy_config,

        provenance_dir / "stage8_full_bri_summary.json":
            stage8_summary,

        provenance_dir / "stage14_graph_summary.json":
            graph_summary,

        provenance_dir / "stage14_representative_summary.json":
            representative_summary,
    }

    for destination, source in copies.items():
        shutil.copy2(
            source,
            destination,
        )

        assert (
            sha256(destination)
            == sha256(source)
        )

    # ========================================================
    # Release manifest.
    # ========================================================

    artifacts = []

    for path in sorted(
        p
        for p in args.output_dir.rglob("*")
        if p.is_file()
    ):
        relative = path.relative_to(
            args.output_dir
        )

        artifacts.append(
            {
                "path":
                    str(relative),

                "bytes":
                    path.stat().st_size,

                "sha256":
                    sha256(path),
            }
        )

    snapshot = (
        root.parent.name
    )

    protocol = (
        root.name
    )

    manifest = {
        "release_name":
            args.output_dir.name,

        "release_version":
            "1.0",

        "snapshot":
            snapshot,

        "protocol":
            protocol,

        "model_scope":
            "model_1",

        "canonical_source":
            str(canonical),

        "canonical_source_sha256":
            sha256(canonical),

        "canonical_input_chain_count":
            578_524,

        "removed_chain_count":
            78_754,

        "retained_chain_count":
            499_770,

        "m1_input_chain_count":
            764,

        "m1_removed_chain_count":
            0,

        "m1_retained_chain_count":
            764,

        "near_duplicate_relation":
            "complete_BRI_L_infinity",

        "distance_representation":
            "exact_integer_milliangstrom",

        "near_duplicate_threshold":
            "d_bri_mA <= 10",

        "representative_policy":
            "comp702_stage14_representative_selection",

        "representative_policy_version":
            "1.0",

        "representative_policy_sha256":
            sha256(args.policy_config),

        "every_removed_chain_has_direct_representative_edge":
            True,

        "connectedness_treated_as_duplicate_equivalence":
            False,

        "m1_deduplication_performed":
            False,

        "stage13_review_subset_used_as_global_edge_set":
            False,

        "old_snapshot_comparison_used":
            False,

        "automatic_transitive_removal":
            False,

        "frozen_upstream_provenance": {
            "stage8_success_sha256":
                sha256(stage8_success),

            "stage13_success_sha256":
                sha256(stage13_success),

            "stage14_graph_success_sha256":
                sha256(graph_success),

            "stage14_representative_success_sha256":
                sha256(representative_success),

            "stage14_representative_mapping_sha256":
                sha256(mapping_path),
        },

        "artifacts":
            artifacts,
    }

    manifest_path = (
        args.output_dir
        / "release_manifest.json"
    )

    manifest_path.write_text(
        json.dumps(
            manifest,
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    # ========================================================
    # Final _SUCCESS marker.
    # ========================================================

    success = {
        "status":
            "PASS",

        "release_name":
            args.output_dir.name,

        "canonical_input_chains":
            578_524,

        "retained_chains":
            499_770,

        "removed_chains":
            78_754,

        "m1_retained":
            764,

        "direct_edge_safety":
            True,

        "connectedness_used_as_equivalence":
            False,

        "old_snapshot_comparison_used":
            False,

        "release_manifest_sha256":
            sha256(manifest_path),
    }

    success_path = (
        args.output_dir
        / "_SUCCESS"
    )

    success_path.write_text(
        json.dumps(
            success,
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    # ========================================================
    # Final report.
    # ========================================================

    print()
    print("===== PDBCLEAN FINAL RELEASE =====")
    print("Release:", args.output_dir)
    print()
    print(
        "Canonical input chains :",
        f"{578_524:,}",
    )
    print(
        "Removed chains         :",
        f"{78_754:,}",
    )
    print(
        "Retained chains        :",
        f"{499_770:,}",
    )
    print(
        "m=1 retained           :",
        f"{764:,}",
    )
    print()
    print(
        "Direct representative edge for every removal:",
        "YES",
    )
    print(
        "Connectedness treated as equivalence:",
        "NO",
    )
    print(
        "Old-snapshot comparison used:",
        "NO",
    )
    print()
    print(
        "Retained dataset:",
        retained_output,
    )
    print(
        "Removed-chain audit:",
        removed_output,
    )
    print(
        "Release manifest:",
        manifest_path,
    )
    print(
        "_SUCCESS:",
        success_path,
    )
    print()
    print(
        "Retained SHA256:",
        sha256(retained_output),
    )
    print(
        "Removed audit SHA256:",
        sha256(removed_output),
    )
    print(
        "Manifest SHA256:",
        sha256(manifest_path),
    )
    print(
        "_SUCCESS SHA256:",
        sha256(success_path),
    )
    print()
    print(
        "FINAL PDBCLEAN RELEASE: PASS"
    )


if __name__ == "__main__":
    main()
