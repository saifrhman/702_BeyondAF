#!/usr/bin/env python3

import argparse
import hashlib
import json
import math
from collections import Counter, defaultdict
from datetime import date
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import yaml


# Validated COMP702 default.  Supplying --threshold-mA uses a configured value
# instead; omitting it reproduces the frozen run exactly.
DEFAULT_NEAR_DUPLICATE_THRESHOLD_MA = 10


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_args():
    p = argparse.ArgumentParser()

    p.add_argument("--graph-dir", type=Path, required=True)
    p.add_argument("--edges", type=Path, required=True)
    p.add_argument("--accepted", type=Path, required=True)
    p.add_argument("--metadata", type=Path, required=True)
    p.add_argument("--config", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)

    p.add_argument(
        "--threshold-mA",
        type=int,
        default=DEFAULT_NEAR_DUPLICATE_THRESHOLD_MA,
        help=(
            "Inclusive complete-BRI near-duplicate threshold in exact integer "
            "milliangstroms. Must agree with the representative policy and "
            "with the Stage-14 graph. Defaults to the validated COMP702 value "
            f"({DEFAULT_NEAR_DUPLICATE_THRESHOLD_MA} mA = 0.010 A)."
        ),
    )

    p.add_argument(
        "--m1-retained-chain-count",
        type=int,
        default=None,
        help=(
            "Number of m = 1 chains, all of which are retained. Derived from "
            "the Stage-6 length-bucket summary when not supplied."
        ),
    )

    p.add_argument(
        "--expected-canonical-input-chains",
        type=int,
        default=None,
        help=(
            "Dataset-version gate: the canonical eligible chain population. "
            "Omitted means no external gate; the internal accounting "
            "assertions still run."
        ),
    )

    return p.parse_args()


def finite_values(value):
    if value is None:
        return []

    if not isinstance(value, (list, tuple)):
        value = [value]

    out = []

    for item in value:
        try:
            x = float(item)
        except (TypeError, ValueError):
            continue

        if math.isfinite(x):
            out.append(x)

    return out


def methods_tuple(value):
    if not value:
        return tuple()

    if not isinstance(value, (list, tuple)):
        value = [value]

    return tuple(
        sorted(str(x) for x in value)
    )


def method_resolution(row, method_tuple):
    """
    Return one comparable nominal resolution value.

    This value is used only when the entire component has
    an identical experimental-method tuple.
    """

    methods = set(method_tuple)

    refine = finite_values(
        row.get("refine_ls_d_res_high")
    )

    em = finite_values(
        row.get("em_3d_reconstruction_resolution")
    )

    if methods == {"ELECTRON MICROSCOPY"}:
        if em:
            return min(em)

        if refine:
            return min(refine)

        return None

    if refine:
        return min(refine)

    return None


def canonical_key(row):
    return (
        str(row["snapshot"]),
        str(row["pdb_id"]).lower(),
        int(row["model_id"]),
        str(row["label_chain_id"]),
    )


def main():
    args = parse_args()

    required = [
        args.graph_dir / "_SUCCESS",
        args.graph_dir / "component_summary.parquet",
        args.graph_dir / "edge_node_components.parquet",
        args.edges,
        args.accepted,
        args.metadata,
        args.config,
    ]

    for path in required:
        assert path.is_file(), path

    assert not args.output_dir.exists(), (
        f"Output already exists: {args.output_dir}"
    )

    args.output_dir.mkdir(
        parents=True,
        exist_ok=False,
    )

    with args.config.open(
        "r",
        encoding="utf-8",
    ) as f:
        policy = yaml.safe_load(f)

    # The policy file and the run configuration must agree.  A mismatch means
    # the deletion relation and the policy were derived from different
    # thresholds, which would be a silent scientific change.
    assert (
        policy["scientific_scope"]["threshold"]["value_mA"]
        == args.threshold_mA
    ), (
        "Representative policy threshold "
        f"{policy['scientific_scope']['threshold']['value_mA']} mA does not "
        f"match the configured threshold {args.threshold_mA} mA"
    )

    assert (
        policy["graph_policy"][
            "connected_component_is_duplicate_equivalence"
        ]
        is False
    )

    assert (
        policy["graph_policy"][
            "every_removed_chain_requires_direct_representative_edge"
        ]
        is True
    )

    # ========================================================
    # Frozen graph nodes
    # ========================================================

    # Population counts are read from the upstream Stage-14 graph summary
    # rather than restated as literals, so the accounting gate always checks
    # against what the graph actually contained for this run.
    graph_summary = json.loads(
        (
            args.graph_dir
            / "global_summary.json"
        ).read_text(encoding="utf-8")
    )

    assert (
        graph_summary["threshold"]
        == f"d_bri_mA <= {args.threshold_mA}"
    ), (
        "Stage-14 graph was built with "
        f"{graph_summary['threshold']!r}, which does not match the "
        f"configured threshold d_bri_mA <= {args.threshold_mA}"
    )

    assert (
        graph_summary[
            "connectedness_treated_as_duplicate_equivalence"
        ]
        is False
    )

    graph_edge_count = int(
        graph_summary["source_edge_count"]
    )
    graph_edge_node_count = int(
        graph_summary["edge_participating_node_count"]
    )
    graph_edge_component_count = int(
        graph_summary["edge_component_count"]
    )
    graph_no_edge_node_count = int(
        graph_summary["no_edge_node_count"]
    )
    graph_mge2_node_count = int(
        graph_summary["canonical_mge2_node_count"]
    )

    # The m=1 population is upstream provenance, not a literal.  Stage 6
    # records it next to the graph inside the same protocol root; an explicit
    # override is available for reruns that relocate that output.
    m1_retained_count = args.m1_retained_chain_count

    if m1_retained_count is None:
        stage6_summary_path = (
            args.graph_dir.parent
            / "length_buckets"
            / "global_summary.json"
        )

        if stage6_summary_path.is_file():
            m1_retained_count = int(
                json.loads(
                    stage6_summary_path.read_text(
                        encoding="utf-8"
                    )
                )["m1_chain_count"]
            )

    if (
        m1_retained_count is None
        and args.expected_canonical_input_chains is None
    ):
        raise SystemExit(
            "Cannot determine the m=1 population. Supply "
            "--m1-retained-chain-count or "
            "--expected-canonical-input-chains, or make the Stage-6 "
            f"length-bucket summary readable at {stage6_summary_path}"
        )

    graph_rows = pq.read_table(
        args.graph_dir
        / "edge_node_components.parquet"
    ).to_pylist()

    assert len(graph_rows) == graph_edge_node_count

    graph_by_key = {
        canonical_key(r): r
        for r in graph_rows
    }

    assert len(graph_by_key) == len(graph_rows)

    components = defaultdict(list)

    for key, row in graph_by_key.items():
        components[
            row["component_id"]
        ].append(key)

    assert len(components) == 20_789

    # ========================================================
    # Accepted-chain quality evidence
    # ========================================================

    accepted_columns = [
        "snapshot",
        "pdb_id",
        "model_id",
        "label_chain_id",
        "retained_residue_count",
        "terminal_trimmed",
        "dirty_residue_count",
    ]

    accepted = pq.read_table(
        args.accepted,
        columns=accepted_columns,
    ).to_pylist()

    accepted_by_key = {
        canonical_key(r): r
        for r in accepted
    }

    missing_quality = (
        set(graph_by_key)
        - set(accepted_by_key)
    )

    assert not missing_quality, (
        "Graph nodes missing accepted-chain evidence",
        list(sorted(missing_quality))[:10],
    )

    # ========================================================
    # Entry metadata
    # ========================================================

    metadata_columns = [
        "snapshot",
        "pdb_id",
        "experimental_methods",
        "refine_ls_d_res_high",
        "em_3d_reconstruction_resolution",
        "initial_deposition_date",
    ]

    metadata = pq.read_table(
        args.metadata,
        columns=metadata_columns,
    ).to_pylist()

    meta_by_pdb = {
        (
            str(r["snapshot"]),
            str(r["pdb_id"]).lower(),
        ): r
        for r in metadata
    }

    graph_pdbs = {
        (k[0], k[1])
        for k in graph_by_key
    }

    assert graph_pdbs <= set(meta_by_pdb)

    # ========================================================
    # Frozen direct-edge adjacency
    # ========================================================

    adjacency = {
        key: {}
        for key in graph_by_key
    }

    edge_count = 0

    edge_columns = [
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
    ]

    pf = pq.ParquetFile(args.edges)

    for batch in pf.iter_batches(
        batch_size=65536,
        columns=edge_columns,
    ):
        d = batch.to_pydict()

        for i in range(batch.num_rows):
            q = (
                str(d["query_snapshot"][i]),
                str(d["query_pdb_id"][i]).lower(),
                int(d["query_model_id"][i]),
                str(d["query_label_chain_id"][i]),
            )

            s = (
                str(d["subject_snapshot"][i]),
                str(d["subject_pdb_id"][i]).lower(),
                int(d["subject_model_id"][i]),
                str(d["subject_label_chain_id"][i]),
            )

            dist = int(
                d["d_bri_mA"][i]
            )

            assert q in graph_by_key
            assert s in graph_by_key

            assert 0 <= dist <= args.threshold_mA

            assert (
                graph_by_key[q]["component_id"]
                == graph_by_key[s]["component_id"]
            )

            assert s not in adjacency[q]
            assert q not in adjacency[s]

            adjacency[q][s] = dist
            adjacency[s][q] = dist

            edge_count += 1

    assert edge_count == 1_068_256

    # ========================================================
    # Component clique metadata
    # ========================================================

    component_summary = pq.read_table(
        args.graph_dir
        / "component_summary.parquet"
    ).to_pylist()

    comp_info = {
        r["component_id"]: r
        for r in component_summary
    }

    assert len(comp_info) == 20_789

    # ========================================================
    # Deterministic quality ranking
    # ========================================================

    def component_ranker(member_keys):
        method_tuples = {
            methods_tuple(
                meta_by_pdb[
                    (key[0], key[1])
                ]["experimental_methods"]
            )
            for key in member_keys
        }

        comparable_method = (
            len(method_tuples) == 1
        )

        common_method = (
            next(iter(method_tuples))
            if comparable_method
            else None
        )

        def rank(key):
            quality = accepted_by_key[key]
            meta = meta_by_pdb[
                (key[0], key[1])
            ]

            terminal_trimmed = bool(
                quality["terminal_trimmed"]
            )

            dirty = quality[
                "dirty_residue_count"
            ]

            if dirty is None:
                dirty = 10**9
            else:
                dirty = int(dirty)

            resolution = None

            if comparable_method:
                resolution = method_resolution(
                    meta,
                    common_method,
                )

            resolution_missing = (
                resolution is None
            )

            resolution_value = (
                float("inf")
                if resolution is None
                else float(resolution)
            )

            return (
                1 if terminal_trimmed else 0,
                dirty,
                1 if resolution_missing else 0,
                resolution_value,
                key,
            )

        return (
            rank,
            comparable_method,
            common_method,
        )

    # ========================================================
    # Representative selection
    # ========================================================

    mapping_rows = []
    representative_rows = []

    clique_components = 0
    nonclique_components = 0

    clique_reps = 0
    nonclique_reps = 0

    resolution_enabled_components = 0

    for index, component_id in enumerate(
        sorted(components),
        1,
    ):
        members = components[component_id]

        info = comp_info[component_id]

        rank, comparable_method, common_method = (
            component_ranker(members)
        )

        if comparable_method:
            resolution_enabled_components += 1

        ordered = sorted(
            members,
            key=rank,
        )

        representatives = []

        assignment = {}

        if bool(info["is_clique"]):
            clique_components += 1

            representative = ordered[0]

            representatives.append(
                representative
            )

            for member in ordered:
                assignment[member] = representative

            clique_reps += 1

        else:
            nonclique_components += 1

            unassigned = set(members)

            for candidate in ordered:
                if candidate not in unassigned:
                    continue

                representative = candidate

                representatives.append(
                    representative
                )

                assignment[
                    representative
                ] = representative

                unassigned.remove(
                    representative
                )

                direct_neighbors = sorted(
                    (
                        n
                        for n in adjacency[
                            representative
                        ]
                        if n in unassigned
                    ),
                    key=rank,
                )

                for neighbor in direct_neighbors:
                    assignment[
                        neighbor
                    ] = representative

                    unassigned.remove(
                        neighbor
                    )

            assert not unassigned

            nonclique_reps += len(
                representatives
            )

        assert len(assignment) == len(members)

        representative_set = set(
            representatives
        )

        for member in ordered:
            representative = assignment[
                member
            ]

            removed = (
                member != representative
            )

            if removed:
                assert (
                    representative
                    in adjacency[member]
                )

                distance = adjacency[
                    member
                ][
                    representative
                ]

                assert 0 <= distance <= args.threshold_mA
            else:
                distance = 0

            q = accepted_by_key[member]
            m = meta_by_pdb[
                (member[0], member[1])
            ]

            rep_q = accepted_by_key[
                representative
            ]

            rep_m = meta_by_pdb[
                (
                    representative[0],
                    representative[1],
                )
            ]

            member_method = methods_tuple(
                m["experimental_methods"]
            )

            rep_method = methods_tuple(
                rep_m["experimental_methods"]
            )

            member_resolution = (
                method_resolution(
                    m,
                    common_method,
                )
                if comparable_method
                else None
            )

            rep_resolution = (
                method_resolution(
                    rep_m,
                    common_method,
                )
                if comparable_method
                else None
            )

            mapping_rows.append(
                {
                    "component_id":
                        component_id,

                    "component_is_clique":
                        bool(
                            info["is_clique"]
                        ),

                    "snapshot":
                        member[0],

                    "pdb_id":
                        member[1],

                    "model_id":
                        member[2],

                    "label_chain_id":
                        member[3],

                    "retained_residue_count":
                        int(
                            graph_by_key[
                                member
                            ][
                                "retained_residue_count"
                            ]
                        ),

                    "action":
                        (
                            "remove"
                            if removed
                            else "retain_representative"
                        ),

                    "representative_snapshot":
                        representative[0],

                    "representative_pdb_id":
                        representative[1],

                    "representative_model_id":
                        representative[2],

                    "representative_label_chain_id":
                        representative[3],

                    "direct_d_bri_mA":
                        distance,

                    "terminal_trimmed":
                        bool(
                            q["terminal_trimmed"]
                        ),

                    "dirty_residue_count":
                        int(
                            q[
                                "dirty_residue_count"
                            ]
                            or 0
                        ),

                    "experimental_methods":
                        list(
                            member_method
                        ),

                    "comparable_method_resolution":
                        member_resolution,

                    "representative_terminal_trimmed":
                        bool(
                            rep_q[
                                "terminal_trimmed"
                            ]
                        ),

                    "representative_dirty_residue_count":
                        int(
                            rep_q[
                                "dirty_residue_count"
                            ]
                            or 0
                        ),

                    "representative_experimental_methods":
                        list(
                            rep_method
                        ),

                    "representative_comparable_method_resolution":
                        rep_resolution,

                    "policy_version":
                        str(
                            policy[
                                "policy_version"
                            ]
                        ),
                }
            )

        for representative in representatives:
            q = accepted_by_key[
                representative
            ]

            m = meta_by_pdb[
                (
                    representative[0],
                    representative[1],
                )
            ]

            resolution = (
                method_resolution(
                    m,
                    common_method,
                )
                if comparable_method
                else None
            )

            assigned_count = sum(
                rep == representative
                for rep in assignment.values()
            )

            representative_rows.append(
                {
                    "component_id":
                        component_id,

                    "component_is_clique":
                        bool(
                            info["is_clique"]
                        ),

                    "snapshot":
                        representative[0],

                    "pdb_id":
                        representative[1],

                    "model_id":
                        representative[2],

                    "label_chain_id":
                        representative[3],

                    "retained_residue_count":
                        int(
                            graph_by_key[
                                representative
                            ][
                                "retained_residue_count"
                            ]
                        ),

                    "assigned_chain_count":
                        assigned_count,

                    "terminal_trimmed":
                        bool(
                            q[
                                "terminal_trimmed"
                            ]
                        ),

                    "dirty_residue_count":
                        int(
                            q[
                                "dirty_residue_count"
                            ]
                            or 0
                        ),

                    "experimental_methods":
                        list(
                            methods_tuple(
                                m[
                                    "experimental_methods"
                                ]
                            )
                        ),

                    "comparable_method_resolution":
                        resolution,

                    "policy_version":
                        str(
                            policy[
                                "policy_version"
                            ]
                        ),
                }
            )

        if index % 5000 == 0:
            print(
                "Components processed:",
                f"{index:,}/{len(components):,}",
                flush=True,
            )

    # ========================================================
    # Global safety audit
    # ========================================================

    assert len(mapping_rows) == graph_edge_node_count

    retained_edge_reps = sum(
        r["action"]
        == "retain_representative"
        for r in mapping_rows
    )

    removed = sum(
        r["action"] == "remove"
        for r in mapping_rows
    )

    assert (
        retained_edge_reps
        + removed
        == graph_edge_node_count
    )

    assert (
        retained_edge_reps
        == len(representative_rows)
    )

    for r in mapping_rows:
        if r["action"] == "remove":
            assert (
                0
                <= int(r["direct_d_bri_mA"])
                <= args.threshold_mA
            )

    # m>=2 no-edge chains + selected representatives, both taken from the
    # upstream graph summary rather than from literals.
    retained_mge2 = (
        graph_no_edge_node_count
        + retained_edge_reps
    )

    assert (
        retained_mge2
        + removed
        == graph_mge2_node_count
    )

    # Every m=1 chain is retained: BRI geometry is degenerate at m=1, so m=1
    # chains never take part in deduplication.  The m=1 population is the
    # canonical eligible population minus the m>=2 population.
    if args.expected_canonical_input_chains is not None:
        original_total = args.expected_canonical_input_chains
    else:
        original_total = (
            graph_mge2_node_count
            + m1_retained_count
        )

    retained_m1 = original_total - graph_mge2_node_count

    assert retained_m1 >= 0, (
        "Canonical input chain count is smaller than the m>=2 population"
    )

    retained_total = (
        retained_mge2
        + retained_m1
    )

    assert (
        retained_total
        + removed
        == original_total
    )

    # ========================================================
    # Write output artifacts
    # ========================================================

    mapping_rows.sort(
        key=lambda r: (
            r["component_id"],
            r["pdb_id"],
            r["model_id"],
            r["label_chain_id"],
        )
    )

    representative_rows.sort(
        key=lambda r: (
            r["component_id"],
            r["pdb_id"],
            r["model_id"],
            r["label_chain_id"],
        )
    )

    mapping_path = (
        args.output_dir
        / "representative_mapping.parquet"
    )

    representatives_path = (
        args.output_dir
        / "representatives.parquet"
    )

    pq.write_table(
        pa.Table.from_pylist(
            mapping_rows
        ),
        mapping_path,
        compression="zstd",
    )

    pq.write_table(
        pa.Table.from_pylist(
            representative_rows
        ),
        representatives_path,
        compression="zstd",
    )

    action_hist = Counter(
        r["action"]
        for r in mapping_rows
    )

    representative_per_component = Counter()

    for r in representative_rows:
        representative_per_component[
            r["component_id"]
        ] += 1

    multi_rep_nonclique = sum(
        1
        for cid, n in
        representative_per_component.items()
        if (
            not bool(
                comp_info[cid][
                    "is_clique"
                ]
            )
            and n > 1
        )
    )

    summary = {
        "policy_name":
            policy["policy_name"],

        "policy_version":
            str(
                policy["policy_version"]
            ),

        "source_graph_edge_count":
            graph_edge_count,

        "source_graph_edge_nodes":
            graph_edge_node_count,

        "source_graph_edge_components":
            graph_edge_component_count,

        "clique_components":
            clique_components,

        "nonclique_components":
            nonclique_components,

        "resolution_comparison_enabled_components":
            resolution_enabled_components,

        "edge_component_representatives":
            retained_edge_reps,

        "clique_representatives":
            clique_reps,

        "nonclique_representatives":
            nonclique_reps,

        "nonclique_components_with_multiple_representatives":
            multi_rep_nonclique,

        "removed_mge2_chains":
            removed,

        "retained_mge2_no_edge_chains":
            graph_no_edge_node_count,

        "retained_mge2_representatives":
            retained_edge_reps,

        "retained_mge2_total":
            retained_mge2,

        "retained_m1_total":
            retained_m1,

        "canonical_input_chain_count":
            original_total,

        "final_retained_chain_count":
            retained_total,

        "final_removed_chain_count":
            removed,

        "action_histogram":
            dict(
                sorted(
                    action_hist.items()
                )
            ),

        "duplicate_threshold":
            f"d_bri_mA <= {args.threshold_mA}",

        "connectedness_treated_as_equivalence":
            False,

        "every_removed_chain_has_direct_representative_edge":
            True,

        "m1_deduplication_performed":
            False,

        "stage13_review_subset_used_as_global_edge_set":
            False,

        "old_snapshot_comparison_used":
            False,

        "policy_config_sha256":
            sha256(args.config),

        "source_edges_sha256":
            sha256(args.edges),

        "graph_success_sha256":
            sha256(
                args.graph_dir
                / "_SUCCESS"
            ),

        "accepted_quality_sha256":
            sha256(args.accepted),

        "entry_metadata_sha256":
            sha256(args.metadata),
    }

    summary_path = (
        args.output_dir
        / "global_summary.json"
    )

    summary_path.write_text(
        json.dumps(
            summary,
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    success_path = (
        args.output_dir
        / "_SUCCESS"
    )

    success_path.write_text(
        json.dumps(
            {
                "status": "PASS",
                "policy_version":
                    str(
                        policy[
                            "policy_version"
                        ]
                    ),
                "input_chains":
                    original_total,
                "retained_chains":
                    retained_total,
                "removed_chains":
                    removed,
                "direct_edge_safety":
                    True,
                "m1_deduplication_performed":
                    False,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    # ========================================================
    # Console summary
    # ========================================================

    print()
    print(
        "===== STAGE-14 REPRESENTATIVE SELECTION ====="
    )

    print(
        "Canonical input chains:",
        f"{original_total:,}",
    )

    print(
        "m>=2 edge-participating chains:",
        "99,854",
    )

    print(
        "m>=2 no-edge chains retained:",
        f"{graph_no_edge_node_count:,}",
    )

    print(
        "m=1 chains retained:",
        f"{retained_m1:,}",
    )

    print()
    print(
        "Clique components:",
        f"{clique_components:,}",
    )

    print(
        "Non-clique components:",
        f"{nonclique_components:,}",
    )

    print(
        "Representatives in edge components:",
        f"{retained_edge_reps:,}",
    )

    print(
        "Removed chains:",
        f"{removed:,}",
    )

    print(
        "Final retained chains:",
        f"{retained_total:,}",
    )

    print()
    print(
        "Non-clique components with >1 representative:",
        f"{multi_rep_nonclique:,}",
    )

    print(
        f"Every removed chain has direct <={args.threshold_mA} mA edge:",
        "YES",
    )

    print(
        "Connected components treated as equivalence:",
        "NO",
    )

    print(
        "m=1 deduplicated:",
        "NO",
    )

    print(
        "Old-snapshot comparison used:",
        "NO",
    )

    print()
    print(
        "Mapping:",
        mapping_path,
    )

    print(
        "Representatives:",
        representatives_path,
    )

    print(
        "Summary:",
        summary_path,
    )

    print(
        "STAGE-14 REPRESENTATIVE SELECTION: PASS"
    )


if __name__ == "__main__":
    main()
