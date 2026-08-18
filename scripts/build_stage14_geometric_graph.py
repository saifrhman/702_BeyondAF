#!/usr/bin/env python3

import argparse
import hashlib
import json
from array import array
from collections import Counter, defaultdict
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq


SCHEMA_VERSION = "1.0"


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--edges", type=Path, required=True)
    p.add_argument("--m1-edges", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--expected-edges", type=int, required=True)
    p.add_argument("--expected-m1-edges", type=int, required=True)
    p.add_argument("--expected-mge2-nodes", type=int, required=True)
    return p.parse_args()


def main():
    args = parse_args()

    assert args.edges.is_file()
    assert args.m1_edges.is_file()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    source_rows = pq.ParquetFile(args.edges).metadata.num_rows
    m1_rows = pq.ParquetFile(args.m1_edges).metadata.num_rows

    assert source_rows == args.expected_edges
    assert m1_rows == args.expected_m1_edges

    # --------------------------------------------------------
    # Dynamic union-find.
    # --------------------------------------------------------

    key_to_id = {}
    id_to_key = []
    node_m = []

    parent = array("I")
    rank = array("B")

    edge_u = array("I")
    edge_v = array("I")

    seen_edges = set()
    distance_hist = Counter()

    def node_id(key, m):
        existing = key_to_id.get(key)

        if existing is not None:
            assert node_m[existing] == m
            return existing

        idx = len(id_to_key)

        assert idx < 2**32

        key_to_id[key] = idx
        id_to_key.append(key)
        node_m.append(m)

        parent.append(idx)
        rank.append(0)

        return idx

    def find(x):
        root = x

        while parent[root] != root:
            root = parent[root]

        while parent[x] != x:
            nxt = parent[x]
            parent[x] = root
            x = nxt

        return root

    def union(a, b):
        ra = find(a)
        rb = find(b)

        if ra == rb:
            return

        if rank[ra] < rank[rb]:
            ra, rb = rb, ra

        parent[rb] = ra

        if rank[ra] == rank[rb]:
            rank[ra] += 1

    # --------------------------------------------------------
    # Stream frozen Stage-8 edges.
    # --------------------------------------------------------

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
    ]

    processed = 0

    print("===== STAGE-14 GEOMETRIC GRAPH BUILD =====", flush=True)
    print("Source edges:", f"{source_rows:,}", flush=True)

    pf = pq.ParquetFile(args.edges)

    for batch in pf.iter_batches(
        batch_size=65536,
        columns=columns,
    ):
        d = batch.to_pydict()
        n = batch.num_rows

        for i in range(n):
            m = int(d["retained_residue_count"][i])
            dist = int(d["d_bri_mA"][i])

            assert m >= 2
            assert 0 <= dist <= 10

            qkey = (
                str(d["query_snapshot"][i]),
                str(d["query_pdb_id"][i]).lower(),
                int(d["query_model_id"][i]),
                str(d["query_label_chain_id"][i]),
            )

            skey = (
                str(d["subject_snapshot"][i]),
                str(d["subject_pdb_id"][i]).lower(),
                int(d["subject_model_id"][i]),
                str(d["subject_label_chain_id"][i]),
            )

            assert qkey != skey

            u = node_id(qkey, m)
            v = node_id(skey, m)

            a = min(u, v)
            b = max(u, v)

            encoded = (a << 32) | b

            assert encoded not in seen_edges, (
                "Duplicate Stage-14 edge",
                qkey,
                skey,
            )

            seen_edges.add(encoded)
            edge_u.append(a)
            edge_v.append(b)

            union(a, b)

            distance_hist[dist] += 1
            processed += 1

        if processed % 250000 < n:
            print(
                "Processed:",
                f"{processed:,}/{source_rows:,}",
                flush=True,
            )

    assert processed == args.expected_edges
    assert len(seen_edges) == args.expected_edges

    # Free duplicate-detection structure before component analysis.
    del seen_edges

    # --------------------------------------------------------
    # Final path compression.
    # --------------------------------------------------------

    for i in range(len(parent)):
        parent[i] = find(i)

    edge_node_count = len(id_to_key)

    assert edge_node_count <= args.expected_mge2_nodes

    no_edge_node_count = (
        args.expected_mge2_nodes
        - edge_node_count
    )

    # --------------------------------------------------------
    # Component node statistics.
    # --------------------------------------------------------

    component_nodes = Counter()
    component_min_key = {}
    component_m = {}

    for idx, key in enumerate(id_to_key):
        root = parent[idx]

        component_nodes[root] += 1

        if (
            root not in component_min_key
            or key < component_min_key[root]
        ):
            component_min_key[root] = key

        if root in component_m:
            assert component_m[root] == node_m[idx]
        else:
            component_m[root] = node_m[idx]

    # --------------------------------------------------------
    # Component edge statistics.
    # --------------------------------------------------------

    component_edges = Counter()

    for u, v in zip(edge_u, edge_v):
        ru = parent[u]
        rv = parent[v]

        assert ru == rv

        component_edges[ru] += 1

    assert sum(component_edges.values()) == args.expected_edges
    assert set(component_edges) == set(component_nodes)

    # --------------------------------------------------------
    # Deterministic component IDs.
    # --------------------------------------------------------

    roots = sorted(
        component_nodes,
        key=lambda r: (
            component_m[r],
            component_min_key[r],
        ),
    )

    root_to_component = {}

    per_length_counter = Counter()

    for root in roots:
        m = component_m[root]
        per_length_counter[m] += 1

        root_to_component[root] = (
            f"M{m:05d}_C{per_length_counter[m]:06d}"
        )

    # --------------------------------------------------------
    # Component summary.
    # --------------------------------------------------------

    component_rows = []

    clique_count = 0
    nonclique_count = 0

    for root in roots:
        n = component_nodes[root]
        e = component_edges[root]

        possible = n * (n - 1) // 2
        is_clique = e == possible

        if is_clique:
            clique_count += 1
        else:
            nonclique_count += 1

        minkey = component_min_key[root]

        component_rows.append(
            {
                "component_id":
                    root_to_component[root],
                "retained_residue_count":
                    component_m[root],
                "node_count":
                    n,
                "edge_count":
                    e,
                "possible_edge_count":
                    possible,
                "edge_density":
                    e / possible,
                "is_clique":
                    is_clique,
                "minimum_member_snapshot":
                    minkey[0],
                "minimum_member_pdb_id":
                    minkey[1],
                "minimum_member_model_id":
                    minkey[2],
                "minimum_member_label_chain_id":
                    minkey[3],
            }
        )

    # --------------------------------------------------------
    # Node -> component mapping for all edge-participating nodes.
    # --------------------------------------------------------

    node_rows = []

    for idx, key in enumerate(id_to_key):
        root = parent[idx]

        node_rows.append(
            {
                "snapshot": key[0],
                "pdb_id": key[1],
                "model_id": key[2],
                "label_chain_id": key[3],
                "retained_residue_count": node_m[idx],
                "component_id":
                    root_to_component[root],
            }
        )

    node_rows.sort(
        key=lambda r: (
            r["retained_residue_count"],
            r["pdb_id"],
            r["model_id"],
            r["label_chain_id"],
        )
    )

    # --------------------------------------------------------
    # Write deterministic graph artifacts.
    # --------------------------------------------------------

    component_path = (
        args.output_dir
        / "component_summary.parquet"
    )

    nodes_path = (
        args.output_dir
        / "edge_node_components.parquet"
    )

    pq.write_table(
        pa.Table.from_pylist(component_rows),
        component_path,
        compression="zstd",
    )

    pq.write_table(
        pa.Table.from_pylist(node_rows),
        nodes_path,
        compression="zstd",
    )

    component_size_hist = Counter(
        r["node_count"]
        for r in component_rows
    )

    summary = {
        "schema_version": SCHEMA_VERSION,
        "relation":
            "paper_faithful_complete_BRI_near_duplicate",
        "threshold":
            "d_bri_mA <= 10",
        "minimum_chain_length":
            2,
        "source_edge_count":
            args.expected_edges,
        "excluded_m1_edge_count":
            args.expected_m1_edges,
        "canonical_mge2_node_count":
            args.expected_mge2_nodes,
        "edge_participating_node_count":
            edge_node_count,
        "no_edge_node_count":
            no_edge_node_count,
        "edge_component_count":
            len(component_rows),
        "singleton_no_edge_component_count":
            no_edge_node_count,
        "total_component_count_including_singletons":
            len(component_rows) + no_edge_node_count,
        "clique_edge_component_count":
            clique_count,
        "nonclique_edge_component_count":
            nonclique_count,
        "distance_histogram_mA": {
            str(k): distance_hist[k]
            for k in sorted(distance_hist)
        },
        "edge_component_size_histogram": {
            str(k): component_size_hist[k]
            for k in sorted(component_size_hist)
        },
        "connectedness_treated_as_duplicate_equivalence":
            False,
        "direct_representative_edge_required":
            True,
        "representative_selection_performed":
            False,
        "removal_decisions_made":
            0,
        "stage13_review_subset_used_as_graph":
            False,
        "old_snapshot_comparison_used":
            False,
        "source_edges_sha256":
            sha256(args.edges),
        "excluded_m1_edges_sha256":
            sha256(args.m1_edges),
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
                "schema_version": SCHEMA_VERSION,
                "edge_count": args.expected_edges,
                "edge_node_count": edge_node_count,
                "edge_component_count": len(component_rows),
                "representative_selection_performed": False,
                "removal_decisions_made": 0,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    # --------------------------------------------------------
    # Report.
    # --------------------------------------------------------

    ranked = sorted(
        component_rows,
        key=lambda r: (
            -r["node_count"],
            -r["edge_count"],
            r["component_id"],
        ),
    )

    print()
    print("===== STAGE-14 GRAPH SUMMARY =====")
    print(
        "m>=2 canonical chains:",
        f"{args.expected_mge2_nodes:,}",
    )
    print(
        "Near-duplicate edges:",
        f"{args.expected_edges:,}",
    )
    print(
        "Chains with >=1 edge:",
        f"{edge_node_count:,}",
    )
    print(
        "Chains with no edge:",
        f"{no_edge_node_count:,}",
    )
    print(
        "Components with edges:",
        f"{len(component_rows):,}",
    )
    print(
        "Clique components:",
        f"{clique_count:,}",
    )
    print(
        "Non-clique components:",
        f"{nonclique_count:,}",
    )

    print()
    print("===== 20 LARGEST EDGE COMPONENTS =====")

    for r in ranked[:20]:
        print(
            r["component_id"],
            f"| m={r['retained_residue_count']}",
            f"| nodes={r['node_count']}",
            f"| edges={r['edge_count']}",
            f"| density={r['edge_density']:.6f}",
            f"| clique={r['is_clique']}",
            f"| min={r['minimum_member_pdb_id']}:"
            f"{r['minimum_member_label_chain_id']}",
        )

    print()
    print("Component summary:", component_path)
    print("Node mapping:", nodes_path)
    print("Global summary:", summary_path)
    print()
    print("m=1 deduplication performed: NO")
    print("Representative selection performed: NO")
    print("Removal decisions made: 0")
    print("Frozen Stage-13 modified: NO")
    print("Old-snapshot comparison used: NO")
    print("STAGE-14 GEOMETRIC GRAPH: PASS")


if __name__ == "__main__":
    main()
