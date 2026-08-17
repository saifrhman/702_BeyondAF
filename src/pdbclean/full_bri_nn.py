"""Exact complete-BRI radius search after Brain filtering.

Scientific cascade
------------------
The input pair relation has already passed the Brain filter.
For each query chain, only its Stage-7 Brain-close subjects are
placed in an Elkin--Kurlin compressed cover tree built on the
flattened complete BRI.

The metric is exactly the L-infinity metric on the complete
m x 9 integer-milliangstrom BRI representation.

The conversion from the Stage-7 pair relation to per-query
reference sets is a COMP702 engineering choice.  It preserves
the Stage-7 candidate relation exactly: no pair rejected by
Brain is introduced into the output search space.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from pdbclean.compressed_cover_tree import CompressedCoverTree


class FullBRINNSearchError(RuntimeError):
    """Raised when exact complete-BRI NN search cannot proceed safely."""


@dataclass(frozen=True)
class CandidateRadiusResult:
    """Hits from a Stage-7 candidate relation."""

    row_indices: np.ndarray
    distances_mA: np.ndarray
    query_neighborhood_count: int
    compressed_tree_count: int


@dataclass(frozen=True)
class AllPairsRadiusResult:
    """Unordered radius hits from one complete BRI bucket."""

    query_indices: np.ndarray
    subject_indices: np.ndarray
    distances_mA: np.ndarray
    unique_point_count: int


def _validate_bri(
    bri: np.ndarray,
) -> np.ndarray:
    array = np.asarray(bri)

    if array.dtype != np.int64:
        raise FullBRINNSearchError(
            "BRI must use exact np.int64 milliangstrom representation"
        )

    if array.ndim != 3:
        raise FullBRINNSearchError(
            "BRI must have shape (chain_count, m, 9)"
        )

    if array.shape[0] < 1:
        raise FullBRINNSearchError(
            "BRI population must contain at least one chain"
        )

    if array.shape[1] < 1:
        raise FullBRINNSearchError(
            "BRI residue count m must be positive"
        )

    if array.shape[2] != 9:
        raise FullBRINNSearchError(
            "BRI final dimension must be 9"
        )

    return array


def _validate_indices(
    values: np.ndarray,
    *,
    name: str,
    chain_count: int,
) -> np.ndarray:
    result = np.asarray(
        values,
        dtype=np.int64,
    )

    if result.ndim != 1:
        raise FullBRINNSearchError(
            f"{name} must be one-dimensional"
        )

    if np.any(result < 0) or np.any(
        result >= chain_count
    ):
        raise FullBRINNSearchError(
            f"{name} contains an out-of-range chain index"
        )

    return result


def search_brain_candidate_pairs(
    bri: np.ndarray,
    query_indices: np.ndarray,
    subject_indices: np.ndarray,
    *,
    radius_mA: int = 10,
) -> CandidateRadiusResult:
    """Search complete BRI only inside the Brain-passed pair relation.

    Each Stage-7 candidate row is an unordered pair represented
    in its existing query/subject orientation.  Rows are grouped
    by query index.  The subjects belonging to that query are the
    reference set for one compressed-cover-tree radius search.

    Identical subject BRI vectors are collapsed before tree
    construction because a metric set cannot contain distinct
    points at zero distance.  Hits are expanded back to all
    original Stage-7 candidate rows afterwards.
    """

    bri = _validate_bri(
        bri
    )

    if (
        not isinstance(radius_mA, (int, np.integer))
        or int(radius_mA) < 0
    ):
        raise FullBRINNSearchError(
            "radius_mA must be a non-negative integer"
        )

    radius_mA = int(radius_mA)

    chain_count = bri.shape[0]

    q_indices = _validate_indices(
        query_indices,
        name="query_indices",
        chain_count=chain_count,
    )
    s_indices = _validate_indices(
        subject_indices,
        name="subject_indices",
        chain_count=chain_count,
    )

    if q_indices.shape != s_indices.shape:
        raise FullBRINNSearchError(
            "Query/subject index shape mismatch"
        )

    pair_count = q_indices.shape[0]

    if pair_count == 0:
        return CandidateRadiusResult(
            row_indices=np.empty(
                0,
                dtype=np.int64,
            ),
            distances_mA=np.empty(
                0,
                dtype=np.int64,
            ),
            query_neighborhood_count=0,
            compressed_tree_count=0,
        )

    if np.any(
        q_indices == s_indices
    ):
        raise FullBRINNSearchError(
            "Candidate relation contains a self-pair"
        )

    flattened = bri.reshape(
        chain_count,
        -1,
    )

    # Stable sorting lets us process each query neighbourhood
    # once while retaining the original Stage-7 row identity.
    order = np.argsort(
        q_indices,
        kind="stable",
    )

    sorted_queries = q_indices[
        order
    ]

    boundaries = np.flatnonzero(
        np.r_[
            True,
            sorted_queries[1:]
            != sorted_queries[:-1],
            True,
        ]
    )

    hit_rows: list[int] = []
    hit_distances: list[int] = []

    neighborhood_count = 0
    tree_count = 0

    for group_number in range(
        boundaries.shape[0] - 1
    ):
        start = int(
            boundaries[group_number]
        )
        stop = int(
            boundaries[group_number + 1]
        )

        rows = order[
            start:stop
        ]

        query_index = int(
            q_indices[rows[0]]
        )

        subjects = s_indices[
            rows
        ]

        if (
            np.unique(subjects).shape[0]
            != subjects.shape[0]
        ):
            raise FullBRINNSearchError(
                "Duplicate subject for one query in "
                "Stage-7 candidate relation"
            )

        neighborhood_count += 1

        subject_points = flattened[
            subjects
        ]

        # np.unique collapses exact complete-BRI duplicates.
        # inverse maps each original candidate subject back to
        # the unique metric point used by the tree.
        unique_points, inverse = np.unique(
            subject_points,
            axis=0,
            return_inverse=True,
        )

        unique_points = np.asarray(
            unique_points,
            dtype=np.int64,
        )
        inverse = np.asarray(
            inverse,
            dtype=np.int64,
        )

        tree = CompressedCoverTree(
            unique_points
        )
        tree_count += 1

        result = tree.radius_neighbors(
            flattened[
                query_index
            ],
            radius_mA,
        )

        if result.indices.shape != (
            result.distances_mA.shape
        ):
            raise FullBRINNSearchError(
                "Compressed-cover-tree result shape mismatch"
            )

        for (
            unique_index,
            distance_mA,
        ) in zip(
            result.indices.tolist(),
            result.distances_mA.tolist(),
        ):
            matching_positions = np.flatnonzero(
                inverse
                == int(unique_index)
            )

            for position in (
                matching_positions.tolist()
            ):
                hit_rows.append(
                    int(rows[position])
                )
                hit_distances.append(
                    int(distance_mA)
                )

    if not hit_rows:
        return CandidateRadiusResult(
            row_indices=np.empty(
                0,
                dtype=np.int64,
            ),
            distances_mA=np.empty(
                0,
                dtype=np.int64,
            ),
            query_neighborhood_count=(
                neighborhood_count
            ),
            compressed_tree_count=tree_count,
        )

    row_array = np.asarray(
        hit_rows,
        dtype=np.int64,
    )
    distance_array = np.asarray(
        hit_distances,
        dtype=np.int64,
    )

    # Restore original Stage-7 row order.
    final_order = np.argsort(
        row_array,
        kind="stable",
    )

    row_array = row_array[
        final_order
    ]
    distance_array = distance_array[
        final_order
    ]

    if np.unique(
        row_array
    ).shape[0] != row_array.shape[0]:
        raise FullBRINNSearchError(
            "A Stage-7 candidate row was emitted more than once"
        )

    if np.any(
        distance_array > radius_mA
    ):
        raise FullBRINNSearchError(
            "Radius search emitted an out-of-radius pair"
        )

    return CandidateRadiusResult(
        row_indices=row_array,
        distances_mA=distance_array,
        query_neighborhood_count=(
            neighborhood_count
        ),
        compressed_tree_count=tree_count,
    )


def search_all_pairs_radius(
    bri: np.ndarray,
    *,
    radius_mA: int = 10,
) -> AllPairsRadiusResult:
    """Find all unordered within-radius pairs in one BRI bucket.

    This is used for a bucket that bypasses Brain, notably m=1.

    Exact duplicate BRI vectors are collapsed for construction of
    the metric tree and expanded back into chain pairs afterwards.
    """

    bri = _validate_bri(
        bri
    )

    if (
        not isinstance(radius_mA, (int, np.integer))
        or int(radius_mA) < 0
    ):
        raise FullBRINNSearchError(
            "radius_mA must be a non-negative integer"
        )

    radius_mA = int(
        radius_mA
    )

    chain_count = bri.shape[0]

    if chain_count < 2:
        return AllPairsRadiusResult(
            query_indices=np.empty(
                0,
                dtype=np.int64,
            ),
            subject_indices=np.empty(
                0,
                dtype=np.int64,
            ),
            distances_mA=np.empty(
                0,
                dtype=np.int64,
            ),
            unique_point_count=chain_count,
        )

    flattened = bri.reshape(
        chain_count,
        -1,
    )

    unique_points, inverse = np.unique(
        flattened,
        axis=0,
        return_inverse=True,
    )

    unique_points = np.asarray(
        unique_points,
        dtype=np.int64,
    )
    inverse = np.asarray(
        inverse,
        dtype=np.int64,
    )

    unique_count = (
        unique_points.shape[0]
    )

    members = [
        np.flatnonzero(
            inverse == i
        ).astype(
            np.int64,
            copy=False,
        )
        for i in range(
            unique_count
        )
    ]

    tree = CompressedCoverTree(
        unique_points
    )

    q_hits: list[int] = []
    s_hits: list[int] = []
    d_hits: list[int] = []

    for unique_query in range(
        unique_count
    ):
        result = tree.radius_neighbors(
            unique_points[
                unique_query
            ],
            radius_mA,
        )

        for (
            unique_subject,
            distance_mA,
        ) in zip(
            result.indices.tolist(),
            result.distances_mA.tolist(),
        ):
            unique_subject = int(
                unique_subject
            )
            distance_mA = int(
                distance_mA
            )

            # Each unordered unique-point relation is
            # expanded exactly once.
            if (
                unique_subject
                < unique_query
            ):
                continue

            left = members[
                unique_query
            ]
            right = members[
                unique_subject
            ]

            if (
                unique_subject
                == unique_query
            ):
                if left.shape[0] < 2:
                    continue

                local_q, local_s = (
                    np.triu_indices(
                        left.shape[0],
                        k=1,
                    )
                )

                for a, b in zip(
                    local_q.tolist(),
                    local_s.tolist(),
                ):
                    q = int(
                        left[a]
                    )
                    s = int(
                        left[b]
                    )

                    q_hits.append(
                        q
                    )
                    s_hits.append(
                        s
                    )
                    d_hits.append(
                        distance_mA
                    )

                continue

            for left_index in (
                left.tolist()
            ):
                for right_index in (
                    right.tolist()
                ):
                    q = int(
                        left_index
                    )
                    s = int(
                        right_index
                    )

                    if q > s:
                        q, s = s, q

                    q_hits.append(
                        q
                    )
                    s_hits.append(
                        s
                    )
                    d_hits.append(
                        distance_mA
                    )

    if not q_hits:
        return AllPairsRadiusResult(
            query_indices=np.empty(
                0,
                dtype=np.int64,
            ),
            subject_indices=np.empty(
                0,
                dtype=np.int64,
            ),
            distances_mA=np.empty(
                0,
                dtype=np.int64,
            ),
            unique_point_count=(
                unique_count
            ),
        )

    q_array = np.asarray(
        q_hits,
        dtype=np.int64,
    )
    s_array = np.asarray(
        s_hits,
        dtype=np.int64,
    )
    d_array = np.asarray(
        d_hits,
        dtype=np.int64,
    )

    order = np.lexsort(
        (
            s_array,
            q_array,
        )
    )

    q_array = q_array[
        order
    ]
    s_array = s_array[
        order
    ]
    d_array = d_array[
        order
    ]

    if np.any(
        q_array >= s_array
    ):
        raise FullBRINNSearchError(
            "All-pairs result is not strictly upper triangular"
        )

    pairs = np.column_stack(
        (
            q_array,
            s_array,
        )
    )

    if np.unique(
        pairs,
        axis=0,
    ).shape[0] != pairs.shape[0]:
        raise FullBRINNSearchError(
            "All-pairs radius search emitted a duplicate pair"
        )

    if np.any(
        d_array > radius_mA
    ):
        raise FullBRINNSearchError(
            "All-pairs radius search emitted an "
            "out-of-radius pair"
        )

    return AllPairsRadiusResult(
        query_indices=q_array,
        subject_indices=s_array,
        distances_mA=d_array,
        unique_point_count=unique_count,
    )


@dataclass(frozen=True)
class ComponentRadiusResult:
    """Complete-BRI hits inside Brain-connected components."""

    row_indices: np.ndarray
    distances_mA: np.ndarray
    participating_chain_count: int
    component_count: int
    compressed_tree_count: int
    non_candidate_radius_hit_count: int


def search_brain_candidate_components(
    bri: np.ndarray,
    query_indices: np.ndarray,
    subject_indices: np.ndarray,
    *,
    radius_mA: int = 10,
) -> ComponentRadiusResult:
    """Search complete BRI inside components of the Brain candidate graph.

    The input edges are the pairs retained by Stage 7.

    One compressed cover tree is built for each connected component
    of that Brain-filtered graph.  Complete-BRI radius hits are then
    intersected with the original Stage-7 edge relation, so component
    grouping can never introduce a pair that failed Brain filtering.

    ``non_candidate_radius_hit_count`` records complete-BRI radius
    hits encountered inside a component that were not Stage-7 edges.
    For the real pipeline this is expected to be zero from
    dBrain <= dBRI.
    """

    bri = _validate_bri(bri)

    if (
        not isinstance(radius_mA, (int, np.integer))
        or int(radius_mA) < 0
    ):
        raise FullBRINNSearchError(
            "radius_mA must be a non-negative integer"
        )

    radius_mA = int(radius_mA)
    chain_count = bri.shape[0]

    q_indices = _validate_indices(
        query_indices,
        name="query_indices",
        chain_count=chain_count,
    )
    s_indices = _validate_indices(
        subject_indices,
        name="subject_indices",
        chain_count=chain_count,
    )

    if q_indices.shape != s_indices.shape:
        raise FullBRINNSearchError(
            "Query/subject index shape mismatch"
        )

    pair_count = q_indices.shape[0]

    if pair_count == 0:
        return ComponentRadiusResult(
            row_indices=np.empty(0, dtype=np.int64),
            distances_mA=np.empty(0, dtype=np.int64),
            participating_chain_count=0,
            component_count=0,
            compressed_tree_count=0,
            non_candidate_radius_hit_count=0,
        )

    if np.any(q_indices == s_indices):
        raise FullBRINNSearchError(
            "Candidate relation contains a self-pair"
        )

    # Canonical unordered pair -> original Stage-7 row.
    edge_to_row: dict[tuple[int, int], int] = {}

    for row, (q, s) in enumerate(
        zip(q_indices.tolist(), s_indices.tolist())
    ):
        q = int(q)
        s = int(s)

        if q > s:
            q, s = s, q

        edge = (q, s)

        if edge in edge_to_row:
            raise FullBRINNSearchError(
                "Duplicate unordered pair in Stage-7 candidate relation"
            )

        edge_to_row[edge] = row

    participating = np.unique(
        np.concatenate(
            (q_indices, s_indices)
        )
    ).astype(np.int64, copy=False)

    # Union-find over canonical chain indices.
    parent = {
        int(index): int(index)
        for index in participating.tolist()
    }
    size = {
        int(index): 1
        for index in participating.tolist()
    }

    def find(index: int) -> int:
        root = index

        while parent[root] != root:
            root = parent[root]

        while parent[index] != index:
            next_index = parent[index]
            parent[index] = root
            index = next_index

        return root

    def union(left: int, right: int) -> None:
        left_root = find(left)
        right_root = find(right)

        if left_root == right_root:
            return

        if size[left_root] < size[right_root]:
            left_root, right_root = (
                right_root,
                left_root,
            )

        parent[right_root] = left_root
        size[left_root] += size[right_root]

    for q, s in zip(
        q_indices.tolist(),
        s_indices.tolist(),
    ):
        union(int(q), int(s))

    component_members: dict[int, list[int]] = {}

    for index in participating.tolist():
        root = find(int(index))

        component_members.setdefault(
            root,
            [],
        ).append(int(index))

    components = [
        np.asarray(
            sorted(members),
            dtype=np.int64,
        )
        for members in component_members.values()
    ]

    components.sort(
        key=lambda values: int(values[0])
    )

    flattened = bri.reshape(
        chain_count,
        -1,
    )

    hit_rows: list[int] = []
    hit_distances: list[int] = []

    non_candidate_hits = 0
    tree_count = 0

    for members in components:
        if members.shape[0] < 2:
            raise FullBRINNSearchError(
                "Brain graph produced a singleton component"
            )

        component_points = flattened[
            members
        ]

        # A compressed cover tree is a metric tree over distinct
        # points. Collapse exact BRI duplicates and expand them
        # back to chain identities after search.
        unique_points, inverse = np.unique(
            component_points,
            axis=0,
            return_inverse=True,
        )

        unique_points = np.asarray(
            unique_points,
            dtype=np.int64,
        )
        inverse = np.asarray(
            inverse,
            dtype=np.int64,
        )

        unique_members = [
            members[
                np.flatnonzero(
                    inverse == unique_index
                )
            ]
            for unique_index in range(
                unique_points.shape[0]
            )
        ]

        tree = CompressedCoverTree(
            unique_points
        )
        tree_count += 1

        for unique_query in range(
            unique_points.shape[0]
        ):
            result = tree.radius_neighbors(
                unique_points[unique_query],
                radius_mA,
            )

            for (
                unique_subject,
                distance_mA,
            ) in zip(
                result.indices.tolist(),
                result.distances_mA.tolist(),
            ):
                unique_subject = int(
                    unique_subject
                )
                distance_mA = int(
                    distance_mA
                )

                # Each unordered unique-point relation once.
                if unique_subject < unique_query:
                    continue

                left_members = unique_members[
                    unique_query
                ]
                right_members = unique_members[
                    unique_subject
                ]

                if unique_subject == unique_query:
                    if left_members.shape[0] < 2:
                        continue

                    local_q, local_s = np.triu_indices(
                        left_members.shape[0],
                        k=1,
                    )

                    chain_pairs = (
                        (
                            int(left_members[a]),
                            int(left_members[b]),
                        )
                        for a, b in zip(
                            local_q.tolist(),
                            local_s.tolist(),
                        )
                    )
                else:
                    chain_pairs = (
                        (
                            int(left),
                            int(right),
                        )
                        for left in left_members.tolist()
                        for right in right_members.tolist()
                    )

                for left, right in chain_pairs:
                    if left > right:
                        left, right = right, left

                    row = edge_to_row.get(
                        (left, right)
                    )

                    if row is None:
                        non_candidate_hits += 1
                        continue

                    hit_rows.append(
                        int(row)
                    )
                    hit_distances.append(
                        distance_mA
                    )

    if not hit_rows:
        return ComponentRadiusResult(
            row_indices=np.empty(
                0,
                dtype=np.int64,
            ),
            distances_mA=np.empty(
                0,
                dtype=np.int64,
            ),
            participating_chain_count=(
                int(participating.shape[0])
            ),
            component_count=len(components),
            compressed_tree_count=tree_count,
            non_candidate_radius_hit_count=(
                non_candidate_hits
            ),
        )

    row_array = np.asarray(
        hit_rows,
        dtype=np.int64,
    )
    distance_array = np.asarray(
        hit_distances,
        dtype=np.int64,
    )

    order = np.argsort(
        row_array,
        kind="stable",
    )

    row_array = row_array[order]
    distance_array = distance_array[order]

    if (
        np.unique(row_array).shape[0]
        != row_array.shape[0]
    ):
        raise FullBRINNSearchError(
            "Stage-7 candidate row emitted more than once"
        )

    if np.any(
        distance_array > radius_mA
    ):
        raise FullBRINNSearchError(
            "Component radius search emitted "
            "an out-of-radius pair"
        )

    return ComponentRadiusResult(
        row_indices=row_array,
        distances_mA=distance_array,
        participating_chain_count=(
            int(participating.shape[0])
        ),
        component_count=len(components),
        compressed_tree_count=tree_count,
        non_candidate_radius_hit_count=(
            non_candidate_hits
        ),
    )
