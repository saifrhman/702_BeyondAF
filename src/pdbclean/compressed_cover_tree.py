"""Exact k-nearest-neighbour search with a compressed cover tree.

Implements the compressed-cover-tree construction and kNN traversal of
Elkin & Kurlin (ICML 2023) for the metric needed by COMP702.

COMP702 uses exact integer milliangstrom complete-BRI vectors with the
L-infinity metric.

Important implementation note
-----------------------------
The published definition of Next(p, i, T) contains an apparent notation
typo in the prose condition.  The linked-child-level description and the
algorithms require Next(p, i, T) to mean the greatest child level j < i.
That operational definition is used here.

Reference:
Y. Elkin and V. Kurlin,
"A New Near-linear Time Algorithm For k-Nearest Neighbor Search Using a
Compressed Cover Tree", ICML 2023.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


class CompressedCoverTreeError(RuntimeError):
    """Raised when compressed-cover-tree invariants cannot be satisfied."""


def linf_distance_mA(a: np.ndarray, b: np.ndarray) -> int:
    """Exact L-infinity distance between integer-mA vectors."""
    a = np.asarray(a, dtype=np.int64)
    b = np.asarray(b, dtype=np.int64)

    if a.shape != b.shape:
        raise ValueError(
            f"L-infinity shape mismatch: {a.shape} != {b.shape}"
        )

    return int(np.max(np.abs(a - b)))


def _level_for_positive_distance(distance: int) -> int:
    """Return maximal x with distance > 2**x.

    Equivalently ceil(log2(distance)) - 1, evaluated exactly for integers.
    """
    distance = int(distance)

    if distance <= 0:
        raise ValueError("Cover-tree reference points must be distinct")

    return (distance - 1).bit_length() - 1


@dataclass(frozen=True)
class KNNResult:
    indices: np.ndarray
    distances_mA: np.ndarray


class CompressedCoverTree:
    """Compressed cover tree over distinct integer metric points."""

    def __init__(self, points: np.ndarray):
        points = np.asarray(points, dtype=np.int64)

        if points.ndim != 2:
            raise ValueError(
                f"points must be a 2-D matrix, got shape {points.shape}"
            )

        if points.shape[0] < 1:
            raise ValueError("Compressed cover tree requires >=1 point")

        # A metric-space reference set cannot contain two distinct vertices
        # at distance zero.  Stage-8 publication will collapse identical BRI
        # vectors before building a tree and later expand their identities.
        if np.unique(points, axis=0).shape[0] != points.shape[0]:
            raise ValueError(
                "Reference points contain duplicate metric vectors; "
                "collapse identical BRI vectors before tree construction"
            )

        self.points = np.ascontiguousarray(points, dtype=np.int64)
        self.n_points = int(points.shape[0])
        self.dimension = int(points.shape[1])

        self.root = 0
        self.parent = np.full(self.n_points, -1, dtype=np.int64)
        self.level = np.zeros(self.n_points, dtype=np.int64)

        # Actual tree children only.  The paper's Child(p) additionally
        # includes p itself; algorithms below add the current reference
        # set explicitly, so self-links are unnecessary.
        self.children: list[list[int]] = [
            [] for _ in range(self.n_points)
        ]

        self._build()

        self.l_max = int(self.level[self.root])
        self.l_min = int(np.min(self.level))

        self._subtree_size = np.ones(
            self.n_points,
            dtype=np.int64,
        )
        self._compute_subtree_sizes(self.root)

    def _distance_indices(self, a: int, b: int) -> int:
        return linf_distance_mA(
            self.points[a],
            self.points[b],
        )

    def _distance_query(self, query: np.ndarray, index: int) -> int:
        return linf_distance_mA(
            query,
            self.points[index],
        )

    def _current_nonroot_levels(self) -> np.ndarray:
        return self.level[1:]

    def _current_l_min(self) -> int:
        return int(np.min(self._current_nonroot_levels()))

    def _current_l_max(self) -> int:
        return 1 + int(np.max(self._current_nonroot_levels()))

    def _children_at_level(
        self,
        node: int,
        level: int,
    ) -> list[int]:
        return [
            child
            for child in self.children[node]
            if int(self.level[child]) == int(level)
        ]

    def _next_level(
        self,
        node: int,
        current_level: int,
        l_min: int,
    ) -> int:
        """Next(p,i,T): greatest Child(p) level strictly below i.

        Elkin--Kurlin define Child(p) to contain p itself in addition
        to its ordinary tree children.  Hence l(p) is also a candidate
        next level whenever l(p) < i.
        """
        candidates: list[int] = []

        node_level = int(self.level[node])

        # Child(p) includes p itself.
        if node_level < int(current_level):
            candidates.append(node_level)

        candidates.extend(
            int(self.level[child])
            for child in self.children[node]
            if int(self.level[child]) < int(current_level)
        )

        if not candidates:
            return int(l_min) - 1

        return max(candidates)

    def _attach(
        self,
        child: int,
        parent: int,
        level: int,
    ) -> None:
        self.parent[child] = int(parent)
        self.level[child] = int(level)
        self.children[parent].append(int(child))

    def _build(self) -> None:
        if self.n_points == 1:
            self.level[self.root] = 0
            return

        # Algorithm 3.5's one-root initial state has no finite child level.
        # The first non-root point therefore attaches directly to the root.
        first = 1
        distance = self._distance_indices(
            self.root,
            first,
        )
        first_level = _level_for_positive_distance(
            distance
        )
        self._attach(
            first,
            self.root,
            first_level,
        )

        for point in range(2, self.n_points):
            self._add_point(point)

        # Algorithm 3.4 final root-level assignment.
        self.level[self.root] = (
            1 + int(np.max(self.level[1:]))
        )

    def _add_point(self, point: int) -> None:
        """Algorithm 3.5: add one point to an existing tree."""
        # During incremental construction the paper keeps the root at
        # l(r)=+infinity.  We realize that symbolically infinite level by
        # choosing, for this insertion, the smallest finite top level that
        # is both compatible with the existing tree and strictly above the
        # level required to cover the new point from the root.
        root_distance = self._distance_indices(
            self.root,
            point,
        )
        required_root_level = (
            _level_for_positive_distance(root_distance)
            + 1
        )

        l_max = max(
            self._current_l_max(),
            required_root_level,
        )
        l_min = self._current_l_min()

        i = l_max - 1

        # R_level stores the R_i sets occurring during this insertion.
        reference_sets: dict[int, list[int]] = {
            l_max: [self.root],
        }

        # eta(i) is the next higher iteration level.
        eta: dict[int, int] = {
            i: l_max,
        }

        while i >= l_min:
            upper = eta[i]
            upper_set = reference_sets[upper]

            candidates: set[int] = set(upper_set)

            for node in upper_set:
                candidates.update(
                    self._children_at_level(
                        node,
                        i,
                    )
                )

            radius = 1 << (i + 1)

            r_i = [
                node
                for node in sorted(candidates)
                if self._distance_query(
                    self.points[point],
                    node,
                )
                <= radius
            ]

            reference_sets[i] = r_i

            if not r_i:
                break

            next_i = max(
                self._next_level(
                    node,
                    i,
                    l_min,
                )
                for node in r_i
            )

            # The published Algorithm 3.5 line is evidently meant to
            # maintain eta(next_i)=i, exactly as Algorithm 4.3 does for
            # the kNN traversal.
            eta[next_i] = i
            i = next_i

        # Elkin--Kurlin AssignParent subprocedure.
        #
        # M is the dictionary of R_i sets saved during AddPoint.
        # Starting from its lowest key, find the first level i for
        # which d(p, R_i) <= 2**i.  Then choose a nearest q in R_i
        # and set l(p) to the maximal x satisfying d(p,q) > 2**x.
        #
        # This guarantees l(p) < l(q) together with the covering
        # condition required by the compressed cover tree.
        for parent_level in sorted(reference_sets):
            parent_set = reference_sets[parent_level]

            if not parent_set:
                continue

            nearest_parent = min(
                parent_set,
                key=lambda node: (
                    self._distance_query(
                        self.points[point],
                        node,
                    ),
                    node,
                ),
            )

            distance = self._distance_indices(
                point,
                nearest_parent,
            )

            # Compare integer-milliangstrom distance with 2**level
            # without using a negative bit shift.  For level < 0,
            # 2**level < 1, so a non-negative integer distance can
            # satisfy the bound only when it is exactly zero.
            if (
                distance == 0
                if parent_level < 0
                else distance <= (1 << parent_level)
            ):
                point_level = _level_for_positive_distance(
                    distance
                )

                self._attach(
                    point,
                    nearest_parent,
                    point_level,
                )
                return

        raise CompressedCoverTreeError(
            f"AssignParent failed for point {point}"
        )

    def _compute_subtree_sizes(
        self,
        node: int,
    ) -> int:
        size = 1

        for child in self.children[node]:
            size += self._compute_subtree_sizes(
                child
            )

        self._subtree_size[node] = size
        return size

    def distinctive_descendant_size(
        self,
        node: int,
        level: int,
    ) -> int:
        """|S_i(node,T)| from the paper.

        Descendant subtrees rooted at children with child-level >= i are
        excluded; children below i remain wholly distinctive.
        """
        size = 1

        for child in self.children[node]:
            if int(self.level[child]) < int(level):
                size += int(
                    self._subtree_size[child]
                )

        return size

    def _collect_subtree(
        self,
        node: int,
        output: list[int],
    ) -> None:
        output.append(int(node))

        for child in self.children[node]:
            self._collect_subtree(
                child,
                output,
            )

    def distinctive_descendants(
        self,
        node: int,
        level: int,
    ) -> list[int]:
        output = [int(node)]

        for child in self.children[node]:
            if int(self.level[child]) < int(level):
                self._collect_subtree(
                    child,
                    output,
                )

        return output

    def _lambda_point(
        self,
        query: np.ndarray,
        candidates: list[int],
        level: int,
        k: int,
    ) -> tuple[int, int]:
        """Algorithm D.8 / lambda_k(q,C)."""
        ranked = sorted(
            (
                self._distance_query(
                    query,
                    node,
                ),
                int(node),
            )
            for node in candidates
        )

        if not ranked:
            raise CompressedCoverTreeError(
                "lambda-point candidate set is empty"
            )

        cumulative = 0
        cursor = 0

        while cursor < len(ranked):
            distance = ranked[cursor][0]
            tied_nodes: list[int] = []

            while (
                cursor < len(ranked)
                and ranked[cursor][0] == distance
            ):
                tied_nodes.append(
                    ranked[cursor][1]
                )
                cursor += 1

            for node in tied_nodes:
                cumulative += (
                    self.distinctive_descendant_size(
                        node,
                        level,
                    )
                )

            if cumulative >= k:
                # Any point at this distance satisfies the lambda-point
                # definition.  Use the smallest index deterministically.
                return min(tied_nodes), int(distance)

        raise CompressedCoverTreeError(
            "lambda-point cumulative descendant count "
            f"{cumulative} is smaller than k={k}"
        )

    def _top_k(
        self,
        query: np.ndarray,
        candidates: list[int] | set[int],
        k: int,
    ) -> KNNResult:
        unique_candidates = sorted(
            set(int(x) for x in candidates)
        )

        ranked = sorted(
            (
                self._distance_query(
                    query,
                    node,
                ),
                node,
            )
            for node in unique_candidates
        )

        if len(ranked) < k:
            raise CompressedCoverTreeError(
                f"Only {len(ranked)} final candidates for k={k}"
            )

        selected = ranked[:k]

        return KNNResult(
            indices=np.asarray(
                [node for _, node in selected],
                dtype=np.int64,
            ),
            distances_mA=np.asarray(
                [distance for distance, _ in selected],
                dtype=np.int64,
            ),
        )

    def knn(
        self,
        query: np.ndarray,
        k: int,
    ) -> KNNResult:
        """Algorithm 4.3: exact k-nearest neighbours."""
        query = np.asarray(
            query,
            dtype=np.int64,
        )

        if query.shape != (self.dimension,):
            raise ValueError(
                f"query shape must be {(self.dimension,)}, "
                f"got {query.shape}"
            )

        k = int(k)

        if not 1 <= k <= self.n_points:
            raise ValueError(
                f"k must be in [1,{self.n_points}], got {k}"
            )

        if self.n_points == 1:
            return self._top_k(
                query,
                [self.root],
                1,
            )

        i = self.l_max - 1

        eta: dict[int, int] = {
            i: self.l_max,
        }

        reference_sets: dict[int, list[int]] = {
            self.l_max: [self.root],
        }

        while i >= self.l_min:
            upper = eta[i]
            upper_set = reference_sets[upper]

            candidates: set[int] = set(
                upper_set
            )

            for node in upper_set:
                candidates.update(
                    self._children_at_level(
                        node,
                        i,
                    )
                )

            candidates_list = sorted(
                candidates
            )

            _, lambda_distance = (
                self._lambda_point(
                    query,
                    candidates_list,
                    i,
                    k,
                )
            )

            prune_bound = (
                lambda_distance
                + (1 << (i + 2))
            )

            r_i = [
                node
                for node in candidates_list
                if self._distance_query(
                    query,
                    node,
                )
                <= prune_bound
            ]

            reference_sets[i] = r_i

            if (
                lambda_distance
                > (1 << (i + 2))
            ):
                final_candidates: set[int] = set()

                for node in r_i:
                    final_candidates.update(
                        self.distinctive_descendants(
                            node,
                            i,
                        )
                    )

                return self._top_k(
                    query,
                    final_candidates,
                    k,
                )

            next_i = max(
                self._next_level(
                    node,
                    i,
                    self.l_min,
                )
                for node in r_i
            )

            eta[next_i] = i
            i = next_i

        # The displayed pseudocode names R_lmin in its final line,
        # while the accompanying complexity proof states explicitly
        # that the final selection is made from R_{eta(i)} after the
        # compressed-level traversal exits.  This distinction matters
        # when Next(...) jumps directly below l_min because no explicit
        # R_lmin iteration is materialized.
        #
        # At loop exit, i < l_min and eta(i) points to the final
        # materialized reference level.
        if i not in eta:
            raise CompressedCoverTreeError(
                "kNN traversal lost final eta(i)"
            )

        final_level = eta[i]

        if final_level not in reference_sets:
            raise CompressedCoverTreeError(
                "kNN traversal lost final reference set"
            )

        return self._top_k(
            query,
            reference_sets[final_level],
            k,
        )

    def radius_neighbors(
        self,
        query: np.ndarray,
        radius_mA: int,
    ) -> KNNResult:
        """Retrieve every reference point within a closed radius.

        The cited method is kNN rather than a radius-query algorithm.
        We therefore run exact kNN with geometrically increasing k until
        the kth distance lies outside the requested radius.  At that point
        every reference point inside the closed radius is necessarily among
        the returned exact nearest neighbours.
        """
        radius_mA = int(radius_mA)

        if radius_mA < 0:
            raise ValueError(
                "radius_mA must be non-negative"
            )

        k = min(8, self.n_points)

        while True:
            result = self.knn(
                query,
                k,
            )

            inside = (
                result.distances_mA
                <= radius_mA
            )

            if (
                result.distances_mA[-1]
                > radius_mA
                or k == self.n_points
            ):
                return KNNResult(
                    indices=result.indices[
                        inside
                    ],
                    distances_mA=(
                        result.distances_mA[
                            inside
                        ]
                    ),
                )

            k = min(
                self.n_points,
                2 * k,
            )

    def validate(self) -> None:
        """Exhaustively validate compressed-cover-tree axioms.

        Intended for unit tests and small validation buckets.
        """
        if self.parent[self.root] != -1:
            raise CompressedCoverTreeError(
                "Root unexpectedly has a parent"
            )

        for child in range(self.n_points):
            if child == self.root:
                continue

            parent = int(
                self.parent[child]
            )

            if parent < 0:
                raise CompressedCoverTreeError(
                    f"Node {child} has no parent"
                )

            if not (
                int(self.level[child])
                < int(self.level[parent])
            ):
                raise CompressedCoverTreeError(
                    f"Level ordering failed for node {child}"
                )

            distance = self._distance_indices(
                child,
                parent,
            )

            covering_radius = (
                1
                << (
                    int(self.level[child])
                    + 1
                )
            )

            if distance > covering_radius:
                raise CompressedCoverTreeError(
                    f"Covering condition failed for node {child}"
                )

        # Separation condition.  With integer-mA distinct vectors and the
        # construction above, l_min cannot be below -1.
        for i in range(
            self.l_min,
            self.l_max + 1,
        ):
            cover = [
                node
                for node in range(
                    self.n_points
                )
                if int(self.level[node]) >= i
            ]

            if i < 0:
                # All distinct integer points have distance >=1 > 2**i.
                continue

            threshold = 1 << i

            for pos, a in enumerate(
                cover
            ):
                for b in cover[
                    pos + 1:
                ]:
                    if (
                        self._distance_indices(
                            a,
                            b,
                        )
                        <= threshold
                    ):
                        raise CompressedCoverTreeError(
                            "Separation condition failed "
                            f"at level {i}: {a}, {b}"
                        )
