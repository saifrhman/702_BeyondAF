import numpy as np

from pdbclean.full_bri_nn import (
    search_all_pairs_radius,
    search_brain_candidate_pairs,
)


def _distance(
    bri: np.ndarray,
    a: int,
    b: int,
) -> int:
    return int(
        np.max(
            np.abs(
                bri[a] - bri[b]
            )
        )
    )


def test_candidate_relation_preserved_with_duplicate_subject_vectors():
    bri = np.zeros(
        (6, 2, 9),
        dtype=np.int64,
    )

    bri[:, 0, 0] = np.array(
        [0, 5, 11, 20, 5, 100],
        dtype=np.int64,
    )

    query = np.array(
        [0, 0, 0, 1, 1, 4],
        dtype=np.int64,
    )
    subject = np.array(
        [1, 2, 4, 2, 3, 2],
        dtype=np.int64,
    )

    result = search_brain_candidate_pairs(
        bri,
        query,
        subject,
        radius_mA=10,
    )

    assert result.row_indices.tolist() == [
        0,
        2,
        3,
        5,
    ]
    assert result.distances_mA.tolist() == [
        5,
        5,
        6,
        6,
    ]


def test_candidate_search_matches_brute_force_randomized():
    for seed in range(50):
        rng = np.random.default_rng(
            20000 + seed
        )

        n = int(
            rng.integers(
                5,
                35,
            )
        )
        m = int(
            rng.integers(
                2,
                7,
            )
        )

        bri = rng.integers(
            -100,
            101,
            size=(n, m, 9),
            dtype=np.int64,
        )

        all_q, all_s = np.triu_indices(
            n,
            k=1,
        )

        chosen = rng.choice(
            all_q.shape[0],
            size=min(
                all_q.shape[0],
                int(
                    rng.integers(
                        1,
                        min(
                            120,
                            all_q.shape[0],
                        )
                        + 1,
                    )
                ),
            ),
            replace=False,
        )

        q = all_q[
            chosen
        ].astype(
            np.int64,
            copy=False,
        )
        s = all_s[
            chosen
        ].astype(
            np.int64,
            copy=False,
        )

        radius = int(
            rng.integers(
                0,
                101,
            )
        )

        result = search_brain_candidate_pairs(
            bri,
            q,
            s,
            radius_mA=radius,
        )

        expected_rows = []
        expected_distances = []

        for row in range(
            q.shape[0]
        ):
            distance = _distance(
                bri,
                int(q[row]),
                int(s[row]),
            )

            if distance <= radius:
                expected_rows.append(
                    row
                )
                expected_distances.append(
                    distance
                )

        assert (
            result.row_indices.tolist()
            == expected_rows
        )
        assert (
            result.distances_mA.tolist()
            == expected_distances
        )


def test_all_pairs_radius_matches_brute_force_randomized():
    for seed in range(50):
        rng = np.random.default_rng(
            30000 + seed
        )

        n = int(
            rng.integers(
                2,
                35,
            )
        )
        m = int(
            rng.integers(
                1,
                5,
            )
        )

        bri = rng.integers(
            -50,
            51,
            size=(n, m, 9),
            dtype=np.int64,
        )

        # Deliberately create exact duplicate complete-BRI
        # vectors in some trials.
        if n >= 4 and seed % 3 == 0:
            bri[2] = bri[0]
            bri[3] = bri[0]

        radius = int(
            rng.integers(
                0,
                61,
            )
        )

        result = search_all_pairs_radius(
            bri,
            radius_mA=radius,
        )

        expected = []

        for q in range(n):
            for s in range(
                q + 1,
                n,
            ):
                distance = _distance(
                    bri,
                    q,
                    s,
                )

                if distance <= radius:
                    expected.append(
                        (
                            q,
                            s,
                            distance,
                        )
                    )

        observed = list(
            zip(
                result.query_indices.tolist(),
                result.subject_indices.tolist(),
                result.distances_mA.tolist(),
            )
        )

        assert observed == expected


def test_component_search_preserves_candidate_edge_relation():
    from pdbclean.full_bri_nn import (
        search_brain_candidate_components,
    )

    bri = np.zeros(
        (3, 2, 9),
        dtype=np.int64,
    )

    # All three chains are within the complete-BRI radius:
    # d(0,1)=5, d(1,2)=5, d(0,2)=10.
    #
    # But only (0,1) and (1,2) are supplied as Brain edges.
    # The component search must therefore not emit (0,2).
    bri[:, 0, 0] = np.array(
        [0, 5, 10],
        dtype=np.int64,
    )

    q = np.array(
        [0, 1],
        dtype=np.int64,
    )
    s = np.array(
        [1, 2],
        dtype=np.int64,
    )

    result = search_brain_candidate_components(
        bri,
        q,
        s,
        radius_mA=10,
    )

    assert result.row_indices.tolist() == [0, 1]
    assert result.distances_mA.tolist() == [5, 5]
    assert result.participating_chain_count == 3
    assert result.component_count == 1
    assert result.compressed_tree_count == 1

    # (0,2) was encountered geometrically but rejected because
    # it was not in the supplied Brain candidate relation.
    assert result.non_candidate_radius_hit_count == 1


def test_component_search_matches_candidate_bruteforce_randomized():
    from pdbclean.full_bri_nn import (
        search_brain_candidate_components,
    )

    for seed in range(75):
        rng = np.random.default_rng(
            40000 + seed
        )

        n = int(
            rng.integers(
                5,
                40,
            )
        )
        m = int(
            rng.integers(
                2,
                7,
            )
        )

        bri = rng.integers(
            -100,
            101,
            size=(n, m, 9),
            dtype=np.int64,
        )

        # Deliberately introduce exact complete-BRI duplicates.
        if n >= 6 and seed % 4 == 0:
            bri[4] = bri[1]
            bri[5] = bri[1]

        all_q, all_s = np.triu_indices(
            n,
            k=1,
        )

        edge_count = int(
            rng.integers(
                1,
                min(
                    all_q.shape[0],
                    160,
                )
                + 1,
            )
        )

        chosen = rng.choice(
            all_q.shape[0],
            size=edge_count,
            replace=False,
        )

        q = all_q[
            chosen
        ].astype(
            np.int64,
            copy=False,
        )
        s = all_s[
            chosen
        ].astype(
            np.int64,
            copy=False,
        )

        radius = int(
            rng.integers(
                0,
                101,
            )
        )

        result = search_brain_candidate_components(
            bri,
            q,
            s,
            radius_mA=radius,
        )

        expected_rows = []
        expected_distances = []

        for row in range(
            q.shape[0]
        ):
            distance = _distance(
                bri,
                int(q[row]),
                int(s[row]),
            )

            if distance <= radius:
                expected_rows.append(row)
                expected_distances.append(
                    distance
                )

        assert (
            result.row_indices.tolist()
            == expected_rows
        )
        assert (
            result.distances_mA.tolist()
            == expected_distances
        )
