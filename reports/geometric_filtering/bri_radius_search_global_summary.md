# PDB707K BRI Geometric Near-Duplicate Search

## Method

- BRI version: `1.2.2.2-comp702`
- Distance metric: L-infinity (Chebyshev)
- Candidate radius: `1.0`
- Final inclusion criterion: **strictly `L_inf_invariant < 1.0`**
- Chains were compared only with chains having the same residue length.
- The cKDTree search was validated against the original BRI pairwise comparison implementation.

## Dataset coverage

| Measure | Count |
|---|---:|
| Total model-1 source records | 461,581 |
| Missing local mmCIF before extraction | 413 |
| Ambiguous duplicate source records omitted | 3 |
| BRI extraction attempts | 461,165 |
| Successful BRI vectors | 460,920 |
| BRI extraction failures | 245 |
| Total source records without a usable BRI vector | 661 |
| Comparable chains | 460,806 |
| Singleton length groups | 114 |
| Equal-length groups compared | 1,113 |
| Completed search tasks | 1,113 |

The 661 omissions consist of:

- 413 chains from PDB entries missing from the local mmCIF archive;
- 218 residue-count mismatches;
- 27 chains where BRI invariants were not generated;
- 3 ambiguous duplicate source records.

Detailed evidence and proposed solutions are available in the
[BRI omission investigation](omission_audit/omission_investigation.md).

## Pairwise search results

| Measure | Count |
|---|---:|
| Possible same-length pairs | 317,507,782 |
| cKDTree candidates at or below 1.0 | 5,264,356 |
| Strict near-duplicate pairs below 1.0 | 5,258,576 |
| Boundary pairs exactly equal to 1.0 removed | 5,780 |
| Exact zero-distance pairs | 13,941 |
| Near-duplicate pairs with equal BRI sequence | 4,098,362 |
| Near-duplicate pairs with different BRI sequence | 1,160,214 |
| Near-duplicate pairs with equal dataset sequence | 4,098,370 |

## Sequence definitions

- **BRI sequence:** the coordinate-derived amino-acid sequence reconstructed
  from the residues used to generate each BRI vector.
- **Dataset sequence:** the sequence recorded for the chain in the original
  PDB707K source table.

The equality categories have the following cross-tabulation:

| BRI sequences equal | Dataset sequences equal | Pair count |
|---|---|---:|
| Yes | Yes | 4,098,362 |
| Yes | No | 0 |
| No | Yes | 8 |
| No | No | 1,160,206 |
| **Total** |  | **5,258,576** |

The dataset-sequence-equal count is eight higher because eight pairs have
identical PDB707K source sequences but different coordinate-derived sequences.

All eight compare `8BFZ-A` with chains from `7Q4M`. The coordinate-derived
sequences differ at residue 22: the `7Q4M` chains contain glutamate (`E`), while
`8BFZ-A` contains glycine (`G`).

Detailed evidence is available in the
[sequence equality explanation](sequence_audit/sequence_equality_explanation.md).

## Output artifacts

- Consolidated pair table: `/users/sgsrehm1/scratch/bri_geometric/invariant_shards_bri1222/radius_search_lt1/all_near_duplicate_pairs_lt1.csv.gz`
- Global JSON summary: `/users/sgsrehm1/scratch/bri_geometric/invariant_shards_bri1222/radius_search_lt1/radius_search_global_summary.json`
- Global CSV summary: `/users/sgsrehm1/scratch/bri_geometric/invariant_shards_bri1222/radius_search_lt1/radius_search_global_summary.csv`
- Per-task outputs: `/users/sgsrehm1/scratch/bri_geometric/invariant_shards_bri1222/radius_search_lt1`

## Interpretation

The 5,258,576 results are pairwise geometric relationships, not the number of
chains that should be removed.

A chain may participate in many near-duplicate pairs, and connected components
may contain indirect relationships where not every pair satisfies the strict
threshold. Final representative selection must therefore account for component
structure rather than removing one chain for every reported pair.

The search currently represents the 460,920 chains with successful BRI vectors.
Recovering omitted chains may change the final pair and cluster counts.
