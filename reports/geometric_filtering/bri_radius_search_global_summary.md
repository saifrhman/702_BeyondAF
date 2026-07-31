# PDB707K BRI Geometric Near-Duplicate Search

## Method

- BRI version: `1.2.2.2-comp702`
- Distance metric: L-infinity (Chebyshev)
- Candidate radius: `1.0`
- Final inclusion criterion: **strictly `L_inf_invariant < 1.0`**
- Chains were compared only with chains having the same residue length.
- The cKDTree method was previously validated against the original BRI comparison implementation.

## Dataset coverage

| Measure | Count |
|---|---:|
| Total dataset chains attempted | 461,165 |
| Successful BRI vectors | 460,920 |
| Comparable chains | 460,806 |
| Singleton chains | 114 |
| Extraction failures | 245 |
| Equal-length groups compared | 1,113 |
| Completed search tasks | 1,113 |

## Pairwise search results

| Measure | Count |
|---|---:|
| Possible same-length pairs | 317,507,782 |
| cKDTree candidates at or below 1 Å | 5,264,356 |
| Strict near-duplicate pairs below 1 Å | 5,258,576 |
| Boundary pairs equal to 1 Å removed | 5,780 |
| Exact zero-distance pairs | 13,941 |
| Near-duplicate pairs with equal BRI sequence | 4,098,362 |
| Near-duplicate pairs with different BRI sequence | 1,160,214 |
| Near-duplicate pairs with equal dataset sequence | 4,098,370 |

## Output artifacts

- Consolidated pair table: `/users/sgsrehm1/scratch/bri_geometric/invariant_shards_bri1222/radius_search_lt1/all_near_duplicate_pairs_lt1.csv.gz`
- Global JSON summary: `/users/sgsrehm1/scratch/bri_geometric/invariant_shards_bri1222/radius_search_lt1/radius_search_global_summary.json`
- Global CSV summary: `/users/sgsrehm1/scratch/bri_geometric/invariant_shards_bri1222/radius_search_lt1/radius_search_global_summary.csv`
- Per-task pair files and summaries: `/users/sgsrehm1/scratch/bri_geometric/invariant_shards_bri1222/radius_search_lt1`

## Interpretation

The 5,258,576 results are pairwise geometric relationships, not the number of chains that should be removed. Connected-component clustering is required to convert these pairwise links into near-duplicate groups and select retained representatives.
