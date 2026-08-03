# Sequence Interpretation of BRI Near-Duplicate Pairs

## Scope

The strict geometric search identified 5,258,576 chain pairs satisfying:

- equal chain length; and
- strict backbone rigid invariant distance `L∞ < 1`.

Each pair was assessed using two sequence representations:

1. **BRI sequence**: the amino-acid sequence reconstructed from residues present
   in the mmCIF coordinates used to generate the BRI vector.
2. **Dataset sequence**: the sequence stored in the original PDB707K source
   dataset.

## Equality cross-tabulation

| BRI sequences equal | Dataset sequences equal | Pair count |
|---|---|---:|
| Yes | Yes | 4,098,362 |
| Yes | No | 0 |
| No | Yes | 8 |
| No | No | 1,160,206 |
| **Total** |  | **5,258,576** |

## Interpretation of the reported counts

### Near-duplicate pairs with equal BRI sequence

There are **4,098,362** geometrically near-duplicate pairs whose
coordinate-derived amino-acid sequences are identical.

These pairs are structurally close under the BRI threshold and contain the same
residue identities in the coordinates used for BRI extraction.

### Near-duplicate pairs with different BRI sequence

There are **1,160,214** geometrically near-duplicate pairs whose
coordinate-derived amino-acid sequences differ:

- 1,160,206 also have different source-dataset sequences;
- 8 have equal source-dataset sequences but different coordinate-derived
  sequences.

These pairs show that geometric near-duplication does not require complete
sequence identity.

### Near-duplicate pairs with equal dataset sequence

There are **4,098,370** geometrically near-duplicate pairs whose original
PDB707K source sequences are identical.

This total is eight larger than the BRI-sequence-equal total:

`4,098,370 - 4,098,362 = 8`

## Explanation of the eight-pair discrepancy

All eight discrepant pairs have:

- equal source-dataset sequences;
- different coordinate-derived BRI sequences;
- one chain classified as `sequence_match`;
- the other chain classified as `minor_substitution`.

Therefore, the source dataset records the two chains as sequence-identical, but
the residues observed in their current mmCIF coordinates differ by at least one
minor substitution.

No pair showed the opposite pattern: there were zero cases where the
coordinate-derived BRI sequences were equal but the source-dataset sequences
were different.

## Conclusion

The BRI and dataset sequence counts measure related but distinct properties:

- BRI sequence equality describes the residues actually represented in the
  coordinates used for geometric comparison.
- Dataset sequence equality describes the sequences recorded in the source
  PDB707K table.

The eight-pair difference is a source-to-coordinate sequence discrepancy, not
an error in the pair-counting calculation.

## Exact source of the eight-pair discrepancy

All eight discrepant pairs compare `8BFZ-A` with one of the following `7Q4M`
chains:

- `7Q4M-B`
- `7Q4M-D`
- `7Q4M-E`
- `7Q4M-F`
- `7Q4M-G`
- `7Q4M-H`
- `7Q4M-I`
- `7Q4M-J`

Every chain has length 31 and starts at residue 12.

The coordinate-derived sequences differ at sequence position 11, corresponding
to residue ID 22:

- the `7Q4M` chains contain `E` (glutamate);
- `8BFZ-A` contains `G` (glycine).

Thus, all eight count discrepancies arise from the same coordinate-level
`E22G` substitution in `8BFZ-A`. The original PDB707K dataset sequences remain
identical for these pairs.

Their geometric distances range from `L∞ = 0.510` to `0.631`, so all eight
still satisfy the strict geometric near-duplicate criterion of `L∞ < 1`.
