# BRI Duplicate Detection — Analysis Handover

Date: 2026-08-18

## Status

The computational BRI duplicate-detection pipeline is complete.

Manual scientific interpretation of individual BRI duplicate and
near-duplicate cases is intentionally stopped at this point and is
handed over to another member of the research group.

No further case-by-case BRI analysis should be required for the
COMP702 pipeline owner before proceeding to the next project stage.

---

## Canonical dataset

PDB snapshot:

- 2026-01-01

Structural model used:

- model 1

Post-cleaning and geometry-valid population:

- 578,524 chains

---

## Paper-faithful geometric search

Implemented sequence:

1. cleaned protein chains
2. geometry validation
3. complete BRI
4. 0.001 Å BRI representation
5. Brain
6. exact retained-chain-length buckets
7. Brain filtering
8. fast nearest-neighbour search on complete BRI
9. exact complete-BRI L-infinity distance
10. near-duplicate threshold <= 0.010 Å

The fast complete-BRI search was validated against the exhaustive
oracle.

### Raw geometric results

Candidate comparisons tested:

- 3,531,895

Exact BRI duplicate chain pairs, dBRI = 0:

- 17,373

All near-duplicate chain pairs, dBRI <= 0.010 Å:

- 1,072,751

Nonzero near-duplicate chain pairs:

- 1,055,378

Pairs above 0.010 Å:

- 2,459,144

These counts describe chain-pair geometric evidence and are not
equivalent to numbers of redundant PDB depositions.

---

## Paper-derived detailed-inspection filtering

Starting geometric near-duplicate chain pairs:

- 1,072,751

After X-ray crystallographic scope and resolution <= 4 Å:

- 95,297

After PanDDA exclusion:

- 95,285

After same-deposition exclusion:

- 11,572 chain pairs

Unique unordered deposition pairs:

- 2,537

Participating PDB depositions:

- 629

Publication:

outputs/pdbclean/20260101/protocol3.2-comp702-v1/\
acta_downstream_investigation_v2/

Manual-review manifest:

outputs/pdbclean/20260101/protocol3.2-comp702-v1/\
acta_manual_review_manifest_v2/

---

## Manual complexity review performed so far

A deposition-level manual review workspace was created to reproduce the
paper's handling of virus and ribosome structures.

Current deposition classifications:

- ribosome_complex: 72
- virus_complex: 32
- neither_complexity_case: 61
- unreviewed: 464

Pairs excluded because at least one deposition was manually classified
as a virus/ribosome complexity case:

- 683

Pairs remaining for detailed review:

- 1,854

These manual classifications are working analysis records only.
They have not modified the frozen geometric publications.

Working files:

- reports/acta_2026_deposition_review.csv
- reports/acta_2026_manual_review.csv
- reports/acta_2026_detailed_review.csv

---

## Exact-BRI exploratory review

Within the 1,854-pair detailed-review population:

Pairs containing exact BRI evidence:

- 101

Pairs containing only nonzero near-duplicate evidence:

- 1,753

Three exact-BRI cases were identified as trivial matches involving only
a retained one-residue chain:

- 3g3p / 8vjw
- 3qxy / 4z33
- 4zkr / 5nwk

These were recorded as:

- trivial_single_residue_match

Substantive exact-BRI deposition pairs remaining:

- 98

For those 98 pairs:

- 119 unique coordinate mmCIF depositions were resolved from the
  frozen 2026 snapshot
- 152 exact matched-chain comparisons were examined
- raw-coordinate signatures equal: 121
- B-factor signatures equal: 118
- occupancy signatures equal: 123

Evidence files:

- reports/acta_2026_exact_bri_deposition_evidence.csv
- reports/acta_2026_exact_bri_pair_evidence.csv

Exploratory evidence grouping:

- all matched raw coordinates and B factors equal: 80 pairs
- all raw coordinates equal but B factors differ: 1 pair
- partial coordinate reuse: 7 pairs
- no raw Cartesian-coordinate reuse: 10 pairs

These categories are diagnostic evidence groupings only.
They are NOT duplicate/non-duplicate scientific classifications.

Of the 98 substantive exact-BRI pairs:

- 43 occur in Acta Table 1
- 55 do not occur in Acta Table 1

The phrase "not in Acta Table 1" must be used instead of "2026-only",
because many of these entries predate the paper snapshot.

The 55 pairs not in Acta Table 1 collapse into 15 connected review
groups. No scientific interpretation of these groups was completed.

---

## Important scientific boundary

BRI geometric similarity alone does not determine whether two PDB
depositions should be treated as redundant, legitimate independent
experiments, reused coordinates, related ligand experiments, or other
scientific cases.

The original Acta work performs manual crystallographic investigation
after candidate discovery.

That manual interpretation is outside the remaining COMP702 work for
the current pipeline owner and is handed over to the research group.

No representative-selection rule has been defined.

No structure should be deleted merely because:

- dBRI = 0;
- dBRI <= 0.010 Å;
- coordinates happen to be identical;
- it belongs to a connected component of the duplicate graph.

---

## Frozen scientific outputs

The following computational results should remain unchanged:

- cleaning
- geometry validation
- complete BRI
- Brain
- length bucketing
- Brain filtering
- complete-BRI nearest-neighbour search
- exhaustive validation oracle
- geometric near-duplicate classification
- downstream metadata
- Acta detailed-inspection candidate filtering
- Acta manual-review manifest generation

The manual-analysis reports under `reports/` are downstream working
records and do not alter these publications.

---

## Handover point

Research-group follow-up may continue with:

- manual inspection of the 1,854 remaining deposition pairs;
- interpretation of exact and nonzero BRI cases;
- examination of crystallographic evidence;
- comparison with the Acta paper;
- classification of duplication mechanisms;
- preparation of figures/tables/publication analysis.

The COMP702 pipeline owner should proceed to the next independent
project stage rather than continue BRI case analysis.

