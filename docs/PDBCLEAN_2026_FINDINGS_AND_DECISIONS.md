# PDBClean 2026 — Final Findings and Scientific Decisions

## Frozen scope

Snapshot: **2026-01-01**
Model scope: **model 1 only**
Final Stage-14 release commit: `6b15ff02f8cf5caf8b5a3021ec64d928d145d2a0`

Stages 1–14 are complete and frozen.

This work deliberately distinguishes:

1. MATCH mathematical definitions.
2. Acta paper methodology.
3. Reference-software behaviour where required for validation.
4. COMP702 implementation and downstream policy choices.

The papers are the scientific authority. Software-package behaviour is not
silently substituted for paper methodology.

No old-snapshot comparison is used to make 2026 scientific decisions.

---

## Canonical pipeline

PDB 2026-01-01, model 1

-> Protocol 3.2-derived cleaning
-> post-cleaning geometric validation
-> complete BRI
-> exact 0.001 A representation
-> Brain
-> exact chain-length buckets
-> Brain filtering
-> fast nearest-neighbour search on complete BRI
-> exact complete-BRI L-infinity distance
-> dBRI <= 0.010 A classification
-> Acta-style detailed investigation
-> COMP702 geometric representative selection
-> final deduplicated release

---

## Stage 1 — Cleaning

Accepted after cleaning:

**578,532 chains**

Cleaning follows Protocol 3.2-derived criteria with the Acta terminal-residue
exception.

Relevant exclusions include:

- non-proteins
- partial occupancy / disorder
- residue-index discontinuity
- missing N/CA/C atoms
- nonstandard amino acids

Incomplete terminal residues may be trimmed when allowed by the protocol.

---

## Stage 2 — Geometric validation

Required:

- six relevant backbone distances >= 0.01 A
- N-CA-C internal angles >= 3 degrees

Results:

- cleaning accepted: 578,532
- geometry-invalid quarantined: 8
- canonical eligible population: **578,524**
- angle violations: 13
- distance violations: 0

Canonical post-geometry source:

`geometric_validation/finalized/eligible.parquet`

---

## Stage 3 — Complete BRI

The paper-defined complete Backbone Rigid Invariant was computed for all
canonical chains.

Canonical BRI population:

**578,524**

---

## Stage 4 — Exact precision representation

Coordinates are represented at 0.001 A precision.

The operational exact representation is:

`BRI_mA = round(1000 * BRI)`

Complete-BRI distance is the L-infinity maximum over the complete invariant.

Important distinction:

- 0.001 A = representation precision
- 0.010 A = geometric near-duplicate threshold

Classification:

- `d_bri_mA = 0`: exact BRI duplicate
- `1 <= d_bri_mA <= 10`: nonzero near-duplicate
- `d_bri_mA > 10`: not near-duplicate

---

## Stage 5 — Brain

Brain follows MATCH Definition 5.1.

Results:

- m >= 2 Brain-defined chains: **577,760**
- m = 1 Brain-undefined chains: **764**

m=1 bypasses Brain.

---

## Stage 6 — Exact chain-length buckets

Only chains with exactly the same final cleaned length m are compared.

Results:

- total chains: 578,524
- distinct m values: 1,308
- m range: 1–4,629
- m=1: 764

---

## Stage 7 — Brain filtering

For same-length chains:

`dBrain = max_j |Brain_j(S) - Brain_j(Q)|`

Retain if:

`dBrain <= 0.010 A`

Production result:

**3,240,429 m>=2 Brain candidate pairs**

Validation:

- 40 buckets
- 86,702 decisions
- false negatives: 0
- false positives: 0

SciPy cKDTree is a COMP702 engineering implementation choice for executing the
paper-defined Brain filter.

---

## Stage 8 — Paper-faithful fast complete-BRI search

Scientific sequence:

exact chain length
-> Brain filtering
-> fast nearest-neighbour search on complete BRI
-> exact complete-BRI L-infinity distance

Production search used an exact compressed cover tree implementation for the
fast nearest-neighbour stage.

Results:

- m>=2 Brain candidates: 3,240,429
- m>=2 near-duplicate pairs: **1,068,256**
- m=1 near-duplicate pairs: **4,495**
- total geometric near-duplicate pairs: **1,072,751**

The fast result exactly matched the preserved exhaustive oracle.

---

## Stage 9 — Distance representation

Authoritative distance:

`d_bri_mA`

Near-duplicate boundary:

`d_bri_mA <= 10`

The boundary is inclusive.

---

## Stage 10 — Geometric classification

Pairs tested by exhaustive validated baseline:

**3,531,895**

Results:

- exact duplicates: **17,373**
- all near-duplicates: **1,072,751**
- nonzero near-duplicates: **1,055,378**
- non-near-duplicates: **2,459,144**

---

## Stage 11 — Acta-style downstream investigation

The full geometric result remains 1,072,751 chain pairs.

Detailed investigation is a separate downstream subset.

Paper-faithful / supplement-supported crystallographic sequence:

1. X-RAY DIFFRACTION
2. `refine_ls_d_res_high <= 4 A`
3. reject PanDDA group-deposition records
4. reject same-deposition hits
5. manual detailed inspection

Counts:

1,072,751
-> **95,297** X-ray <=4 A
-> **95,285** after PanDDA rejection
-> **11,572** cross-deposition chain pairs
-> **2,537** unique deposition pairs

The older EM-inclusive Stage-11 v1 interpretation is preserved only for
provenance and is scientifically superseded by v2.

---

## Stage 12 — Scientific validation gates

Validated properties include:

- exact complete-BRI metric
- Brain lower-bound behaviour
- exact 10 mA threshold boundary
- fast search == exhaustive oracle
- no use of connectedness as duplicate equivalence
- no old-snapshot comparison in 2026 decisions

All required gates passed.

---

## Stage 13 — Detailed investigation

Manual complexity review excluded 683 pair rows involving structures considered
impractical for detailed inspection, including ribosome/whole-virus complexes.

Detailed-review population:

**1,854 deposition pairs**

Evidence classes:

- exact BRI: **101**
- nonzero 1–10 mA: **1,753**

All 1,854 were reviewed.

### Exact-BRI findings

Observed situations include:

- ligand-model reuse
- legacy ligand-model reuse
- problematic redeposition
- probable redeposition
- clear redundant deposition
- alternative processing of the same model
- metadata/model conflict
- replicate same model
- partial protomer reuse
- entire model-1 structural reuse
- polymer-model reuse with non-polymer differences
- partial polymer-model reuse
- geometry-origin indeterminate
- trivial m=1 matches

Exact BRI equality was not automatically interpreted as erroneous experimental
duplication.

Stage-13 labels were not used as a global automatic deletion rule.

### Nonzero findings

Final nonzero categories:

- **1,423** trivial single-residue near matches
- **319** nonzero near-BRI geometry cases
- **11** canonical-alignment-resolved near-BRI cases

For the 11 resolved cases, canonical accepted-residue alignment gave
approximately 7e-6 A RMSD, with maximum atom deviations around
1.3e-5–1.5e-5 A.

Stage 13 is an investigation stage, not the global deduplication graph.

---

## Stage 14 — COMP702 representative-selection policy

This stage is explicitly a **COMP702 downstream choice**.

The Acta paper does not prescribe a representative-selection algorithm.

Global graph:

- m>=2 nodes: **577,760**
- frozen m>=2 near-duplicate edges: **1,068,256**
- edge-participating chains: **99,854**
- no-edge m>=2 chains: **477,906**
- edge components: **20,789**
- clique components: **20,494**
- non-clique components: **295**

m=1 is excluded from deduplication.

All **764 m=1 chains are retained**.

Connected-component membership is **not** treated as duplicate equivalence.

### Representative policy v1

Ranking:

1. prefer `terminal_trimmed = false`
2. prefer lower dirty-residue count
3. compare experimental resolution only when experimental-method tuples are
   directly comparable
4. canonical chain key as deterministic tie-break

No arbitrary global method preference such as “X-ray always beats EM” is used.

Clique component:

- one highest-ranked representative

Non-clique component:

- deterministic quality-ordered greedy direct-edge cover
- multiple representatives may exist

Critical invariant:

**Every removed chain must have its own direct `d_bri_mA <= 10` edge to the
specific representative to which it is assigned.**

Transitive connectedness alone is never sufficient for removal.

Results:

- edge-component representatives: **21,100**
- removed chains: **78,754**
- m>=2 no-edge retained: **477,906**
- m=1 retained: **764**
- final retained population: **499,770**

Independent audit confirmed all **78,754** required removal edges exactly.

---

## Final release

Release:

`outputs/releases/PDBClean-20260101-protocol3.2-comp702-v1-dedup-v1`

Population:

- canonical input: **578,524**
- removed: **78,754**
- retained: **499,770**

Retained dataset SHA256:

`8ae52ad96586c2552f74083b480350973c86bdcca41ae1f30f7353472d769c8b`

Removed-chain audit SHA256:

`4cb3bea6c6a61f27de60818d097cf72c0c047f603d13f76c8286bbae647d3360`

Release-manifest SHA256:

`1e6d6b249b6530fb501351fe6bd8d78647d3dad549db67e3fde486c2e3f8b918`

_SUCCESS SHA256:

`945c6c34358b127ea07365384f6f50429af315a26d9878233e4638cacf34c400`

---

## Interpretation boundary

The final product should be described as:

**a geometrically deduplicated PDB chain dataset under COMP702 representative
policy v1**

It should not be described as proof that every removed chain corresponds to an
incorrect or invalid experiment.

---

## Superseded artifacts

These remain provenance only:

- exhaustive Stage-8 search: validation oracle
- Stage-11 v1 EM-inclusive interpretation: superseded by v2
- `before_*.csv`: temporary manual-review checkpoints
- Stage-13 detailed subset: not the global Stage-14 graph

Generated large pipeline outputs are not duplicated into Git history. Their
paths, hashes, counts, and producer provenance are recorded instead.
