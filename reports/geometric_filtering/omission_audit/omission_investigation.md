# BRI Omission Investigation

## Scope

The model-1 PDB707K source contains 461,581 chains. BRI extraction produced 460,920 usable invariant vectors, leaving 661 omissions.

## Omission classification

| Omission class | Chains |
|---|---:|
| Missing local mmCIF | 413 |
| BRI residue-count mismatch | 218 |
| BRI invariant not generated | 27 |
| Ambiguous duplicate chain key | 3 |
| **Total** | **661** |

The counts were reproduced by comparing the full model-1 source records against the successful BRI extraction metadata and recorded extraction failures.

## 1. Missing local mmCIF files

### Finding

The 413 omitted chains correspond to 66 PDB entries that were absent from the local mmCIF archive used during BRI extraction:

`/users/sgsrehm1/scratch/bri_geometric/mmcif_raw`

Representative examples were:

- `1E7J`: 1 affected chain
- `8D35`: 2 affected chains
- `8T8O`: 204 affected chains

All 66 missing PDB entries returned HTTP status `200` from the official RCSB mmCIF download endpoint.

### Interpretation

These structures are not unavailable or withdrawn. The omission was caused by an incomplete local mmCIF mirror.

### Suggested solution

1. Download the 66 missing mmCIF files from RCSB.
2. Place them in the expected divided-directory layout.
3. Rerun BRI extraction for the affected 413 chains.
4. Add successful vectors to the geometric search.
5. Record any remaining failures under their new extraction failure reason.

Evidence:

- `missing_mmcif_rcsb_status.csv`
- `all_model1_bri_omissions.csv`

## 2. BRI residue-count mismatch

### Extraction behaviour

For each source record, the extraction script selects coordinates satisfying:

- matching model ID;
- matching label chain ID;
- residue ID greater than or equal to `start_residue`;
- residue ID below `start_residue + chain_length`.

BRI then uses `ATOM` records for the backbone atoms `N`, `CA`, and `C`.

The extraction is rejected when the number of generated invariant rows differs from the source `chain_length`.

### Representative sample findings

#### 1AAM, chain A

- Source interval: residues 1–396
- Expected residues: 396
- Coordinate residue IDs present: 395
- BRI invariant rows: 395
- Missing backbone atoms among present residues: 0
- Deficit: 1 residue

#### 7CLR, chain FA

- Source interval: residues 22–232
- Expected residues: 211
- Coordinate residue IDs present: 200
- BRI invariant rows: 200
- Missing backbone atoms among present residues: 0
- Deficit: 11 residues

#### 8RX0, chain H

- Source interval: residues 1–176
- Expected residues: 176
- Coordinate residue IDs present: 71
- Observed coordinate range: residues 15–92
- BRI invariant rows: 71
- Missing backbone atoms among present residues: 0
- Deficit: 105 residues

### Verified interpretation

For these representative samples, the mismatch was not caused by missing `N`, `CA`, or `C` atoms within residues that were present.

Instead, the source dataset expected a continuous residue-number interval, but some residue IDs in that interval were entirely absent from the deposited coordinate records.

The severe `8RX0-H` mismatch may additionally indicate that the source chain-range record does not correspond cleanly to the coordinate chain currently present in the mmCIF file.

### Investigation still required

The following must be checked across the 218 failures:

1. Exact missing residue IDs.
2. Whether missing positions are terminal or internal.
3. Whether the source interval represents polymer sequence length rather than observed coordinate length.
4. Whether chain identifiers or residue ranges changed between dataset preparation and the current mmCIF version.
5. Whether some mismatches can be recovered using coordinate-derived chain boundaries.

### Preliminary solution options

- For small terminal losses, derive the usable coordinate range and document the adjustment.
- For internal coordinate gaps, avoid treating the chain as one uninterrupted geometric sequence without further validation.
- For large mismatches, verify chain mapping and dataset metadata before including the structure.
- Do not silently replace the source length with the observed length, because this changes the dataset definition.

## 3. BRI invariant not generated

There are 27 chains in this class.

Direct reruns and coordinate-level diagnostics identified four root-cause
subclasses:

- 11 chains with duplicate backbone atom records;
- 8 chains with source residue-range mismatches;
- 5 non-protein nucleic-acid chains;
- 3 chains with author-chain/label-chain namespace mismatches.

The detailed evidence and affected chains are documented in the later sections
of this report.

## 4. Ambiguous duplicate chain keys

Three omitted records belong to `8RX0`:

| Chain | Start residue | Chain length |
|---|---:|---:|
| B | 1 | 104 |
| F | 1 | 77 |
| I | 1 | 111 |

The source CSV contains repeated PDB/model/chain identifiers with different residue ranges.

The production manifest retained the alternative records:

- `8RX0-B`: start 8, length 148
- `8RX0-F`: start 1, length 140
- `8RX0-I`: start 1, length 75

### Interpretation

A key consisting only of PDB ID, model ID, and chain ID is insufficient to distinguish these source records.

### Suggested solution

Use a composite record identifier that includes at least:

- PDB ID
- model ID
- label chain ID
- start residue
- chain length

The retained and omitted records should also be checked against label-chain, author-chain, entity and polymer metadata in the mmCIF file before the final dataset is produced.

## Current conclusion

The 661 omissions are not one homogeneous failure group.

- 413 are recoverable local archive omissions.
- 218 are discrepancies between expected source residue intervals and observed coordinates.
- 27 are BRI computation failures requiring atom-level investigation.
- 3 arise from ambiguous duplicated source identifiers.

No omitted chains should be removed permanently until the recoverable cases and metadata inconsistencies have been investigated.

### Residue-gap localisation

The exact absent residue positions were examined for the three representative
residue-count mismatch cases.

| Chain | Observed range | Internal gaps | Terminal omissions |
|---|---:|---|---|
| 1AAM-A | 1–396 | Residue 246 | None |
| 7CLR-FA | 22–232 | Residues 146–156 | None |
| 8RX0-H | 15–92 | Residues 48–54 | Residues 1–14 and 93–176 |

#### Interpretation

`1AAM-A` and `7CLR-FA` cover the expected terminal boundaries but contain
unresolved internal coordinate gaps. Their source chain lengths therefore count
polymer positions for which no atomic coordinates were deposited.

`8RX0-H` is different: only residues 15–92 are observed, with an additional
internal gap at residues 48–54. Of the 105 absent positions, 98 are terminal and
7 are internal. This is evidence of substantial experimental truncation or a
mismatch between the PDB707K source interval and the coordinate chain.

These samples show that residue-count mismatches should be separated into at
least:

1. isolated internal coordinate gaps;
2. contiguous internal unresolved regions; and
3. major terminal truncation or source-to-coordinate range mismatch.

Chains containing internal gaps should not be silently converted into continuous
BRI vectors because residues on opposite sides of an unresolved region are not
necessarily geometrically contiguous. Major truncations require source metadata
and chain-mapping validation before recovery.

### Invariant-generation failure messages

All 27 chains in the `bri_invariant_not_generated` class recorded the same
exception:

`AttributeError("'NoneType' object has no attribute 'copy'")`

This is not the original BRI failure. The `MiniChain.invariant` property catches
the underlying exception raised by `get_chain_invariant()`, emits a warning, and
returns `None`. The extraction script then calls `.copy()` on this `None` value,
producing the recorded `AttributeError`.

Consequently, the production failure table does not preserve the actual
atom-level or numerical cause of these 27 failures.

The affected records range from chains of length 1 to chains of length 500.
Therefore, insufficient chain length can explain some records, but it cannot
explain the entire class.

Representative cases were rerun by bypassing the
exception-swallowing `invariant` property and calling
`get_chain_invariant()` directly. This recovered the original exceptions and
enabled the subsequent root-cause classification.

### Direct investigation of invariant-generation failures

The representative chains were rerun using
`MiniChain.get_chain_invariant()` directly, bypassing the property that hides
the original exception.

#### 5CGH-CA

- Source start residue: 4
- Source chain length: 1
- Coordinate residues selected: 0
- Backbone atoms selected: 0
- Root exception: `IndexError: positional indexers are out-of-bounds`

The requested chain segment is empty. This is not primarily a chain-length
failure: the source model, label-chain identifier or residue range does not
select any coordinates in the current mmCIF file.

#### 7T3M-G

- Source start residue: 1
- Source chain length: 110
- Coordinate residues selected: 0
- Backbone atoms selected: 0
- Root exception: `IndexError: positional indexers are out-of-bounds`

Despite the expected length of 110 residues, the selected segment is empty.
This indicates a source-to-mmCIF chain mapping or residue-range mismatch.

#### 5U8X-C

- Source start residue: 27
- Source chain length: 500
- Coordinate residues selected: 500
- Backbone atom rows: 1,503
- Problematic residue: 282
- Atom counts at residue 282: two `N`, two `CA`, and two `C` records
- Root exception:
  `ValueError: Index contains duplicate entries, cannot reshape`

This failure is caused by duplicate backbone atom records rather than absent
coordinates. Possible causes include alternate conformations or duplicated atom
records that were not resolved before BRI reshaping.

### Refined classification

The broad `bri_invariant_not_generated` category contains at least two verified
subclasses:

1. **Empty source-selected coordinate segment**
   - Example: `5CGH-CA`
   - Example: `7T3M-G`
   - Likely solution: verify label-chain, author-chain, entity, model and residue
     mappings against the current mmCIF file.

2. **Duplicate backbone atom records**
   - Example: `5U8X-C`, residue 282
   - Likely solution: inspect alternate-location identifiers and occupancies,
     select one conformer deterministically, and rerun BRI extraction.

The original extraction output obscured these causes because the underlying
exception was replaced by a `.copy()` error after BRI returned `None`.
Future extraction runs should explicitly test for `None` and preserve the
original exception.

### Chain-mapping investigation

#### 5CGH-CA

`CA` exists as a label-chain identifier and maps to author chain `c`. However,
the chain contains one observed coordinate residue numbered 3, whereas the
source record requests start residue 4 with length 1.

The empty coordinate segment is therefore caused by a source-to-mmCIF residue
numbering mismatch, rather than a missing chain.

#### 7T3M-G

No label chain named `G` exists in model 1. Instead, `G` is an author-chain
identifier mapped to label chain `D`, entity 1.

The extractor compared the source value `G` against the mmCIF label-chain field,
so it selected zero atoms. This verifies an author-chain versus label-chain
identifier namespace mismatch.

#### Implication

Source `chain_id` values cannot always be assumed to use the mmCIF
`label_asym_id` namespace. Extraction should first determine whether each source
identifier is a label-chain ID or an author-chain ID and, when necessary, map
the author identifier to its corresponding label chain.

Residue ranges must also be validated independently, because a valid chain
mapping can still produce an empty segment when source and current mmCIF
residue numbering disagree.

### Full root-cause classification of invariant-generation failures

All 27 chains were rerun using `MiniChain.get_chain_invariant()` directly and
examined for chain mapping, residue-range selection and backbone atom counts.

| Root-cause class | Chains |
|---|---:|
| Duplicate backbone atom records | 11 |
| Source residue-range mismatch | 8 |
| Other BRI computation error | 5 |
| Author-chain/label-chain namespace mismatch | 3 |
| **Total** | **27** |

#### Duplicate backbone atom records: 11 chains

Affected chains:

- `2AAA-A`
- `2Z1Y-B`
- `5U8X-A`
- `5U8X-C`
- `5U8X-D`
- `6TVW-B`
- `6XND-A`
- `6XND-B`
- `6XND-C`
- `6XND-D`
- `8RXA-B`

These segments contain more than one copy of at least one backbone atom (`N`,
`CA` or `C`) for the same residue ID. BRI then fails while reshaping the
residue-by-atom table, typically with:

`ValueError: Index contains duplicate entries, cannot reshape`

Some records also have fewer coordinate residues than the source length, such
as `6TVW-B` and `8RXA-B`. Therefore, duplicate atoms are the immediate BRI
failure, but residue-count discrepancies may also need correction.

The likely recovery procedure is to inspect alternate-location identifiers and
occupancies, select one conformer deterministically, and rerun BRI extraction.

#### Source residue-range mismatch: 8 chains

Affected chains:

- `5CGH-CA`
- `5CGH-DA`
- `6NSR-E`
- `6QNX-B`
- `6QNX-C`
- `7XYG-G`
- `7Y6K-C`
- `8G94-G`

For each record, the source chain identifier exists or is related to a chain in
the mmCIF file, but the requested source residue interval selects zero
coordinates.

This indicates disagreement between the PDB707K start/length fields and the
current mmCIF residue numbering or chain mapping.

The likely recovery procedure is to validate the source interval against
label-sequence IDs, author residue IDs, entity IDs and observed coordinate
ranges before extraction.

#### Author-chain versus label-chain namespace mismatch: 3 chains

Affected chains:

- `7T3M-G`, mapped to label chain `D`
- `7T3M-H`, mapped to label chain `E`
- `7T3M-I`, mapped to label chain `F`

The source identifiers are author-chain IDs, but the extraction script compares
them directly against the mmCIF label-chain field. This results in empty
coordinate selections.

These chains may be recoverable by applying the author-to-label chain mapping
before selecting their coordinates. Their residue ranges must still be checked
after remapping.

#### Other BRI computation errors: 5 chains

Affected chains:

- `7XYF-E`
- `7XYF-F`
- `7XYG-E`
- `7XYG-F`
- `8GWG-I`

These chains contain selected coordinates and no simple missing or duplicated
backbone atoms were detected, but direct BRI calculation fails with an
`IndexError`.

`8GWG-I` is additionally notable because only 23 coordinate residues were
selected from an expected source length of 113.

These five cases require inspection of their exact tracebacks, atom ordering,
residue continuity and chain-identifier mapping before a solution can be
recommended.

### Implication for the extraction pipeline

The original class `bri_invariant_not_generated` was too broad and its stored
`.copy()` exception hid the actual failure.

Future extraction should:

1. call BRI in a way that preserves the root exception;
2. explicitly detect empty coordinate selections;
3. resolve author-chain and label-chain namespaces;
4. detect duplicate backbone atoms before invariant calculation;
5. record observed residue count and range for every failure.

### Non-protein chains in the BRI input

The five previously unresolved BRI failures were examined at the residue and
atom level.

| Chain | Polymer content | Observed residues |
|---|---|---:|
| 7XYF-E | DNA | 106 |
| 7XYF-F | DNA | 105 |
| 7XYG-E | DNA | 106 |
| 7XYG-F | DNA | 105 |
| 8GWG-I | RNA | 23 |

The four DNA chains contain nucleotide residue labels such as `DA`, `DC`, `DG`
and `DT`. The RNA chain contains `A`, `C`, `G` and `U`.

All five consist of valid `ATOM` records, but none contain the protein backbone
atom set `N`, `CA` and `C` required by the protein BRI implementation.
Consequently, `get_feature()` returns an empty table and the invariant code
fails when it assumes that protein coordinate columns are present.

These are not unexplained numerical BRI errors. They are non-protein polymer
chains that were passed into a protein-specific geometric representation.

`8GWG-I` also contains only 23 observed coordinate residues despite a source
chain length of 113. It may therefore have an additional source-range or
chain-mapping discrepancy, but its immediate BRI failure is caused by its RNA
polymer type.

#### Suggested solution

The extraction manifest should be filtered by polymer type before BRI
calculation. Only polypeptide chains with an appropriate amino-acid residue
alphabet and usable `N`, `CA` and `C` backbone atoms should enter the protein
BRI pipeline.

DNA, RNA and other non-polypeptide polymers should be excluded explicitly with
their polymer type recorded as the reason, rather than being reported as BRI
computation failures.

### Final classification of the 27 invariant-generation failures

| Root-cause class | Chains |
|---|---:|
| Duplicate backbone atom records | 11 |
| Source residue-range mismatch | 8 |
| Non-protein nucleic-acid chain | 5 |
| Author-chain/label-chain namespace mismatch | 3 |
| **Total** | **27** |

All 27 records now have a concrete root-cause class.

### Detailed investigation of ambiguous 8RX0 records

The source dataset contains repeated model-1 records for label chains `B`, `F`
and `I`, using the same PDB/model/chain identifier but different residue
ranges.

| Chain | Source interval | Coordinate coverage | Production outcome |
|---|---:|---:|---|
| B | start 1, length 104 | 97/104 | Omitted as ambiguous |
| B | start 8, length 148 | 148/148 | Successful BRI vector |
| F | start 1, length 77 | 77/77 | Omitted as ambiguous |
| F | start 1, length 140 | 140/140 | Successful BRI vector |
| I | start 1, length 75 | 75/75 | Successful BRI vector |
| I | start 1, length 111 | 75/111 | Omitted as ambiguous |

#### Chain B

The deposited label chain contains residues 8–155, totalling 148 observed
residues. The retained record, start 8 and length 148, matches the complete
observed chain exactly.

The omitted start-1, length-104 record expects residues 1–104, but residues 1–7
are absent. It also represents only part of the deposited chain.

#### Chain F

The deposited chain contains residues 1–140.

Both source intervals select valid coordinate segments:

- start 1, length 77 selects a valid 77-residue prefix;
- start 1, length 140 selects the complete observed chain.

The shorter record is therefore not structurally invalid, but it is a redundant
subchain representation under the same PDB/model/chain identifier. The
production manifest retained the complete 140-residue chain.

#### Chain I

The deposited chain contains residues 1–75. The retained length-75 record
matches it exactly.

The omitted length-111 record overstates the observed chain by 36 residues and
has coordinate coverage of only 75/111.

#### Interpretation

The `ambiguous_duplicate_chain_key` category combines two related but distinct
problems:

1. source records whose residue ranges disagree with the observed chain
   (`8RX0-B` length 104 and `8RX0-I` length 111);
2. a valid but redundant subsegment of a longer record for the same chain
   (`8RX0-F` length 77).

A source key consisting only of PDB ID, model ID and chain ID is insufficient
when multiple residue ranges are present.

#### Suggested solution

For final dataset construction:

1. use a composite source-record key including start residue and chain length;
2. compare every competing record with the observed coordinate range;
3. prefer the complete fully covered chain when one record is a strict prefix
   of another;
4. retain partial segments only when the scientific objective explicitly
   requires domain or fragment-level records;
5. record the disambiguation decision rather than silently dropping duplicates.

## Recommended recovery plan

| Omission class | Chains | Recommended action |
|---|---:|---|
| Missing local mmCIF | 413 | Download the 66 available RCSB files and rerun extraction |
| BRI residue-count mismatch | 218 | Classify internal gaps, terminal truncations and metadata mismatches before deciding whether to split, correct or exclude each chain |
| Duplicate backbone atom records | 11 | Resolve alternate locations or duplicated atoms deterministically, then rerun BRI |
| Source residue-range mismatch | 8 | Reconcile source intervals with current mmCIF label and author residue numbering |
| Non-protein nucleic-acid chain | 5 | Exclude explicitly from the protein BRI pipeline and record polymer type |
| Author/label chain namespace mismatch | 3 | Map author-chain identifiers to label-chain identifiers and validate residue ranges |
| Ambiguous duplicate source record | 3 | Prefer the complete, fully covered chain and record each disambiguation decision |
| **Total** | **661** | |

### Priority

1. Recover the 413 chains missing only because the local archive was incomplete.
2. Repair deterministic metadata and atom-selection issues.
3. Re-extract all recoverable chains.
4. Recalculate the geometric near-duplicate search using the expanded successful set.
5. Preserve unresolved exclusions in a final omission table with evidence and reason codes.

### Reproducibility improvements

Future extraction runs should record:

- source and resolved chain identifier namespaces;
- expected and observed residue ranges;
- observed residue count;
- missing internal and terminal residue positions;
- duplicate backbone atoms and alternate-location selection;
- polymer type;
- the original BRI exception rather than the secondary `.copy()` error;
- the recovery or exclusion decision for every failed chain.
