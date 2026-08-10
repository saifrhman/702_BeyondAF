# PDBClean Pipeline Specification

## 1. Purpose

This project builds a reproducible, auditable, restartable data pipeline that
converts a fixed wwPDB snapshot into:

1. A quality-controlled protein-chain dataset.
2. BRI representations for accepted chains.
3. Geometric near-duplicate pairs identified using cKDTree and strict
   L-infinity distance.
4. A final geometrically deduplicated protein-chain dataset.

The historical PDB707K CSV is treated as the output of an older PDB snapshot,
not as a permanent source dataset. The new pipeline regenerates an equivalent
dataset from an explicitly versioned snapshot.

## 2. Dataset Versioning

Every dataset release must be identified by:

- PDB snapshot date.
- Cleaning protocol version.
- Pipeline Git commit.
- Configuration checksum.

Example release name:

PDBClean-20260101-protocol3.2-comp702-v1

Published releases are immutable. Changing the snapshot, cleaning rules,
configuration, or representative-selection policy creates a new release rather
than modifying an existing release.

## 3. Data Layers

### Bronze Layer

The Bronze layer stores immutable source metadata for the selected PDB
snapshot.

It contains:

- Snapshot date.
- PDB ID.
- S3 object key.
- Compressed file size.
- ETag.
- Manifest generation timestamp.

No structural filtering occurs in this layer.

### Silver Layer

The Silver layer is the deterministic in-memory parsed representation used
between immutable source mmCIF objects and Gold quality decisions.

It is not published as a separate persisted dataset in this pipeline version.
Instead, Silver observations are reconstructed from the exact immutable source
object identified by the Bronze manifest S3 key, byte size and ETag.

The parsed representation preserves, for every deposited chain and model:

- PDB ID.
- Model ID.
- Entity ID.
- Label chain ID.
- Author chain ID.
- Polymer type.
- Residue identifiers and names.
- Full atom observations required by the parser.
- Raw occupancy information.
- Alternate-location information.
- Source-object lineage.

Persisting a second copy of the structural snapshot is unnecessary for
reproducibility because the immutable manifest identifies the exact source
bytes and the versioned parser defines the deterministic Bronze-to-Silver
transformation.

### Gold Layer

The Gold layer contains the validated research datasets.

It includes:

- Quality-accepted chains.
- Quality-rejected chains with explicit reasons.
- BRI vectors.
- BRI extraction failures.
- Strict L-infinity near-duplicate pairs.
- Graph components.
- Representative-selection decisions.
- Final retained and removed chains.

No chain may disappear silently. Every candidate chain must be recorded as
accepted, rejected, or failed during processing.

## 4. Snapshot Input

The initial source dataset is the complete wwPDB snapshot dated 2026-01-01.

Source prefix:

20260101/pub/pdb/data/structures/divided/mmCIF/

Verified snapshot manifest:

- 246,905 compressed mmCIF files.
- 85,079,649,893 total bytes.
- 79.24 GiB compressed.

A dated S3 prefix is not automatically treated as a complete structural
snapshot. The pipeline must verify that the divided mmCIF directory exists and
contains a plausible complete set of coordinate files before processing begins.

The snapshot manifest is immutable and records the exact objects used by the
pipeline, including their S3 keys, sizes, and ETags.

## 5. Canonical Identifiers

The canonical structural-chain key is:

(snapshot, pdb_id, model_id, label_chain_id)

The label chain ID is canonical because atom-site records use the label
namespace consistently for structural processing.

The pipeline must also preserve:

- Author chain ID.
- Entity ID.
- Source mmCIF key.
- Source ETag.

The initial release processes model 1 only. This is a versioned configuration
choice and must not be hard-coded into the parser.

### 5.1 Model Selection

Coordinate mmCIF parsing preserves all deposited models in Silver.

Model selection is a separate, configuration-controlled processing step and
must not be hard-coded into the parser.

For `protocol3.2-comp702-v1`, the configured policy is:

- `selection.models.policy: first_model`;
- `selection.models.model_id: 1`.

Therefore only parsed chains whose `model_id` equals the configured
`selection.models.model_id` proceed to Protocol 3.2 candidate selection and
quality cleaning.

Chains from other models are outside the processing scope of this release.
Their parsed observations exist only in the in-memory Silver representation
during source processing and remain reproducible from the immutable source
mmCIF. They are not emitted as Protocol 3.2 non-candidates, rejections, or
accepted chains.

The selected model identifier must remain part of the canonical chain key and
all downstream Gold and BRI provenance.

### 5.2 Candidate Accounting

Candidate accounting is evaluated after configured model selection and before
Q001.

A selected-model chain is a Protocol 3.2 candidate when its BRI-style
backbone projection is non-empty; that is, selecting `ATOM` rows whose atom
names are `N`, `CA`, or `C` produces at least one row.

Therefore:

- `candidate_chain_count` is the number of selected-model chains with a
  non-empty Protocol 3.2 backbone projection;
- `candidate_entry_count` is the number of unique `pdb_id` values having at
  least one candidate chain;
- selected-model chains with an empty projection are counted as
  non-candidates instead.

Candidate status is established before Q001. An entry may therefore count as
a candidate entry and subsequently have its candidate chains rejected by the
entry-level Q001 protein check.

For an entry with multiple candidate chains, the entry contributes one to
`candidate_entry_count` and one per candidate chain to
`candidate_chain_count`.

## 6. Protocol 3.2 Quality Cleaning

The scientific cleaning contract is Protocol 3.2 of Anosova et al. The pinned
BRI v1.2.2 implementation is the executable specification for operational
details that are not fully specified by the paper.

The rules in this protocol version are therefore not independent optional
filters. They form one stateful cleaning procedure. Residues can be marked
dirty and removed while a remaining continuous chain segment is retained.

For auditability, the pipeline preserves:

- the immutable source-object identity needed to reproduce the original
  deposited observations through the deterministic Silver parser; and
- all residue removals and terminal chain outcomes in persisted Gold outputs.

### 6.1 Backbone Projection

Before the chainwise cleaning functions are applied, BRI obtains:

`get_feature("features", HETATM=False)`

This restricts the working coordinates to `ATOM` records whose atom names are:

- N;
- CA;
- C.

Side-chain atoms and HETATM records are therefore outside the operational
backbone-cleaning input.

The in-memory Silver observations remain unchanged; the Protocol 3.2
projection is derived without mutating them.

A Silver chain becomes a Protocol 3.2 chain candidate only when this
BRI-style projection is non-empty.

Therefore:

- HETATM-only chains remain represented in the in-memory Silver parse but
  are not Protocol 3.2 candidates;
- ATOM chains containing no N, CA, or C rows remain represented in the
  in-memory Silver parse but are not Protocol 3.2 candidates;
- these non-candidate chains are not Q001-Q006 scientific rejections and are
  not processing errors.

This reproduces BRI v1.2.2 `on_entry()`, which skips chains when
`get_feature("features", HETATM=False)` returns no target backbone rows.

### Q001: Entry Contains a Polypeptide

Protocol 3.2 excludes non-protein entries.

In pinned BRI v1.2.2, `field_check()` implements this as an entry-level test.
An entry is considered peptide-containing when any `_entity_poly.type` begins
with:

`polypeptide`

This is not a per-chain whitelist. In a mixed entry containing both a
polypeptide and another polymer type, Q001 passes at entry level and the
available ATOM chains continue to the chainwise cleaning stage.

For Q001 record:

- PDB ID;
- entry-level peptide classification;
- observed entity polymer types;
- terminal reason when the entry contains no polypeptide entity.

### Q002: Disorder

The paper identifies chains containing atoms with occupancy below 1 as
disordered.

The pinned BRI v1.2.2 operational implementation applies
`disorder_check()` to the N/CA/C backbone projection.

A backbone residue is marked dirty when either:

1. more than one row has the same
   `(residue_id, chain_id, model_id, atom)` identity; or
2. the raw occupancy token does not start with `"1"`, except that the literal
   `"."` token is explicitly tolerated.

Consequences of the executable semantics include:

- `"1.00"` is accepted;
- `"0.50"` is dirty;
- `"."` is accepted;
- `"?"` is dirty;
- alternate-location labels are not an independent rejection criterion;
- partial-occupancy side-chain atoms are outside this check because only
  N/CA/C rows are supplied to it.

All residues marked dirty by Q002 are removed from the working chain.

For every Q002 dirty residue record:

- residue identifier;
- deposited residue name;
- affected backbone atom;
- raw occupancy token;
- duplicate-backbone status;
- alternate-location identifier when present;
- Q002 reason.

### Q003: Residue Continuity

BRI maps `_atom_site.label_seq_id` to its internal `residue_id`.

`residue_continuity_check()` constructs the integer range from the minimum
observed residue ID through the maximum and tests for absent positions.

Missing residue ID `0` is explicitly ignored by the BRI implementation.

A continuity failure is a chain-level failure: the current working chain does
not produce a clean retained chain.

Q003 is not applied only once. It is rerun after residue-removal stages as
described in Section 6.2.

For every Q003 chain-break record:

- current retained start residue;
- current retained end residue;
- missing residue identifiers;
- cleaning stage at which the break was detected;
- Q003 reason.

### Q004: Backbone Completeness

Protocol 3.2 requires the main backbone atoms N, CA and C.

In pinned BRI v1.2.2, `residue_completeness_check()` operates after the chain
has already been projected to N/CA/C rows. It counts rows per `residue_id` and
marks the residue dirty whenever:

`backbone row count != 3`

Duplicate backbone atom records are normally identified earlier by Q002.
Operationally, however, Q004 reproduces the BRI row-count test exactly.

Q004 dirty residues are removed from the working chain together with Q006
dirty residues before continuity is checked again.

For every Q004 dirty residue record:

- residue identifier;
- deposited residue name;
- observed N/CA/C row count;
- available backbone atom names;
- Q004 reason.

### Q005: Backbone Clash

Protocol 3.2 removes backbone clashes below 0.01 angstrom.

Pinned BRI v1.2.2 uses:

`clash_check(chain, 0.01)`

with the local atom set:

`N(i), CA(i), C(i), N(i+1)`

and six pair types:

- N(i) - CA(i);
- N(i) - C(i);
- N(i) - N(i+1);
- CA(i) - C(i);
- CA(i) - N(i+1);
- C(i) - N(i+1).

These are six comparison types, not six total measurements per chain.

For a continuous chain containing `L` residues, the implementation evaluates:

`3L + 3(L - 1) = 6L - 3`

distances.

For example, a two-residue chain produces nine evaluated distances.

The threshold is strict:

`distance < 0.01 angstrom`

is dirty, while a distance exactly equal to 0.01 angstrom is accepted.

Residues returned by the clash check are marked dirty and removed, after which
continuity is checked once more.

For every Q005 dirty record:

- residue identifier returned by the clash implementation;
- compared atom labels;
- measured distance;
- threshold;
- Q005 reason.

### Q006: Standard Amino-Acid Mapping

Protocol 3.2 removes non-standard amino acids.

Pinned BRI v1.2.2 does not implement this as a direct three-letter whitelist.
Before `standard_residue_check()` is called, deposited residue labels are
mapped using the package's bundled CCD-derived:

`amino_acid_short`

mapping.

The mapped one-letter code must belong to the canonical 20-code set:

`A R N D C Q E G H I L K M F P S T W Y V`

This means the pinned executable implementation includes mappings such as:

- LLP -> K;
- MSE -> M;

which pass the standard-residue check, while mappings such as:

- SEC -> U;
- PYL -> O;
- UNK -> X;

do not pass.

The original deposited residue name must always be preserved in Gold
diagnostics even when a mapped one-letter code is used for the cleaning
decision.

Q006 dirty residues are removed together with Q004 dirty residues.

There is a known reference discrepancy concerning LLP: direct correspondence
described LLP as an example that would be removed, whereas the pinned BRI
v1.2.2 CCD mapping maps LLP to K and therefore accepts it. The pipeline follows
the pinned executable behaviour for reproducibility and records this
discrepancy explicitly rather than silently changing the implementation.

### 6.2 Stateful Execution Order

The Protocol 3.2 cleaning implementation must reproduce the stateful order of
BRI v1.2.2.

For a candidate entry and chain:

1. Apply Q001 at entry level.
2. Project the chain to ATOM N/CA/C rows.
3. Map deposited residue labels through the BRI CCD mapping.
4. Detect Q002 disorder residues.
5. Record and remove Q002 dirty residues.
6. Apply Q003 continuity to the remaining chain.
7. Detect Q006 non-standard residues.
8. Detect Q004 incomplete residues.
9. Record and remove the union of Q006 and Q004 dirty residues.
10. Apply Q003 continuity again.
11. Detect Q005 clashes using the strict 0.01 angstrom threshold.
12. Record and remove Q005 dirty residues.
13. Apply Q003 continuity a final time.
14. Mark the remaining chain clean.

The order is scientifically important. Q002, Q004, Q005 and Q006 do not
automatically reject the whole original chain.

For example, if a defective terminal residue is removed and residues remaining
from the new minimum to maximum identifier are consecutive, the trimmed chain
can still be retained.

If removal of a defective internal residue creates a missing position inside
the remaining residue range, the subsequent Q003 check detects a chain break
and no clean chain is emitted from that working chain.

If all working residues are removed, no clean chain is emitted.

### 6.3 Paper Criteria Versus Executable Details

The paper is the scientific source of truth for the Protocol 3.2 criteria and
reported dataset outcome.

BRI v1.2.2 is used as the executable source of truth for details including:

- the N/CA/C feature projection;
- raw occupancy-token handling;
- duplicate-backbone detection;
- `label_seq_id` continuity semantics;
- special handling of residue ID 0;
- row-count implementation of backbone completeness;
- CCD-based residue mapping;
- clash pair construction;
- strict clash threshold;
- stateful residue removal and continuity rechecking.

Where direct correspondence, the paper and the pinned executable
implementation appear to disagree, the discrepancy must be documented rather
than resolved by an undocumented local rule.

## 7. Quality-Cleaning Outputs

Quality-task shards are intermediate stage outputs, distinct from the final
published release.

The quality-stage output root is derived from configuration and runtime
provenance as:

`<storage.output_root>/<snapshot>/<release.protocol_version>/quality/`

Neither the snapshot identifier nor protocol version is hard-coded into the
quality output root is derived from configuration and runtime
provenance as:

`<storage.output_root>/<snapshot>/<release runner.

Within that stage root, each SLURM task writes independent Gold output shards.

Clean retained chains:

`accepted/task_<task_id>.parquet`

Chain-level rejections:

`rejected/task_<task_id>.parquet`

Chains outside the Protocol 3.2 candidate boundary:

`non_candidates/task_<task_id>.parquet`

Dirty-residue lineage:

`dirty_residues/task_<task_id>.parquet`

Processing errors:

`errors/task_<task_id>.parquet`

Task summaries:

`summaries/task_<task_id>.json`

After successful completion of all tasks, these shards are merged into:

- `accepted_chains.parquet`;
- `rejected_chains.parquet`;
- `non_candidate_chains.parquet`;
- `dirty_residues.parquet`;
- `processing_errors.parquet`;
- `quality_summary.json`.

### 7.1 Accepted Chains

An accepted-chain record describes the clean remainder after the full stateful
Protocol 3.2 procedure.

It must include:

- canonical source chain key;
- original Protocol 3.2 backbone-projection start and end
  `label_seq_id` values;
- retained start and end `label_seq_id` values;
- retained residue count;
- retained residue identifiers in surviving backbone row order;
- retained sequence using the pinned BRI residue mapping;
- whether terminal trimming occurred;
- number of dirty residues removed;
- dirty rule identifiers encountered before acceptance;
- source mmCIF key;
- source ETag;
- cleaning protocol version;
- pipeline Git commit.

An accepted chain may therefore also have records in `dirty_residues.parquet`.

### 7.2 Dirty Residues

Dirty-residue records preserve every residue-level cleaning decision made
before a final chain outcome is known.

Each record must include:

- canonical source chain key;
- residue identifier;
- deposited residue name;
- mapped residue code when applicable;
- dirty rule identifier;
- dirty type compatible with the BRI concepts
  (`disordered`, `incomplete`, `non-standard`, or `clash`);
- structured rule-specific details;
- cleaning stage;
- source lineage.

Dirty-residue records are not synonymous with rejected-chain records.

A chain can contain dirty terminal residues, have them removed, remain
continuous and still be accepted.

### 7.3 Rejected Chains

A rejected-chain record represents a candidate chain for which no clean
retained chain is emitted.

Examples include:

- entry rejected by Q001;
- Q003 chain break before or after residue removal;
- all working residues removed during cleaning.

Each record must include:

- canonical source chain key where available;
- terminal status;
- terminal reason;
- stage at which processing terminated;
- relevant missing residue identifiers for a chain break;
- number and types of dirty residues already recorded;
- source mmCIF key;
- source ETag.

### 7.4 Non-Candidate Chains

A non-candidate record represents a parsed Silver chain whose BRI-style
Protocol 3.2 backbone projection is empty.

Specifically, a Silver chain is outside the Protocol 3.2 candidate boundary
when selecting `ATOM` rows whose atom names are `N`, `CA`, or `C` produces no
rows.

Examples include:

- HETATM-only chains;
- ATOM chains containing no N/CA/C backbone rows.

Non-candidate chains are preserved for accounting but are not scientific
Protocol 3.2 rejections, because the chainwise cleaning rules were never
applied to them.

Each non-candidate record must include:

- canonical source chain key;
- terminal status `non_candidate`;
- terminal reason;
- candidate-selection stage;
- source mmCIF key;
- source ETag;
- cleaning protocol version;
- pipeline Git commit.

### 7.5 Processing Errors

Processing errors are reserved for technical failures rather than scientific
rejections or non-candidate classifications.

Examples include:

- unreadable or malformed mmCIF;
- missing required structural fields;
- unexpected parser failure;
- task-level I/O failure.

Scientific filtering outcomes and candidate-boundary classifications must
never be recorded as processing errors.

### 7.6 Complete Accounting

No parsed Silver chain may disappear silently.

Every parsed Silver chain must have exactly one terminal chain-level outcome:

- non-candidate;
- accepted;
- rejected; or
- processing error.

Every Protocol 3.2 candidate chain must therefore terminate as accepted,
rejected, or processing error. A non-candidate chain must never be counted as
a scientific rejection.

Residue-level dirty records are additional lineage and do not replace the
terminal chain outcome.

The quality summary must reconcile:

- Silver parsed chains;
- non-candidate chains with empty Protocol 3.2 backbone projections;
- candidate entries;
- candidate chains;
- Q001 entry-level exclusions;
- accepted clean chains;
- accepted chains that were trimmed;
- rejected chains by terminal reason;
- dirty residues by rule/type;
- processing errors.

### 7.6.1 Processing-Error Stage Vocabulary

Quality-stage processing errors use stable stage identifiers:

- `source_download_verify`: the source object could not be downloaded and
  verified against its manifest metadata. This includes transport failures,
  byte-count mismatches, missing response ETags, and ETag mismatches. These
  are source-level errors and therefore have null `model_id` and
  `label_chain_id`.
- `mmcif_parse`: verified source bytes could not be parsed into Silver chain
  observations. These are source-level errors and therefore have null
  `model_id` and `label_chain_id`.
- `quality_cleaning`: an individual selected-model chain encountered an
  unexpected exception during Protocol 3.2 cleaning or Gold
  materialization. These are chain-level errors and carry both `model_id`
  and `label_chain_id`.

The downloader performs transfer and manifest verification as one operation,
so download and verification failures are deliberately not classified by
parsing exception-message text.

### 7.7 Task Summary Contract

Each quality-cleaning task must atomically write:

`summaries/task_<task_id>.json`

The task summary is an observability and accounting record. It is not a
replacement for the Parquet Gold outputs.

The JSON object must include:

- summary schema name and version;
- task identifier;
- snapshot identifier;
- cleaning-protocol version;
- pipeline Git commit;
- start and completion timestamps in UTC;
- runtime in seconds;
- SLURM job identifier when available;
- SLURM array-task identifier when available;
- peak memory usage when available.

The summary must include the following source-processing counts:

- input source-object count;
- successfully processed source-object count;
- failed source-object count.

The source-object accounting invariant is:

`input_source_object_count = successful_source_object_count + failed_source_object_count`

The summary must include the following chain and residue counts:

- parsed Silver chain count across all deposited models;
- selected Silver chain count after configured model selection;
- candidate entry count;
- candidate chain count;
- non-candidate chain count;
- accepted chain count;
- accepted trimmed-chain count;
- rejected chain count;
- dirty-residue count;
- processing-error count;
- chain-level processing-error count;
- source/entry-level processing-error count;
- total Gold record count.

A processing-error record is chain-level only when both `model_id` and
`label_chain_id` are available. Errors without complete canonical chain
identity are counted as source/entry-level errors.

The selected-chain accounting invariant is:

`selected_silver_chain_count = accepted_chain_count + rejected_chain_count + non_candidate_chain_count + chain_level_processing_error_count`

`parsed_silver_chain_count` is deliberately not used in this equation because
the parser preserves chains from models outside the configured processing
scope.

Dirty-residue records are additional lineage and are deliberately excluded
from this equation.

`total_gold_record_count` is the physical record count across all Gold shard
types and therefore includes dirty-residue and processing-error records.

The task summary must additionally contain deterministic count mappings for:

- rejected chains by terminal reason;
- rejected chains by terminal stage;
- dirty residues by rule ID;
- dirty residues by dirty type;
- processing errors by processing stage;
- processing errors by error type.

The summary must contain explicit Boolean validation results for the
source-object and selected-chain accounting invariants. A task must not be
marked successful when either required accounting invariant is false.

Task-summary JSON must be written to a temporary path and atomically renamed
only after the summary has been fully serialized.

### 7.6.2 Deterministic Manifest Partitioning

Quality-cleaning task identifiers are zero-based integer partition indices.

For a validated immutable source manifest with `N` rows and configured
`batch_size = B`, task `t` owns the half-open manifest row interval:

`[t * B, min((t + 1) * B, N))`

Therefore:

- valid task identifiers are `0` through `ceil(N / B) - 1`;
- every manifest row belongs to exactly one task;
- no two tasks overlap;
- no source row is skipped;
- the final task may contain fewer than `B` source objects;
- an out-of-range task identifier is an execution error rather than an empty
  successful task.

Partitioning is performed against the deterministic row order of the
validated immutable manifest. The canonical manifest generator sorts rows by
`pdb_id` before publication, so rerunning the same snapshot/configuration
produces the same task-to-source assignment.

The number of quality partitions is derived at runtime from the validated
manifest row count and configured `batch_size`:

`partition_count = ceil(manifest_row_count / batch_size)`

No snapshot-specific partition count is hard-coded. When a future snapshot
contains a different number of coordinate mmCIF objects, its task count is
recomputed automatically from that snapshot's validated manifest.

### 7.7.1 Task Publication Order and Completion

A quality task must publish outputs in the following order:

1. Materialize all Gold tables and build the task summary in memory.
2. Validate both required accounting invariants.
3. Atomically write the five task-level Parquet shards:
   accepted, rejected, non-candidates, dirty residues, and processing errors.
4. Atomically write `summaries/task_<task_id>.json` last.

The task-summary file is the task-level completion marker. A task is complete
only when its five expected Parquet shards exist and its valid task-summary
JSON has been published.

If execution fails after one or more Parquet shards have been written but
before the summary is published, the task is incomplete. Rerunning the same
task must replace the deterministic task shard paths atomically rather than
append duplicate records.

The task-summary completion marker is distinct from the release-level
`_SUCCESS` marker. `_SUCCESS` is written only after all tasks and release-wide
validation have completed successfully.

Task execution timing uses the following semantics:

- `started_at_utc` is captured immediately before source processing begins;
- elapsed runtime is measured with a monotonic high-resolution clock;
- the five Parquet shards are included in task runtime;
- `completed_at_utc`, final `runtime_seconds`, and peak-memory usage are
  captured after all five Parquet shards have been successfully published
  and immediately before the task-summary completion marker is written.

Therefore, a task-summary record describes the completed processing and
Parquet-publication work that precedes that summary.

## 8. BRI Generation

BRI representations are generated only for chains accepted by the
quality-cleaning stage.

For each accepted chain, the pipeline must:

1. Retrieve the exact mmCIF object recorded in the source manifest.
2. Verify the source object against its recorded size and ETag.
3. Select the recorded model and canonical label chain.
4. Reconstruct the validated ordered N, CA, C backbone.
5. Generate the backbone rigid invariant.
6. Store the BRI vector with complete chain and source lineage.

Each BRI record must include:

- Snapshot date.
- PDB ID.
- Model ID.
- Label chain ID.
- Author chain ID.
- Chain length.
- Amino-acid sequence.
- BRI vector location or partition.
- Source mmCIF key.
- Source ETag.
- Pipeline version.
- Configuration version.

The following accounting invariant must hold:

accepted chain count
=
successful BRI count
+
explicitly recorded BRI failure count

No accepted chain may disappear silently during BRI generation.

A BRI failure should be investigated as either:

- A quality-cleaning rule that failed to detect invalid input.
- A parser or chain-selection defect.
- An unsupported BRI edge case.
- A software or infrastructure failure.

Confirmed structural failures should be added to the regression-test fixtures.

## 9. cKDTree and Strict L-Infinity Search

Chains must be partitioned by residue count before geometric comparison because
BRI vectors from chains of different lengths have different dimensions.

For each equal-length partition, the pipeline must:

1. Load the BRI vectors and their canonical chain metadata.
2. Build a `scipy.spatial.cKDTree`.
3. Query candidate pairs using radius 1.0 and Chebyshev distance.
4. Recalculate the exact L-infinity distance for every returned candidate pair.
5. Retain only pairs satisfying the strict condition:

L_inf < 1.0

The exact post-filter is required because a cKDTree radius query includes pairs
on the boundary where:

L_inf = 1.0

Singleton length groups contain no comparable pairs and are retained
automatically.

Each retained pair record must include:

- First canonical chain key.
- Second canonical chain key.
- Chain length.
- Exact L-infinity distance.
- Whether the amino-acid sequences are equal.
- Sequence difference count.
- BRI partition identifier.
- Pipeline and configuration versions.

## 10. Near-Duplicate Graph

Represent the geometric search result as a graph.

In this graph:

- Each accepted chain is a node.
- Each strict L-infinity pair is an edge.
- An edge exists only when the exact distance satisfies L_inf < 1.0.

For every connected component, record:

- Component identifier.
- Number of chains.
- Number of edges.
- Maximum possible number of edges.
- Edge density.
- Whether the component is a clique.

Connected components must not automatically be treated as complete duplicate
groups.

For example, a component may contain:

- Chain A near Chain B.
- Chain B near Chain C.
- Chain A not near Chain C.

Therefore, a chain may only be removed when it is directly within the strict
L-infinity threshold of a retained representative.

Chains with no geometric edges must remain represented as singleton components
or explicitly recorded no-edge chains.

## 11. Representative Selection

Representative selection must be deterministic, versioned, and auditable.

Candidate chains may be ranked using:

1. Coordinate completeness.
2. Experimental resolution.
3. Experimental method.
4. Alternate-location burden.
5. Sequence consistency.
6. Deposition metadata.
7. Canonical chain key as the final tie-break.

The exact ranking policy must be stored in configuration and assigned its own
version.

For clique components, one highest-ranked representative may be retained.

For non-clique components, connected-component membership alone is not
sufficient justification for removal.

A chain may only be removed when it is directly connected to a retained
representative by a strict L-infinity edge.

Every removal decision must record:

- Removed canonical chain key.
- Retained representative key.
- Exact L-infinity distance to the representative.
- Component identifier.
- Ranking evidence.
- Representative-selection policy version.

## 12. Execution Model

The pipeline uses SLURM arrays and explicit job dependencies.

### 12.1 Source Retries and Concurrency

Quality-task source handling is controlled by the versioned execution
configuration.

`execution.max_retries` is the number of retries permitted after the initial
download attempt. Therefore `max_retries: 3` permits at most four transfer
attempts for one source object.

Retries apply only to explicitly classified transient source-transfer failures
(`SnapshotTransportError`). Immutable-source verification failures, including
byte-count mismatches, missing response ETags and ETag mismatches, are not
retried. Other non-transport `SnapshotError` failures are also not retried.
Unexpected exceptions are not converted into source failures and must
propagate so that the task fails visibly.

`execution.download_concurrency` is the maximum number of source-processing
workers active within one quality task. A worker performs the source download,
verification, parsing and chain cleaning for one manifest row. This permits
network transfers to overlap while bounding task-level concurrency.

Concurrency must not change logical output ordering. Source results are
consumed and aggregated in the original validated manifest order even when
individual workers finish out of order. Gold records, processing-error lineage
and task accounting therefore remain deterministic with respect to the
immutable manifest.

The processing stages are:

1. Discover and validate a complete PDB snapshot.
2. Generate the immutable source manifest.
3. Parse and quality-filter mmCIF batches.
4. Merge and validate quality-cleaning outputs.
5. Generate BRI vectors for accepted chains.
6. Merge and validate BRI outputs.
7. Run cKDTree searches for each chain-length partition.
8. Consolidate strict L-infinity pairs.
9. Build graph components and diagnostics.
10. Select representatives.
11. Validate and publish the final release.

Each array task must write to a temporary output path first.

Example:

task_000123.parquet.tmp

After successful validation, the task atomically renames the file to:

task_000123.parquet

A stage writes a `_SUCCESS` marker only after:

- All expected tasks have completed.
- Output schemas have been validated.
- Record accounting checks have passed.
- No required shard is missing.

Downstream stages must use SLURM `afterok` dependencies so that a failed stage
prevents incomplete data from being published.

## 13. Data-Engineering Guarantees

### Idempotency

Running the same task with the same snapshot, configuration, pipeline version,
and input partition must produce the same logical output.

Completed shards must not be duplicated when a task is rerun.

### Restartability

A failed task must be rerunnable independently without repeating successful
tasks from the same stage.

### Data Lineage

Every output record must be traceable to:

- Snapshot date.
- Source S3 object.
- Source ETag.
- Pipeline Git commit.
- Configuration version.
- Processing stage.
- SLURM job and task identifiers.

### Schema Enforcement

All output schemas, data types, required fields, and nullable fields must be
defined explicitly.

The pipeline must not depend on automatic type inference when merging shards.

### Full Record Accounting

Every input structure and candidate chain must be counted.

For each stage, the following relationship must hold:

input records
=
successful output records
+
rejected records
+
explicitly recorded failures

### Observability

Each task must record:

- Start and completion timestamps.
- Input record count.
- Output record count.
- Rejection counts by rule.
- Failure counts by category.
- Runtime.
- Peak memory usage where available.
- SLURM job and task identifiers.

### Validation

A release may only be published after confirming:

- Source-manifest uniqueness.
- No duplicate canonical chain keys.
- Complete record accounting.
- Valid output schemas.
- No missing required shards.
- Successful output checksum generation.

## 14. Incremental Snapshot Processing

When a new complete PDB snapshot is released, the pipeline must compare its
manifest with the previous successful release.

Each source entry is classified as:

- Unchanged.
- New.
- Changed.
- Removed.

An entry is unchanged when its PDB ID, source key, size, and ETag match the
previous snapshot manifest.

When the pipeline version and cleaning configuration are compatible, unchanged
entries may reuse their previous quality-cleaning results.

New and changed entries must be parsed and processed again.

Entries absent from the new snapshot must not appear in the new release.

The incremental quality-cleaning relationship is:

new snapshot
=
reused unchanged entries
+
processed new entries
+
processed changed entries
-
removed entries

BRI vectors may also be reused for unchanged accepted chains when:

- The source ETag is unchanged.
- The chain-selection policy is unchanged.
- The quality-cleaning protocol is unchanged.
- The BRI implementation and configuration are unchanged.

Geometric deduplication requires additional care because new or changed chains
can create new strict L-infinity edges with previously processed chains.

Therefore, the initial implementation may reuse parsing and BRI outputs while
rebuilding the global cKDTree search, graph, and representative decisions for
each new release.

A future optimized implementation may update geometric results incrementally,
but it must produce results equivalent to a complete global recomputation.

## 15. Automatic Snapshot Watcher

A scheduled watcher may inspect the public PDB snapshot bucket for new dated
prefixes.

A dated prefix must not be treated as a complete coordinate snapshot unless:

- The divided mmCIF directory exists.
- The directory contains a plausible number of `.cif.gz` files.
- The total compressed size passes plausibility checks.
- The manifest remains stable across repeated observations.
- No completed release already exists for the snapshot.
- No active pipeline run already owns the snapshot.

The watcher must maintain persistent state containing:

- Latest complete snapshot observed.
- Latest successfully published snapshot.
- Active pipeline run, if any.
- Pipeline version.
- Configuration version.

When a valid new snapshot is detected, the watcher must:

1. Generate and validate the source manifest.
2. Compare it with the previous successful release.
3. Prepare incremental processing inputs.
4. Submit the SLURM dependency chain.
5. Record the submitted job identifiers.
6. Monitor stage-completion markers.
7. Publish the release only after all validation checks pass.

The initial 20260101 release will be triggered manually.

Automatic execution must only be enabled after the manual pipeline has
demonstrated:

- Idempotent reruns.
- Safe task restart behaviour.
- Correct SLURM dependency handling.
- Complete record accounting.
- Protection against duplicate concurrent runs.

## 16. Release Outputs

Each completed release must be stored in a versioned directory.

Example:

outputs/releases/PDBClean-20260101-protocol3.2-comp702-v1/

The release must contain:

- source_manifest.parquet
- accepted_chains.parquet
- rejected_chains.parquet
- processing_errors.parquet
- quality_summary.json
- bri_metadata.parquet
- bri_failures.parquet
- near_duplicate_pairs.parquet
- length_group_summaries.parquet
- graph_components.parquet
- retained_chains.parquet
- removed_chains.parquet
- representative_decisions.parquet
- release_summary.json
- run_manifest.json
- SHA256SUMS
- _SUCCESS

The `_SUCCESS` marker must be created last.

Its presence means that:

- All required pipeline stages completed successfully.
- Every expected output exists.
- Schemas were validated.
- Record-accounting checks passed.
- Checksums were generated.
- The release is safe for downstream use.

A release directory without `_SUCCESS` must be treated as incomplete.

## 17. Run Manifest

Every pipeline run must produce a machine-readable run manifest.

The manifest must include:

- Dataset release identifier.
- Snapshot date.
- Snapshot source prefix.
- Source-manifest checksum.
- Cleaning-protocol version.
- BRI implementation version.
- Geometric-search configuration.
- Representative-selection policy version.
- Pipeline Git commit.
- Configuration-file checksums.
- Python version.
- Important package versions.
- Start and completion timestamps.
- SLURM job identifiers.
- Output record counts.
- Validation results.
- Final run status.

The run manifest provides the primary provenance record for reproducing a
release.

## 18. Configuration Management

Scientific and operational decisions must be stored in version-controlled
configuration files rather than scattered through Python or SLURM scripts.

The initial configuration must define:

- Snapshot date.
- Snapshot source prefix.
- Selected model policy.
- Canonical chain namespace.
- Supported polymer types.
- Occupancy threshold.
- Occupancy atom scope.
- Residue-continuity policy.
- Required backbone atoms.
- Minimum consecutive backbone distance.
- Permitted amino-acid residue names.
- BRI implementation settings.
- L-infinity threshold.
- Strict or inclusive threshold behaviour.
- SLURM batch size.
- Download concurrency.
- Retry limits.
- Output locations.
- Representative-selection policy.

Every run must save a copy of its resolved configuration and its checksum.

Invalid or incomplete configuration must cause the pipeline to stop before
processing begins.

## 19. Temporary Data and Storage Policy

Compressed mmCIF files should not be duplicated permanently unless required for
a specific reproducibility decision.

During batch processing, each task should:

1. Read its assigned manifest partition.
2. Download a small bounded number of mmCIF objects.
3. Verify downloaded object size and identity.
4. Parse and process the structures.
5. Write validated output shards.
6. Remove temporary downloads after successful processing.

Temporary files must be stored under a run-specific directory.

Example:

${TMPDIR}/pdbclean/<release_id>/<job_id>/<task_id>/

Shared output directories must not be used for partially downloaded source
objects.

The pipeline must clean stale temporary files safely without deleting outputs
belonging to active jobs.

## 20. Download Reliability

Source downloads must support:

- Bounded retries.
- Exponential backoff.
- Connection timeouts.
- Size verification.
- Detection of truncated gzip files.
- Structured error reporting.
- Safe resumption at task level.

A failed download must not be classified as a structural rejection.

It must be recorded as an infrastructure or source-access failure and remain
eligible for retry.

The pipeline must distinguish among:

- Object not found.
- Network timeout.
- HTTP error.
- Size mismatch.
- Invalid gzip stream.
- Invalid mmCIF content.
- Parser failure.

## 21. Atomicity and Concurrency Control

Each output shard must be written to a unique temporary path before being
validated and atomically renamed.

Only one active run may own a given combination of:

- Snapshot.
- Pipeline version.
- Cleaning configuration.
- Release identifier.

The pipeline must use a lock or run-ownership record to prevent duplicate
concurrent submissions.

A stale lock may only be removed after confirming that its associated SLURM jobs
are no longer active.

Merge stages must reject:

- Missing shards.
- Duplicate task identifiers.
- Duplicate canonical chain keys.
- Shards produced by incompatible configurations.
- Shards produced from different source manifests.

## 22. Testing Strategy

Testing is divided into unit, integration, regression, and system tests.

### Unit Tests

Unit tests must cover:

- Snapshot-prefix validation.
- Manifest parsing.
- Canonical chain-key construction.
- Polymer-type classification.
- Occupancy filtering.
- Residue-gap detection.
- Missing backbone-atom detection.
- Duplicate backbone-atom detection.
- Backbone-distance calculation.
- Amino-acid validation.
- Exact L-infinity calculation.
- Strict removal of boundary pairs where L_inf equals 1.0.
- Deterministic representative ranking.

### Integration Tests

Integration tests must cover:

- Downloading and parsing a small mmCIF set.
- Producing accepted, rejected, and error outputs.
- Generating BRI vectors for accepted chains.
- Running cKDTree radius search.
- Building graph components.
- Selecting representatives.
- Merging task shards.
- Repeating a task without creating duplicate records.

### Regression Tests

Known structures from the previous BRI investigation must be preserved as
regression cases.

Initial examples include:

- 1AAM chain A for an internal missing residue.
- 7CLR chain FA for a multi-residue internal gap.
- 8RX0 chains for range and duplicate-source-record edge cases.
- 7XYF chains E and F for DNA rejection.
- 7XYG chains E and F for DNA rejection.
- 8GWG chain I for RNA rejection.
- 7T3M chains G, H, and I for author-label namespace mapping.
- 5CGH chain CA for source-range mismatch.
- 5U8X chain C for duplicate backbone atoms.

Each regression fixture must document:

- Expected canonical chain.
- Expected acceptance or rejection.
- Expected rule identifiers.
- Relevant structural evidence.

### System Tests

Before full production processing, the pipeline must pass:

1. A hand-selected structural fixture set.
2. A 20-entry end-to-end pilot.
3. A 500-entry manifest pilot.
4. A rerun of selected completed tasks.
5. A forced task-failure and restart test.
6. A small BRI and cKDTree equivalence test.
7. A merge and publication dry run.

## 23. Scientific Validation

The new quality-cleaning pipeline must be compared with the historical PDB707K
dataset where entries overlap.

The comparison should report:

- Historical chains present in the new snapshot.
- Historical chains absent from the new snapshot.
- Chains accepted by both pipelines.
- Chains accepted historically but rejected by the new pipeline.
- Chains rejected historically but accepted by the new pipeline.
- Differences caused by changed source structures.
- Differences caused by clarified cleaning rules.
- Differences caused by parser or namespace handling.

BRI and geometric-search validation must include:

- Successful BRI accounting.
- Comparison against the previous pairwise implementation on a small set.
- Equality of strict L-infinity pair results.
- Explicit counting of boundary pairs removed at L_inf equal to 1.0.
- Graph component and clique diagnostics.
- Direct representative-to-removed-chain distance validation.

## 24. Security and Operational Safety

The pipeline must not store credentials in source code, configuration files, job
logs, or release manifests.

Public snapshot downloads should require no embedded credentials.

Paths supplied through configuration must be validated before deletion or
cleanup operations.

Cleanup code must refuse to recursively delete:

- Empty paths.
- Filesystem roots.
- User home directories.
- Shared scratch roots.
- Paths outside the current run directory.

Log files must avoid recording unnecessary environment variables or secrets.

## 25. Repository Layout

The intended repository structure is:

src/pdbclean/
    __init__.py
    cli.py
    config.py
    schemas.py
    manifest.py
    snapshot.py
    downloader.py
    mmcif_parser.py
    chain_builder.py
    quality_filters.py
    quality_runner.py
    merge.py
    bri_extract.py
    bri_storage.py
    kdtree_search.py
    graph.py
    representatives.py
    validation.py
    audit.py
    logging_utils.py

config/pdbclean/
    protocol_3_2_comp702_v1.yaml

scripts/pdbclean/
    create_manifest.py
    run_quality_array.sbatch
    merge_quality_outputs.sbatch
    run_bri_array.sbatch
    merge_bri_outputs.sbatch
    run_kdtree_array.sbatch
    build_graph.sbatch
    publish_release.sbatch
    watch_snapshots.sbatch

tests/pdbclean/
    fixtures/
    test_manifest.py
    test_parser.py
    test_quality_filters.py
    test_bri.py
    test_kdtree.py
    test_graph.py
    test_representatives.py
    test_idempotency.py

docs/pdbclean/
    pipeline_spec.md
    protocol_decisions.md
    schema_reference.md
    operations_runbook.md

Generated outputs must not be committed to Git unless they are small,
intentional fixtures or summary reports.

## 26. Implementation Order

The pipeline will be implemented incrementally.

### Phase 1: Foundation

- Define configuration schema.
- Define accepted, rejected, non-candidate, dirty-residue, and error schemas.
- Implement structured logging.
- Implement release and canonical chain identifiers.

### Phase 2: Snapshot Manifest

- Convert the existing manifest-generation logic into reusable code.
- Validate snapshot completeness.
- Write CSV and Parquet manifests.
- Add uniqueness and checksum validation.

### Phase 3: Parser and Quality Rules

- Implement mmCIF parsing.
- Implement model and chain selection.
- Implement Q001 through Q006.
- Produce structured accepted, rejected, non-candidate, and dirty-residue records.
- Validate against known regression structures.

### Phase 4: Batch Execution

- Partition the source manifest.
- Implement download and retry handling.
- Add atomic task outputs.
- Add SLURM array scripts.
- Add task summaries and restart behaviour.

### Phase 5: BRI Integration

- Connect accepted chains to the existing BRI implementation.
- Store BRI vectors by chain-length partition.
- Record all extraction failures.
- Enforce accepted-chain accounting.

### Phase 6: Geometric Search

- Build cKDTree searches by chain length.
- Apply exact strict L-infinity post-filtering.
- Store pair and length-group summaries.
- Validate against the previous implementation.

### Phase 7: Graph and Representatives

- Build graph components.
- Calculate density and clique diagnostics.
- Implement deterministic representative selection.
- Validate every removal against a direct edge.

### Phase 8: Publication

- Merge all outputs.
- Run release-wide validation.
- Generate checksums and run manifest.
- Write `_SUCCESS`.
- Produce research summary reports.

### Phase 9: Automation

- Implement complete-snapshot discovery.
- Implement persistent watcher state.
- Add duplicate-run protection.
- Submit the SLURM dependency chain automatically.
- Test automatic execution using a simulated new snapshot.

## 27. Initial Production Acceptance Criteria

The first production release is accepted only when:

- The 20260101 source manifest contains exactly 246,905 unique mmCIF objects.
- Every source object is processed or explicitly recorded as failed.
- Every parsed Silver chain is classified as non-candidate, accepted,
  rejected, or explicitly recorded as an error.
- Every Protocol 3.2 candidate chain is accepted, rejected, or explicitly
  recorded as an error.
- No canonical chain key appears more than once in a final category.
- Q001 through Q006 have passing unit and regression tests.
- Every accepted chain is accounted for during BRI generation.
- cKDTree results match exact validation on test partitions.
- All retained geometric pairs satisfy strict L_inf < 1.0.
- No pair with L_inf equal to 1.0 is retained.
- Every removed chain has a direct strict edge to its retained representative.
- Rerunning completed tasks does not duplicate records.
- Failed tasks can be restarted independently.
- Release schemas and checksums validate successfully.
- The run manifest contains complete provenance.
- The final release directory contains `_SUCCESS`.

## 28. Immediate Next Step

After approving this specification, implementation begins with:

1. The versioned YAML configuration.
2. Explicit Python output schemas.
3. Snapshot-manifest validation code.
4. Unit tests for configuration and manifest handling.

The quality parser and structural filters will only be implemented after these
foundation components are validated.
