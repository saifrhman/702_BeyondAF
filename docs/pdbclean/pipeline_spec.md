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

The Silver layer stores parsed structural observations before acceptance or
rejection decisions are applied.

It contains one record per candidate chain, including:

- PDB ID.
- Model ID.
- Entity ID.
- Label chain ID.
- Author chain ID.
- Polymer type.
- Residue identifiers.
- Residue names.
- Atom observations.
- Occupancy information.
- Alternate-location information.
- Source-object lineage.

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

## 6. Quality-Cleaning Rules

All cleaning rules are configurable, versioned, and auditable.

### Q001: Protein Polymer

Retain supported polypeptide chains.

Reject:

- DNA chains.
- RNA chains.
- Carbohydrate polymers.
- Other unsupported polymer types.

For every rejected chain, record the detected polymer type and the Q001 rule
identifier.

### Q002: Full Occupancy

Reject a chain when any atom covered by the configured occupancy policy has an
occupancy value below 1.0.

For each affected chain, record:

- Minimum occupancy.
- Number of partial-occupancy atoms.
- Affected atom names.
- Residue identifiers.
- Alternate-location identifiers.
- The Q002 rule identifier.

The exact atom scope used by this rule must be defined in the versioned
configuration.

### Q003: Consecutive Residues

Observed polymer residue positions must be consecutive according to
`label_seq_id`.

Reject a chain when the observed residue sequence contains:

- Internal gaps.
- Discontinuous coordinate segments.
- Missing positions inside the retained residue range.

For each affected chain, record:

- Observed start position.
- Observed end position.
- Missing `label_seq_id` values.
- Number of missing positions.
- The Q003 rule identifier.

### Q004: Complete and Unique Backbone

Every retained residue must contain exactly one selected atom for each of:

- N.
- CA.
- C.

Reject a chain when any retained residue contains:

- A missing N, CA, or C atom.
- More than one selected N, CA, or C atom.
- Unresolved alternate-location backbone atoms.

For each affected chain, record:

- Residue identifier.
- Atom name.
- Observed atom count.
- Alternate-location identifiers.
- The Q004 rule identifier.

### Q005: Backbone-Distance Validity

Construct the ordered backbone as:

N1, CA1, C1, N2, CA2, C2, and so on.

Reject a chain when the distance between any consecutive backbone atoms is less
than 0.01 angstrom.

For each affected chain, record:

- First residue identifier.
- First atom name.
- Second residue identifier.
- Second atom name.
- Measured distance.
- The Q005 rule identifier.

### Q006: Standard Amino Acids

Retain only residues allowed by the configured amino-acid policy.

The initial protocol permits the standard 20 amino acids.

Reject a chain when it contains:

- Unsupported amino-acid residue names.
- Modified residues without an explicit approved mapping.
- Non-amino-acid residues inside the selected protein chain.

For each affected chain, record:

- Residue identifier.
- Residue name.
- Number of unsupported residues.
- The Q006 rule identifier.

Any future mapping of modified residues to standard amino acids requires a new
protocol version.

## 7. Quality-Cleaning Outputs

Each SLURM task writes independent output shards.

Accepted chains:

accepted/task_<task_id>.parquet

Rejected chains:

rejected/task_<task_id>.parquet

Processing errors:

errors/task_<task_id>.parquet

Task summaries:

summaries/task_<task_id>.json

After all tasks complete successfully, the shards are merged into:

accepted_chains.parquet

rejected_chains.parquet

processing_errors.parquet

quality_summary.json

Every rejected-chain record must include:

- Canonical chain key.
- First failing rule.
- All failing rules.
- Structured failure details.
- Source mmCIF key.
- Source ETag.

Every candidate chain must appear in exactly one final category:

- Accepted.
- Rejected.
- Processing error.

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
- Define accepted, rejected, and error schemas.
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
- Produce structured accepted and rejected records.
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
- Every candidate chain is accepted, rejected, or explicitly recorded as an
  error.
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
