# PDBClean Development Status

Last updated: 2026-08-10

Branch: `feature/pdbclean-pipeline`

Current checkpoint:

`5d64dc8` — Add 1AAM revision differential regression

Latest full test suite:

`245 passed`

---

# 1. Scientific Contract

PDBClean reproduces Protocol 3.2 from:

Olga Anosova, Alexey Gorelov, William Jeffcott, Ziqiu Jiang, Vitaliy Kurlin,
"A Complete and Bi-Continuous Invariant of Protein Backbones under Rigid Motion",
MATCH 94(1), 97-134, 2025.

The project uses the following source-of-truth hierarchy:

1. Published Protocol 3.2 defines the scientific cleaning criteria.
2. Pinned BRI v1.2.2 defines executable operational behaviour where the paper
   does not specify implementation details.
3. Minhao Wu's dissertation provides historical and downstream compatibility
   context for PDB707K.
4. Direct author correspondence provides clarification. If correspondence
   conflicts with the paper or executable implementation, the discrepancy must
   be documented rather than silently resolved.

Pinned BRI source:

`tools/backbone_rigid_invariant-1.2.2.2-comp702/src/bri/`

PDBClean must not introduce biological filters that are absent from the
published protocol or pinned BRI implementation.

---

# 2. Historical Protocol 3.2 Reference

The published Protocol 3.2 experiment used the PDB available on 4 May 2024.

Reported input:

- 213,191 PDB entries
- 1,091,420 chains

Reported Protocol 3.2 output:

- 104,688 entries
- 707,410 chains
- runtime: 4 h 48 min 11 s

Reported filtering categories:

1. 4,513 non-proteins
2. 178,153 disordered chains containing atoms with occupancy below 1
3. 201,648 chains with non-consecutive residue indices
4. 108 chains with backbone clashes below 0.01 Angstrom
5. 9,941 incomplete chains missing one of N, CA, or C
6. 4,364 chains containing non-standard amino acids

These category counts are not assumed to be mutually exclusive sequential
removals.

Final-output validation is therefore more important than attempting to sum the
individual category counts.

---

# 3. Production Snapshot

The dissertation production snapshot is:

`20260101`

Official PDB Core coordinate prefix:

`20260101/pub/pdb/data/structures/divided/mmCIF/`

Validated inventory:

- 246,905 coordinate mmCIF objects
- 85,079,649,893 compressed bytes
- approximately 79.24 GiB compressed

Real Bronze manifest location:

`outputs/pdbclean/20260101/bronze/`

Artifacts:

- `source_manifest.csv`
- `source_manifest.parquet`
- `source_manifest_summary.json`

Validated manifest properties:

- 246,905 rows
- 246,905 unique PDB IDs
- 246,905 unique S3 keys
- total bytes: 85,079,649,893

CSV SHA256:

`69b4f7e41d8ef5c7c6ebb31e1905fbb96943ec3f71f3afcc12f037d03fc91ea5`

Parquet SHA256:

`c5a79b170eb63242c6ace0b70bf3915f043a3c4c87b64e67d7b84fb9ebe494b5`

The real Bronze manifest must be retained.

The historical PDB707K chain membership is not expected to be exactly
reproduced from a newer PDB snapshot because individual PDB entries may have
been revised since the historical cleaning run.

---

# 4. Data Architecture

## 4.1 Bronze

Bronze is a persisted, snapshot-specific source manifest.

It is independent of Protocol 3.2 and can be reused by later processing stages.

Each source object is identified using immutable snapshot provenance including:

- PDB ID
- S3 object key
- expected compressed byte count
- expected ETag

Bronze is the authoritative inventory for all downstream task partitioning.

## 4.2 Silver

Silver is a deterministic in-memory parsed representation.

Silver is not published as a separate persisted dataset.

Primary objects:

- `ChainObservation`
- `AtomObservation`

The mmCIF parser preserves:

- all deposited `_atom_site` rows
- ATOM and HETATM records
- all deposited models
- label and author chain namespaces
- label and author residue namespaces
- residue names
- atom names
- occupancies
- raw occupancy tokens
- Cartesian coordinates

Scientific filtering happens after this preservation layer.

## 4.3 Gold

Gold contains Protocol 3.2 processing outcomes.

Quality-stage Gold tables are:

- accepted chains
- rejected chains
- non-candidate chains
- dirty residues
- processing errors

Gold records contain sufficient provenance to trace decisions back to the
source snapshot object and pipeline Git revision.

---

# 5. Protocol 3.2 Backbone Projection

Pinned BRI applies the equivalent of:

`get_feature("features", HETATM=False)`

before the integrated chain cleaner.

Therefore the operational Protocol 3.2 backbone input contains only rows where:

- `group_PDB == ATOM`
- atom name is one of `N`, `CA`, `C`

HETATM rows are preserved in Silver but do not enter the Protocol 3.2 cleaner.

A selected chain with no projected ATOM N/CA/C rows is classified as a
non-candidate rather than a scientific rejection or processing error.

---

# 6. Protocol 3.2 Quality Rules

## 6.1 Q001 — Protein Entry Check

Pinned BRI determines whether an entry is protein-capable at entry level.

An entry passes Q001 when any entity has an `entity_poly_type` beginning with:

`polypeptide`

This is not implemented as an independent per-chain polymer whitelist.

Mixed entries may therefore contain selected chains that are not themselves
protein polymers while the entry still passes Q001.

## 6.2 Q002 — Disorder Check

Q002 operates on the projected ATOM N/CA/C backbone rows.

Pinned BRI identifies disorder using two mechanisms.

First, duplicate backbone observations are identified using:

`(residue_id, chain_id, model_id, atom)`

Second, occupancy is inspected using the raw mmCIF occupancy token.

An occupancy is treated as defective when:

- it does not start with `1`
- and it is not `.`

Therefore:

- `1.00` passes
- `1` passes
- `.` is tolerated
- `0.50` fails
- `?` fails

There is no additional explicit altloc-rejection rule in PDBClean.

Residues affected by Q002 are removed before the next continuity check.

## 6.3 Q003 — Residue Continuity

Q003 uses `label_seq_id`.

For a projected chain with minimum and maximum label residue IDs, the expected
set is the integer range from minimum through maximum.

Missing residue ID `0` is specially ignored to match pinned BRI behaviour.

Any other internal missing label sequence ID produces a chain break.

Q003 is applied multiple times during integrated cleaning because earlier
quality rules may remove residues and thereby create a new internal gap.

## 6.4 Q004 — Backbone Completeness

Q004 operates after projection to N/CA/C.

A residue is complete when exactly three projected backbone rows remain.

Therefore a residue is defective when its backbone row count is not exactly 3.

Normally duplicate backbone atoms have already been identified by Q002 before
Q004 is reached.

## 6.5 Q005 — Backbone Clash Check

The strict clash threshold is:

`distance < 0.01 Angstrom`

A distance exactly equal to `0.01` passes.

Pinned BRI evaluates six backbone pair types:

1. N_i to CA_i
2. N_i to C_i
3. N_i to N_(i+1)
4. CA_i to C_i
5. CA_i to N_(i+1)
6. C_i to N_(i+1)

These are six pair types, not six total measurements.

For a chain of length L, the number of evaluated distances is:

`6L - 3`

The dirty residue attribution follows the first atom/residue of the tested
pair, matching pinned BRI behaviour.

## 6.6 Q006 — Standard Amino-Acid Check

Pinned BRI first maps deposited three-letter residue names through its
CCD-derived `amino_acid_short` mapping.

Examples from the pinned implementation include:

- `ALA -> A`
- `LYS -> K`
- `LLP -> K`
- `MSE -> M`
- `SEC -> U`
- `PYL -> O`
- `UNK -> X`

The mapped one-letter code is then checked against the canonical 20 amino acids:

`A R N D C Q E G H I L K M F P S T W Y V`

Therefore, under pinned BRI behaviour:

- LLP is accepted as K when it reaches Q006
- MSE is accepted as M
- SEC is non-standard
- PYL is non-standard
- UNK is non-standard

When Q006 and Q004 identify the same residue, Q006 non-standard attribution
takes precedence because of the pinned BRI dirty-residue ordering.

---

# 7. Integrated Cleaning Order

PDBClean reproduces the stateful order of pinned BRI
`integrated_chainwise_filter`.

The operational sequence is:

1. Map deposited residue labels through the pinned residue mapping.
2. Snapshot the full projected sequence.
3. Run Q002 disorder detection.
4. Remove Q002 dirty residues.
5. If empty, return the dirty lineage.
6. Run Q003 continuity.
7. If broken, reject the chain.
8. Run Q006 non-standard detection.
9. Run Q004 completeness detection on the same current chain.
10. Combine Q006 and Q004 dirty residues with Q006 precedence on overlap.
11. Remove those dirty residues.
12. If empty, return the dirty lineage.
13. Run Q003 continuity again.
14. If broken, reject the chain.
15. Run Q005 clash detection.
16. Add Q005 dirty residues.
17. Remove accumulated dirty residue IDs.
18. If empty, return the dirty lineage.
19. Run the final Q003 continuity check.
20. If broken, reject the chain.
21. Otherwise return the remaining backbone rows as clean.

Consequences:

- terminal defective residues may be removed while the remaining chain survives
- removal of an internal residue may create a chain break and reject the chain
- the cleaner is stateful; rules cannot be evaluated as independent filters


---

# 8. Model Selection and Candidate Accounting

The parser preserves all deposited models.

Production processing selects only the configured model:

    selection:
      models:
        policy: first_model
        model_id: 1

The selection policy is applied after parsing and before Protocol 3.2 quality
processing.

Only selected-model chains participate in Gold accounting.

A selected chain is a Protocol 3.2 candidate when its ATOM N/CA/C backbone
projection is non-empty.

Therefore:

- parsed chains may include models that are not selected
- selected chains may become non-candidates when their backbone projection is
  empty
- candidate chains then proceed through Q001-Q006

The selected-chain accounting identity is:

`selected Silver chains = accepted + rejected + non-candidates + chain-level processing errors`

Parsed Silver chain count is deliberately excluded from this equation because
the parser preserves all deposited models.

---

# 9. Production Configuration

Current production configuration:

`config/pdbclean/protocol_3_2_comp702_v1.yaml`

Release configuration:

    release:
      dataset_name: PDBClean
      protocol_version: protocol3.2-comp702-v1

Snapshot configuration:

    snapshot:
      mode: fixed
      snapshot_id: '20260101'
      bucket_url: https://pdbsnapshots.s3.us-west-2.amazonaws.com
      expected_mmcif_count: 246905
      expected_total_bytes: 85079649893

Execution configuration:

    execution:
      batch_size: 500
      download_concurrency: 4
      max_retries: 3
      connection_timeout_seconds: 60
      atomic_writes: true
      write_success_markers: true

Storage configuration:

    storage:
      temporary_root: ${TMPDIR}/pdbclean
      output_root: outputs/pdbclean
      retain_downloaded_mmcif: false

Configuration validation enforces valid snapshot settings, model selection,
execution settings, release provenance, and storage configuration.

---

# 10. Dynamic Task Partitioning

Task count must never be hard-coded.

The runtime partition rule is:

`partition_count = ceil(validated_manifest_row_count / configured_batch_size)`

For the current snapshot:

- manifest rows: 246,905
- batch size: 500
- resulting task count: 494

The value 494 is not a pipeline constant.

All task bounds, validation, SLURM arrays, and merge expectations must derive
dynamically from the validated Bronze manifest.

Relevant implementation:

- `manifest_partition_count`
- `select_manifest_partition`
- `resolve_manifest_snapshot`

---

# 11. Source Download and Verification

Every source object is downloaded using its validated Bronze manifest record.

The downloader verifies:

1. compressed byte count equals the manifest value
2. HTTP response ETag equals the manifest ETag

ETag is treated as immutable snapshot object identity.

It is not interpreted as a content MD5 checksum.

Snapshot exception hierarchy:

- `SnapshotError`
- `SnapshotTransportError`
- `SnapshotVerificationError`

Retry policy is deliberately narrow.

Only `SnapshotTransportError` is retried.

Deterministic verification failures are not retried, including:

- byte-count mismatch
- missing response ETag
- ETag mismatch
- invalid expected manifest metadata

Unexpected programming errors propagate instead of being converted into
ordinary source failures.

---

# 12. Concurrent Source Processing

Quality batches use:

`ThreadPoolExecutor(max_workers=download_concurrency)`

Current configured concurrency:

`4`

Concurrency applies at source-processing level:

- download
- verification
- mmCIF parsing
- model selection
- cleaning
- Gold materialization

Results are consumed in original manifest order.

Parallel execution therefore does not change deterministic Gold or error
record ordering.

---

# 13. Processing Error Taxonomy

Stable processing stages are:

## 13.1 `source_download_verify`

Used for expected source-level download or verification failures.

Examples:

- exhausted transport retries
- compressed byte-count mismatch
- missing response ETag
- ETag mismatch

These records have null model and chain identifiers.

## 13.2 `mmcif_parse`

Used when verified source bytes cannot be parsed into the required coordinate
mmCIF representation.

These are source-level failures.

## 13.3 `quality_cleaning`

Used for unexpected exceptions during processing of an individual selected
chain.

These are chain-level failures.

Exception-message strings are not parsed to determine processing stage.

---

# 14. Quality Task Publication

Quality-stage root:

`<storage.output_root>/<snapshot>/<release.protocol_version>/quality/`

Production path:

`outputs/pdbclean/20260101/protocol3.2-comp702-v1/quality/`

Each task publishes deterministic Parquet shards:

- `accepted/task_<ID>.parquet`
- `rejected/task_<ID>.parquet`
- `non_candidates/task_<ID>.parquet`
- `dirty_residues/task_<ID>.parquet`
- `errors/task_<ID>.parquet`

Each output uses an explicit Arrow schema.

Task publication order:

1. process the selected manifest partition
2. validate accounting in memory
3. atomically write all five Parquet shards
4. record completion timing and peak memory
5. write the task summary JSON last

The task summary is the task-level completion marker.

If shards exist but the summary does not exist, the task is incomplete.

A rerun replaces deterministic task paths instead of appending duplicate
records.

---

# 15. Quality Task Accounting

Each task summary records information including:

- task ID
- snapshot
- cleaning protocol
- pipeline Git commit
- start and completion timestamps
- runtime
- input source count
- successful source count
- failed source count
- parsed Silver chain count
- selected Silver chain count
- candidate entry count
- candidate chain count
- accepted chain count
- rejected chain count
- non-candidate chain count
- dirty residue count
- processing error count
- SLURM identifiers when available
- peak process memory

Important accounting checks:

- `source_object_accounting_valid`
- `selected_chain_accounting_valid`

Source-object identity:

`input sources = successful sources + failed sources`

Selected-chain identity:

`selected chains = accepted + rejected + non-candidates + chain-level processing errors`

A task must satisfy its accounting checks before publication is valid.

---

# 16. Git Provenance

Production execution requires a clean Git worktree.

`resolve_clean_git_commit()`:

- resolves the full current commit SHA
- checks Git porcelain status
- rejects modified files
- rejects staged but uncommitted files
- rejects untracked files

Generated pipeline data below `outputs/` is Git-ignored.

This ensures each production task is tied to an exact executable repository
state.

---

# 17. Quality Task CLI

Production task entry point:

`scripts/pdbclean/run_quality_task.py`

Required arguments:

- `--config`
- `--manifest`
- `--task-id`

The CLI performs:

1. configuration loading and validation
2. manifest loading
3. snapshot resolution
4. full manifest validation
5. fixed-snapshot expected count and byte validation
6. dynamic partition-count calculation
7. task partition selection
8. clean Git commit resolution
9. dynamic quality output-root construction
10. quality task execution
11. task publication

There is no hard-coded maximum task ID.

---

# 18. Real Bronze Manifest

The production Bronze manifest has already been generated successfully.

Command used:

`python scripts/pdbclean/create_manifest.py --config config/pdbclean/protocol_3_2_comp702_v1.yaml --output-dir outputs/pdbclean/20260101/bronze`

Observed result:

- snapshot: `20260101`
- coordinate mmCIF objects: 246,905
- compressed bytes: 85,079,649,893
- approximately 79.24 GiB compressed

Independent validation confirmed:

- summary row count = 246,905
- Parquet row count = 246,905
- summary total bytes = 85,079,649,893
- Parquet byte sum = 85,079,649,893
- unique PDB IDs = 246,905
- unique S3 keys = 246,905
- no temporary publication files remained

This real Bronze manifest must be retained.

---

# 19. Real End-to-End Smoke Tests

Temporary smoke runs exercised the real production path:

`manifest -> S3 -> verification -> parsing -> model selection -> cleaning -> Gold -> summary`

Temporary smoke outputs have been removed.

The production Bronze manifest remains.

## 19.1 Infrastructure Smoke

PDB IDs:

- 100D
- 101D

Results:

- input sources: 2
- successful sources: 2
- failed sources: 0
- parsed Silver chains: 11
- selected Silver chains: 11
- candidate chains: 0
- non-candidate chains: 11
- processing errors: 0

Both accounting checks passed.

This primarily validated infrastructure and non-candidate handling.

## 19.2 Protein Cleaning Smoke

PDB IDs:

- 1AAM
- 1UBQ
- 2OLO

Results:

- input sources: 3
- successful sources: 3
- failed sources: 0
- parsed Silver chains: 10
- selected Silver chains: 10
- candidate entries: 3
- candidate chains: 3
- accepted chains: 2
- rejected chains: 1
- non-candidate chains: 7
- dirty residues: 4
- processing errors: 0

Both accounting checks passed.

### 1UBQ chain A

Accepted after terminal Q002 trimming.

Original range:

`1-76`

Dirty residues:

- 73 LEU, occupancy 0.45
- 74 ARG, occupancy 0.45
- 75 GLY, occupancy 0.25
- 76 GLY, occupancy 0.25

Retained range:

`1-72`

Retained residue count:

`72`

Dirty rule:

`Q002`

`terminal_trimmed = true`

### 2OLO chain A

Accepted unchanged.

Original and retained range:

`2-394`

Residue count:

`393`

Dirty residue count:

`0`

`terminal_trimmed = false`

### 1AAM chain A

Rejected for Q003 internal continuity.

Missing projected label sequence ID:

`246`

Terminal stage:

`Q003_after_Q002`

Terminal reason:

`internal_label_seq_id_gaps:246`

The cause of this result was subsequently verified against historical PDB
revisions and pinned BRI.

---

# 20. 1AAM Revision Investigation

Historical PDB707K contains:

`1AAM,1,1,A,1,396,396,...`

Historical cleaning therefore retained all 396 residues.

Two official 1AAM revision files were compared directly.

Permanent fixtures:

`tests/pdbclean/fixtures/1aam_revisions/`

Fixture SHA256 values:

v1.3:

`bfe37d729a1feac8d4b2d8a57b7d1d8d088bd8f34d0817c5788699042ec1b012`

v2.0:

`57dd5b1c7bd03b38b1fcc3348805b553d762c3dfb0d2bd30e3aa13678a42d06f`

## 20.1 1AAM v1.3

Residue 246 is represented as:

- LYS
- `group_PDB = ATOM`
- label sequence ID 246
- author sequence ID 258

Pinned BRI:

- 396 projected residues
- residue 246 present
- accepted

PDBClean:

- 396 projected residues
- residue 246 present
- accepted

Retained residue IDs are identical.

## 20.2 1AAM v2.0

The revised entry represents position 246 as LLP.

The residue is present in the deposited mmCIF, but its coordinate rows are:

- LLP
- `group_PDB = HETATM`
- label sequence ID 246
- author sequence ID 258

The LLP record contains N, CA, and C coordinates, but they are HETATM.

Because pinned BRI uses `HETATM=False`, those coordinates do not enter the
Protocol 3.2 backbone projection.

Pinned BRI:

- projected range: 1-396
- unique projected residues: 395
- missing residue: 246
- outcome: `chain-break`

PDBClean:

- projected range: 1-396
- unique projected residues: 395
- missing residue: 246
- outcome: rejected
- terminal stage: `Q003_after_Q002`
- terminal reason: `internal_label_seq_id_gaps:246`

The projected residue IDs are identical between both implementations.

## 20.3 Interpretation

This is a PDB revision effect, not a PDBClean implementation discrepancy.

LLP maps to K in pinned BRI when it reaches Q006.

In 1AAM v2.0, LLP never reaches Q006 because its coordinates are HETATM and are
removed by the earlier backbone projection.

Therefore:

- do not special-case 1AAM
- do not alter HETATM handling to reproduce historical membership
- historical PDB707K and a modern snapshot may legitimately differ

The demonstrated principle is:

`same Protocol 3.2 implementation + different PDB revision = potentially different cleaned membership`

---

# 21. 1AAM Differential Regression

Permanent regression:

`test_1aam_revision_projection_and_outcome_match_pinned_bri`

For v1.3 it verifies:

- identical BRI and PDBClean projected residue IDs
- residue 246 present
- 396 projected residues
- both implementations accept
- identical retained residue IDs

For v2.0 it verifies:

- identical projected residue IDs
- residue 246 absent
- 395 projected residues
- pinned BRI reports chain-break
- PDBClean rejects at `Q003_after_Q002`
- missing label sequence ID is 246

The fixture SHA256 values are asserted by the test.

---

# 22. Current Quality-Stage Validation State

Latest full suite:

`245 passed`

Validation now includes:

- unit tests
- synthetic rule tests
- pinned-BRI differential tests
- residue-mapping differential tests
- threshold-boundary tests
- cleaning-order tests
- parser tests
- schema tests
- download and verification tests
- task-accounting tests
- task-publication tests
- dynamic partition tests
- Git provenance tests
- CLI tests
- real snapshot smoke tests
- real 1AAM revision differential regression

The Protocol 3.2 cleaning engine is considered sufficiently validated to move
forward.

No production cleaning behaviour should be changed merely to force equality
with historical PDB707K membership when the PDB source record itself has been
revised.


---

# 23. Historical PDB707K BRI Baseline

A previous BRI investigation was completed using:

`/users/sgsrehm1/COMP702_BeyondAF/data/checked_PDB707K_cleaned_chains_sequences_19Feb2025.csv`

Historical counts:

- total records: 707,292
- model-1 records: 461,581
- BRI attempts: 461,165
- BRI failures: 245
- BRI successes: 460,920
- comparable BRI vectors: 460,806
- singleton length groups: 114

Same-length raw pairwise comparisons:

`317,507,782`

Strict near-duplicate pairs satisfying:

`L-infinity < 1`

count:

`5,258,576`

cKDTree candidate pairs using radius:

`L-infinity <= 1`

count:

`5,264,356`

Pairs lying exactly on the boundary:

`L-infinity = 1`

count:

`5,780`

These boundary pairs must be removed.

A pilot comparison between the historical pairwise implementation and the
cKDTree implementation produced:

- pairwise strict pairs: 1,007
- cKDTree strict pairs after exact filtering: 1,007

This established equality for the pilot subset.

---

# 24. Near-Duplicate Distance Semantics

Near-duplicate detection uses Chebyshev distance on comparable BRI vectors.

The scientific threshold is strict:

`L-infinity < 1`

The cKDTree stage may use radius:

`<= 1`

only as a candidate-generation step.

Every returned candidate pair must then undergo an exact L-infinity
calculation.

Pairs satisfying exactly:

`L-infinity = 1`

must be removed.

Therefore the operational sequence is:

1. group comparable chains by compatible BRI shape or length
2. build cKDTree using Chebyshev geometry
3. retrieve radius-1 candidates
4. calculate exact L-infinity for every candidate
5. retain only pairs with strict distance below 1

The cKDTree result itself is not the final scientific duplicate set.

---

# 25. BRI Generation Gate

BRI generation must occur only after the Protocol 3.2 Gold accepted-chain
dataset has been produced and globally validated.

Before using cKDTree, the new BRI-generation implementation must be compared
against pinned BRI on the exact same retained backbone coordinates.

The differential gate must verify:

- identical retained chain length
- identical residue order
- identical atom order expected by BRI
- identical BRI matrix dimensions
- tight numerical equality of BRI values
- identical compact-vector construction
- identical Chebyshev distance behaviour
- identical strict threshold semantics

No large-scale cKDTree search should begin until this gate passes.

---

# 26. Historical Duplicate Graph Results

The historical strict duplicate graph contained:

- chains with at least one duplicate edge: 279,910
- chains with no duplicate edge: 181,010
- connected components with edges: 45,492

Component structure:

- clique components: 36,472
- non-clique components: 9,020

A provisional one-per-connected-component removal would have removed:

`234,418`

and retained:

`226,502`

However, this strategy is not scientifically safe for non-clique components.

Connectedness is transitive as a graph property.

The strict near-duplicate relation itself is not automatically transitive.

Therefore connected components must not be treated as duplicate equivalence
classes without additional checks.

---

# 27. Representative Selection Safety Rule

Representative selection must be deterministic.

More importantly, every removed chain must be directly justified.

Required safety rule:

`every removed chain must have a direct L-infinity < 1 edge to its retained representative`

It is not sufficient for two chains merely to belong to the same connected
component.

For example, if:

`A -- B -- C`

but A and C do not satisfy the strict duplicate threshold directly, retaining A
does not justify removing C solely because C is connected to A through B.

Representative selection must therefore preserve direct-edge evidence for
every removal.

The ranking policy itself must also be deterministic so that repeated runs over
the same immutable inputs produce the same retained set.

---

# 28. Expected Final Auditability

The final dataset must support reconstruction of every major decision.

For quality cleaning, audit evidence includes:

- source snapshot object
- model and chain identity
- acceptance or rejection
- terminal rule or stage
- dirty residues
- decision-time evidence
- retained residue range
- retained sequence where applicable

For duplicate cleaning, audit evidence must include:

- both chain identifiers for each strict duplicate edge
- exact L-infinity distance
- representative identifier
- deterministic representative rank information
- direct representative-to-removed edge
- retained or removed status
- component or local graph context where useful

The final release must make it possible to answer:

`Why was this chain removed?`

without rerunning the entire pipeline.

---

# 29. Remaining Pipeline Stages

The Protocol 3.2 cleaning engine is sufficiently validated.

The remaining stages are listed below in execution order.

## 29.1 Quality Shard Merge and Global Validation

Implement a stage that discovers all expected quality task outputs and validates
the complete quality run.

Requirements:

- derive expected task count dynamically from the validated Bronze manifest
- discover task summaries
- detect missing tasks
- detect unexpected task IDs
- detect duplicate or conflicting task publication
- verify all five expected shards per completed task
- validate Parquet schemas
- validate task summary schema
- revalidate task-level accounting
- validate task snapshot
- validate cleaning protocol
- validate pipeline Git commit consistency
- aggregate global counts
- detect duplicate Gold chain identities
- validate accepted/rejected/non-candidate/error disjointness where applicable
- produce deterministic merged outputs
- produce a global validation summary
- write stage-level `_SUCCESS` only after all validation passes

A partial task output without its summary must not count as complete.

Rerunning the merge stage over unchanged valid task outputs must be idempotent.

## 29.2 Dynamic SLURM Quality Array

Create production SLURM execution for the quality stage.

Requirements:

- derive partition count from the actual validated manifest
- never hard-code 494
- submit one task per dynamic manifest partition
- pass task ID to `run_quality_task.py`
- use the production configuration and production Bronze manifest
- write logs to deterministic locations
- preserve task rerun capability
- avoid node-specific hard-coding
- use appropriate Barkla resources

The full 20260101 quality run should begin only when merge/global validation is
ready to consume and verify its outputs.

## 29.3 Full 20260101 Protocol 3.2 Run

Process all coordinate mmCIF objects in the fixed snapshot.

After completion:

- verify every expected task summary exists
- run the global merge/validator
- inspect failure counts and categories
- investigate unexpected processing errors
- publish quality-stage `_SUCCESS` only after validation

## 29.4 BRI Generation

Generate BRI representations from accepted Gold chains only.

BRI input must use the exact retained Protocol 3.2 backbone.

Record sufficient chain identity and provenance to join BRI vectors back to
Gold accepted chains.

## 29.5 BRI Differential Validation

Compare the new BRI generation against pinned BRI.

Do not proceed to production near-duplicate search until numerical and semantic
equivalence is established.

## 29.6 cKDTree Candidate Search

Group comparable BRI vectors.

Use cKDTree Chebyshev radius 1 to generate candidates.

Then perform exact strict filtering.

Final duplicate edges must satisfy:

`L-infinity < 1`

## 29.7 Duplicate Graph Construction

Build an undirected graph from strict near-duplicate edges.

Record:

- node identifiers
- edge endpoints
- exact distances
- component identifiers
- graph statistics

## 29.8 Deterministic Representative Selection

Apply a documented deterministic ranking policy.

Every removed chain must have a direct strict duplicate edge to its chosen
representative.

Non-clique connected components require particular care.

## 29.9 Final Dataset and Audit Tables

Produce:

- retained chains
- removed chains
- strict duplicate edges
- representative mapping
- component summaries
- quality-cleaning audit tables
- processing-error tables
- relevant provenance

## 29.10 Final Release

Target release structure:

`outputs/releases/PDBClean-<snapshot>-<protocol>/`

The final release should include:

- retained cleaned dataset
- removed-chain audit
- strict duplicate edges
- representative mapping
- source manifest
- pipeline configuration
- Git provenance
- validation summaries
- relevant checksums
- final `_SUCCESS`

The exact release layout must remain deterministic and versioned.

---

# 30. Immediate Next Engineering Task

The next implementation task is:

`quality shard merge and global validation`

Do not start the full 20260101 SLURM array before this stage exists.

The merge/validator should be designed first because it defines the completion
contract for the expensive distributed cleaning run.

The first design questions to resolve are:

1. exact merged output paths
2. expected task-summary discovery rules
3. task and shard validation rules
4. global accounting identities
5. duplicate-record detection
6. Git/config/snapshot consistency rules
7. stage-level summary schema
8. `_SUCCESS` publication semantics
9. rerun and partial-output behaviour

---

# 31. Development Working Rules

The following rules must remain in force during development.

- Do not hard-code snapshot-specific task counts.
- Do not add biological filters absent from Protocol 3.2 or pinned BRI.
- Do not special-case individual PDB entries.
- Preserve immutable source provenance.
- Distinguish paper semantics from executable BRI behaviour.
- Explicitly document any discrepancy between paper, code, dissertation, and
  author correspondence.
- Preserve all deposited data in Silver before scientific filtering.
- Keep Silver in memory unless architecture is deliberately revised.
- Use explicit Arrow schemas for published tabular outputs.
- Write deterministic output paths.
- Use atomic publication.
- Treat summaries and `_SUCCESS` markers as explicit completion contracts.
- Preserve idempotent rerun behaviour.
- Use explicit Git staging paths.
- Never use `git add -A`.
- Production runs require a clean Git worktree.
- Keep the real Bronze manifest.
- Temporary smoke and revision-check outputs may be removed after validation.
- Do not modify the cleaning implementation merely to match historical PDB707K
  membership when the source PDB revision differs.
- Do not move to cKDTree until BRI equivalence has been established.
- Do not treat connected components automatically as duplicate equivalence
  classes.
- Every removed duplicate chain must have a direct strict edge to its retained
  representative.

---

# 32. Current Checkpoint

Current Git checkpoint before this status document:

`5d64dc8`

Commit message:

`Add 1AAM revision differential regression`

Latest validated full test suite:

`245 passed`

The next code-development stage is the quality shard merge/global validator.

---

# 33. Barkla Long-Partition QOS Limits

Production submission testing identified an additional Barkla scheduling
constraint applied by the `long` partition.

Partition configuration:

- allowed user QOS: `normal`
- partition QOS: `longlimits`

The effective `longlimits` values are:

- `MaxSubmitJobsPU = 100`
- `MaxJobsPU = 4`

A single Slurm array containing all 494 logical quality partitions was therefore
rejected with:

`QOSMaxSubmitJobPerUserLimit`

The Slurm array concurrency suffix does not solve this issue because an array
such as `0-493%24` still represents 494 submitted array elements.

The pipeline must therefore distinguish:

- logical PDBClean task IDs, derived dynamically from the Bronze manifest
- physical Slurm worker IDs, constrained by the cluster scheduler

The planned production mapping preserves all 494 logical task IDs and their
deterministic `task_<id>` outputs while distributing them across a smaller
number of physical Slurm workers.

Current planned physical execution:

- up to 64 Slurm array workers
- physical array concurrency capped at 4
- each worker processes logical task IDs using deterministic striding
- final merge remains dependent on successful completion of all physical
  workers

For the current 20260101 manifest this implies approximately 7-8 logical
quality partitions per physical worker.

This scheduling adaptation changes only orchestration. It does not change
Protocol 3.2 semantics, manifest partitioning, Gold schemas, logical task
identities, or scientific outputs.

---

# 34. Snapshot Selection Is Not Hard-Coded

The production PDBClean configuration is no longer tied to the
`20260101` PDB snapshot.

The generic production configuration now uses:

- `snapshot.mode: latest_complete`
- the official PDB snapshots bucket
- no fixed `snapshot_id`
- no snapshot-specific expected mmCIF count
- no snapshot-specific expected byte total

The existing snapshot subsystem already supports both:

- `fixed`: explicitly select a concrete `YYYYMMDD` snapshot
- `latest_complete`: discover available snapshots and resolve the newest
  complete coordinate snapshot

A dynamically selected snapshot is resolved to a concrete snapshot ID before
dataset construction. That concrete ID is recorded in the Bronze manifest and
is then used by downstream quality/merge stages, preserving reproducibility.

The existing `20260101` Bronze manifest remains a valid fixed historical
dataset and must not be deleted or modified merely because the generic
production configuration now uses `latest_complete`.

Additional configuration validation now prevents `latest_complete` from
carrying:

- `snapshot_id`
- `expected_mmcif_count`
- `expected_total_bytes`

This prevents snapshot-specific constants from being accidentally reintroduced
into the generic production configuration.

Regression status after this change:

- targeted snapshot/manifest/quality/merge/submission tests: 75 passed
- full repository test suite: 297 passed

---

# 35. Latest-Snapshot Resolution Must Avoid Recursive Crawls

Testing `snapshot.mode: latest_complete` against the live PDB snapshots bucket
identified a performance problem in the original resolver.

Snapshot discovery returned 33 dated prefixes, with the newest including:

- 20260415
- 20260101
- 20250101
- 20240101

The `20260415` prefix does not expose the canonical divided-mmCIF coordinate
archive expected for a complete PDB coordinate snapshot.

Previously, when the canonical check failed, `latest_complete` fell back to a
recursive S3 scan of the entire dated prefix. This caused snapshot resolution
to stall for a long time.

The resolver was corrected so that:

- `latest_complete` checks the canonical divided-mmCIF layout only;
- a non-canonical candidate is skipped immediately;
- the resolver then tries the next discovered dated snapshot;
- `fixed` mode retains the recursive-layout fallback for explicitly requested
  unusual snapshot layouts.

Live read-only verification after the fix resolved:

- selection mode: `latest_complete`
- concrete snapshot: `20260101`
- layout: `canonical_divided_mmcif`
- source prefix:
  `20260101/pub/pdb/data/structures/divided/mmCIF/`

Thus `20260101` is now the result of dynamic snapshot discovery, not a
hard-coded production snapshot.

Regression status after the resolver fix:

- snapshot tests: 19 passed
- full repository suite: 297 passed

---

# 36. Production Quality-Worker Memory Requirement

The first complete production quality run on the resolved `20260101` snapshot
showed that the original 8 GiB Slurm memory request was insufficient for a
small subset of physical workers.

Production array job:

- Slurm job: `10186820`
- physical workers: `0-63`
- concurrency: 4
- requested memory per worker: 8 GiB

Of the 64 physical workers, 59 completed successfully and five were terminated
for out-of-memory conditions:

- worker 20
- worker 24
- worker 31
- worker 37
- worker 41

Those five workers left seven logical quality tasks incomplete:

- 415
- 425
- 468
- 472
- 479
- 485
- 489

The seven missing logical tasks were resubmitted in recovery job `10193699`
with a 16 GiB memory request. All seven completed successfully with exit code
`0:0`.

After recovery, quality-task coverage was validated as:

- completed logical tasks: 494
- missing logical tasks: 0
- unexpected logical tasks: 0

The production quality-worker Slurm request has therefore been increased from:

`#SBATCH --mem=8G`

to:

`#SBATCH --mem=16G`

This is an execution-resource correction only. It does not change the Bronze
manifest, logical task partitioning, Protocol 3.2 cleaning semantics, Gold
schemas, provenance identities, or scientific results.

A regression test now requires the production quality worker to request
16 GiB, preventing an accidental return to the insufficient 8 GiB setting.

Validation after the change:

- Slurm submission/worker tests: 10 passed
- full repository test suite: 304 passed
