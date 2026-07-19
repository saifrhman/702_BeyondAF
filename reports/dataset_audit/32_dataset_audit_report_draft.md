# COMP702 Dataset Audit Report

## 1. Scope and objective

This audit examines the datasets and data-processing steps used in the inherited OpenFold retraining workflow for the COMP702 project.

The main objectives are to:

1. Compare the PDB707K dataset with the OpenFold data used in the inherited workflow.
2. Determine what data entered the OpenFold retraining pipeline.
3. Reconstruct, where possible, how the retraining subset was produced.
4. Identify discrepancies between reported dataset counts and the files that are currently preserved.
5. Document any missing provenance that prevents exact reproduction of the original dataset-selection procedure.

The audit is based only on material transferred from the previous project work and subsequently centralised under `COMP702_BeyondAF`, together with the preserved PDB707K dataset, inherited OpenFold scripts and notebooks, training configuration files, and available Git history.

Historical paths referring to `/users/sgmwu14/...` are treated only as provenance describing where data originally existed. Those directories are not accessible from the current account and are not used as live data sources.

## 2. PDB707K dataset audit

### 2.1 Preserved dataset

The preserved dataset is:

`data/checked_PDB707K_cleaned_chains_sequences_19Feb2025.csv`

Direct verification gives:

- Total physical lines: 707,293
- Header lines: 1
- Data rows: 707,292

The dataset contains the following columns:

`pdb_id, entity_id, model_id, chain_id, start_residue, chain_length, auth_chain_length, auth_chain_id, auth_seq_id_start, auth_seq_id_end, seq`

### 2.2 Difference from the original reported PDB707K dataset

The inherited BRI notebook shows that the original file
`PDB707K_cleaned_chains_sequences_19Feb2025.csv`
contained 707,410 rows.

The same notebook reconstructs the creation of the preserved checked dataset:

- Original dataset: 707,410 rows
- After removal of rows associated with the duplicated/problematic subset: 707,174 rows
- Corrected/reconstructed rows added back: 118
- Final checked dataset: 707,292 rows

Therefore:

`707,410 - 236 + 118 = 707,292`

The previously observed net difference of 118 rows is therefore explained by a correction workflow in which 236 original rows were removed and 118 corrected/reconstructed rows were restored.

The corrected replacement dataframe was produced by inner-merging `af_clean.csv` with `duplicated_PDB707K_cleaned_chains_sequences_19Feb2025.csv` using:

- `pdb_id`
- `model_id`
- `chain_id`
- `chain_length`

The intermediate files `af_clean.csv`,
`duplicated_PDB707K_cleaned_chains_sequences_19Feb2025.csv`, and
`checked_duplicated_PDB707K_cleaned_chains_sequences_19Feb2025.csv`
are not present in the transferred material. Therefore, the numerical transformation can be reconstructed from notebook code and saved outputs, but the exact identities of all 236 removed rows and 118 replacement rows cannot currently be enumerated.

### 2.3 Model-1 chain preparation

The preserved preparation script filters the checked PDB707K dataset to `model_id = 1`, normalises each chain identifier as:

`lowercase_pdb_ID + "_" + uppercase_chain_ID`

and then removes duplicate normalised chain keys.

Direct computation gives:

- All checked PDB707K rows: 707,292
- `model_id = 1` rows: 461,581
- Unique normalised model-1 chain keys: 461,578
- Duplicate model-1 chain-key rows removed: 3
- Unique PDB entries represented among model-1 rows: 104,586

The current prepared PDB707K chain list therefore contains 461,578 unique model-1 chain identifiers.

## 3. Comparison with the OpenFold data

### 3.1 Reported OpenFold and PDB707K comparison

The inherited project documentation reports the following counts for the comparison between processed OpenFold identifiers and PDB707K:

- Processed OpenFold alignment identifiers: 130,413
- Reported PDB707K/OpenFold intersection: 65,656
- Reported PDB707K-only chains: 641,754
- Reported OpenFold-only chains: 67,181

The comparison is described as operating on normalised chain identifiers and using `model_id = 1` for PDB707K.

The preserved PDB707K preparation script uses identifiers of the form:

`lowercase_pdb_ID + "_" + uppercase_chain_ID`

However, the original script that performed the OpenFold/PDB707K matching has not been located in the transferred material.

### 3.2 Unresolved OpenFold-only count discrepancy

The reported PDB707K side is arithmetically consistent:

`707,410 - 65,656 = 641,754`

This exactly matches the reported PDB707K-only count.

In contrast, the reported OpenFold counts do not correspond to a simple unique-set subtraction:

`130,413 - 65,656 = 64,757`

This differs from the reported OpenFold-only count of 67,181 by 2,424.

Therefore, the four reported values cannot currently be reproduced as a simple two-set comparison using unique chain identifiers.

Possible explanations include:

- different processing stages being counted;
- duplicate-chain handling;
- identifier normalisation or transformation;
- additional filtering;
- repeated identifiers;
- or a discrepancy in the reported counts.

The transferred material does not contain the original 130,413-identifier OpenFold list, the exact 65,656-chain intersection list, or the original matching/filtering script. Consequently, the source of the 2,424-count difference cannot currently be resolved exactly.

### 3.3 Reproducibility status of the comparison

The PDB707K side of the comparison can be reproduced from the preserved checked dataset and preparation script.

The OpenFold side cannot currently be reproduced chain-for-chain because the following original intermediate artefacts were not transferred or are not present in accessible storage:

- the list of 130,413 processed OpenFold identifiers;
- the exact list of 65,656 reported intersection chains;
- the original OpenFold/PDB707K matching script;
- the original `alignment_subset/alignments` directory;
- the associated `lists/` directory used in the subset-packaging workflow.

Therefore, the reported comparison can currently be documented and partially reconstructed, but its exact chain-level membership cannot be independently verified from the preserved artefacts.

## 4. OpenFold retraining pipeline

### 4.1 Retraining inputs identified from inherited scripts

The inherited H100 training script associated with the retraining run references the following inputs:

- Structural mmCIF files:
  `pdb_data/mmcif/mmcif_files`

- Alignment database:
  `alignment_subset/alignment_dbs`

- Alignment database index:
  `alignment_subset/alignment_dbs/alignment_db.index`

- Chain metadata cache:
  `pdb_data/data_caches/chain_data_cache.json`

- mmCIF/template metadata cache:
  `pdb_data/data_caches/mmcif_cache.json`

The training command used the OpenFold `initial_training` configuration with:

- Maximum template date: `2021-10-10`
- Random seed: `42`
- Precision: `bf16-mixed`
- Hardware: one H100 GPU in the identified H100 training script

The output directory referenced by this script was:

`openfold_runs/run_subset_gpu32_h100`

This corresponds to the retraining run associated with the inherited epoch-32/33 checkpoint material.

### 4.2 Reconstruction of the downstream subset-processing pipeline

The inherited OpenFold notebook shows that the training pipeline started from an already-created directory:

`alignment_subset/alignments`

The subsequent processing was:

`alignment_subset/alignments`
→ `create_alignment_db_sharded.py`
→ `alignment_subset/alignment_dbs/alignment_db.index`
→ `alignment_data_to_fasta.py`
→ `all-seqs.fasta`
→ MMseqs clustering at 40% sequence identity
→ `all-seqs_clusters-40.txt`
→ `generate_chain_data_cache.py`
→ training cache files
→ OpenFold retraining

The 40% sequence-identity cluster file was therefore generated after the alignment subset already existed. It does not explain how the original subset itself was selected.

### 4.3 Missing upstream subset-selection step

The transferred scripts and notebook consistently consume an already-existing:

`alignment_subset/alignments`

No preserved script has been found that creates this directory by matching the full OpenFold alignment data against PDB707K.

The packaging scripts show that the original workflow packaged:

- `alignment_subset`
- `lists`
- `pdb_data`

However, neither the original `alignment_subset`, the `lists` directory, nor the packaged subset archive is present in accessible storage.

Therefore, the downstream retraining pipeline can be reconstructed with reasonable confidence, but the upstream step that selected the exact chains placed into `alignment_subset/alignments` is not currently reproducible from the transferred artefacts.

## 5. OpenFold duplicate-chain handling

OpenFold uses the auxiliary file:

`duplicate_pdb_chains.txt`

The inherited OpenFold documentation states that this file groups PDB chains that are sequence-identical.

The inherited OpenFold source code in `create_alignment_db_sharded.py` shows how these groups are handled when constructing the alignment database index.

For each duplicate group, OpenFold:

1. Searches the chains in the group for the first chain that already has alignment data in the database index.
2. Uses that chain as the alignment representative.
3. Maps the remaining chain identifiers in the duplicate group to the same underlying alignment database entry.

Therefore, duplicate chain identifiers are not necessarily removed from the alignment index. Instead, sequence-identical chains can share the same underlying alignment-data representation, reducing redundant alignment storage.

This is distinct from the proposed geometric de-duplication stage:

- OpenFold duplicate handling: based on sequence-identical chains and shared alignment data.
- Proposed COMP702 geometric de-duplication: based on geometric similarity measured using BRI invariants and L-infinity distances.

The OpenFold duplicate-chain mapping therefore provides a useful sequence-based redundancy baseline for the later comparison with BRI-derived geometric near-duplicate groups.

The original `duplicate_pdb_chains.txt` file used in the inherited run is not present in the transferred material, although the OpenFold source code and documentation describing its role are preserved.

## 6. Audit conclusions

### 6.1 Comparison of the OpenFold data and the data used in the inherited project

The preserved evidence indicates that the inherited project did not retrain OpenFold directly on the full available OpenFold alignment collection.

Instead, the workflow used a pre-existing subset stored under:

`alignment_subset/alignments`

This subset was subsequently converted into an OpenFold alignment database, FASTA file, 40% sequence-identity clusters, and training cache files before retraining.

The project documentation reports that 130,413 processed OpenFold alignment identifiers were compared with PDB707K and that 65,656 identifiers formed an intersection. However, the original identifier lists and the script that performed this comparison are not present in the transferred material.

The PDB707K dataset used in the current project is a checked version containing 707,292 rows. The inherited BRI notebook shows that it was derived from an original 707,410-row PDB707K dataset through a correction workflow in which 236 rows were removed and 118 corrected/reconstructed rows were added back.

Therefore, the current preserved PDB707K dataset differs from the original reported PDB707K dataset by a documented data-correction step rather than by the current model-1 filtering script.

### 6.2 What entered the OpenFold retraining pipeline

The retraining pipeline used:

- structural coordinates from the OpenFold/PDB mmCIF dataset;
- chain metadata from `chain_data_cache.json`;
- template metadata from `mmcif_cache.json`;
- precomputed alignment data from `alignment_subset/alignment_dbs`;
- the corresponding `alignment_db.index`;
- sequence clusters generated at 40% sequence identity from the selected alignment subset.

The critical chain-selection step occurred before `alignment_subset/alignments` entered the preserved notebook workflow.

The exact chain membership of that original subset cannot currently be reconstructed because the original alignment subset, chain lists, and matching script were not included in the transferred artefacts.

### 6.3 Overall reproducibility assessment

The following parts of the inherited dataset workflow are reproducible or reconstructable from the preserved material:

- the current checked PDB707K dataset and its row count;
- the PDB707K model-1 filtering and chain-ID normalisation;
- the numerical transformation from the original 707,410-row PDB707K dataset to the checked 707,292-row dataset;
- the downstream OpenFold alignment-database construction;
- the 40% sequence clustering stage;
- the generation of OpenFold training caches;
- the retraining input paths and major training configuration;
- OpenFold's handling of sequence-identical duplicate chains.

The following parts are not currently reproducible at chain level:

- the original 130,413 processed OpenFold identifier list;
- the reported 65,656-chain PDB707K/OpenFold intersection;
- the reported 67,181 OpenFold-only set;
- the exact script used to match OpenFold identifiers against PDB707K;
- the exact contents of the original `alignment_subset/alignments`;
- the original `lists/` directory used during subset packaging.

Consequently, the downstream OpenFold retraining pipeline is substantially reconstructable, but the exact upstream procedure used to define the training subset is not fully reproducible from the transferred artefacts.

The next dataset stage should therefore treat the preserved 461,578 unique model-1 PDB707K chains as the reproducible starting point for the new geometric de-duplication workflow, while clearly separating this new dataset construction from the incompletely preserved historical subset-selection procedure.

## 7. Evidence and provenance

The main preserved artefacts used in this audit were:

- `data/checked_PDB707K_cleaned_chains_sequences_19Feb2025.csv`
  - Preserved checked PDB707K dataset.

- `scripts/01_prepare_pdb707k.py`
  - Reproduces model-1 filtering, chain-ID normalisation, and chain-key de-duplication.

- `tools/BRI_code_1.2.2/Backbone_rigid_invariant/showcase.ipynb`
  - Preserves the workflow and saved outputs showing the transformation from the original 707,410-row PDB707K dataset to the 707,292-row checked dataset.

- `code/COMP390_code/docs/OpenFold_Documentation_Notebook.ipynb`
  - Documents the downstream processing of `alignment_subset/alignments`, including alignment-database creation, FASTA generation, 40% sequence clustering, and cache generation.

- `code/COMP390_code/task_scripts_original_backups/Task_02_Training_sbatch/openfold_gpu32_h100.sbatch.orig_minhao_paths`
  - Identifies the structural, alignment, cache, and configuration inputs used for the H100 OpenFold retraining run.

- `code/COMP390_code/task_scripts_original_backups/Task_04_Packing_sbatch/pack_openfold_subset_to_scratch.sbatch`
  - Shows that the historical workflow packaged `alignment_subset`, `lists`, and `pdb_data`.

- `code/COMP390_code/openfold/scripts/alignment_db_scripts/create_alignment_db_sharded.py`
  - Shows the exact OpenFold behaviour for mapping sequence-identical duplicate chains to shared alignment database entries.

- `code/COMP390_code/openfold/docs/source/Aux_seq_files.md`
  - Documents the meaning and format of `duplicate_pdb_chains.txt`.

Git history and searches across the transferred and centralised project material were also inspected for the original OpenFold/PDB707K matching script and chain lists. No preserved copy of those upstream selection artefacts was located.
