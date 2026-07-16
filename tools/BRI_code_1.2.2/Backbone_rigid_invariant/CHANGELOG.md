# Changelog

All notable changes to this project will be documented in this file.

## [1.2.2] - 2026-01-07

### 🚀 Features

- Add invariant BTP(bond-torsion point)
- 1 vs group invariant compare
- Add pipeline option; More info from cleaning
- Update bri pipeline
- Folder as input for bri, clean pipeline

### 🐛 Bug Fixes

- Duplicate pipeline

### 📚 Documentation

- Add docs for pipelines

### ⚡ Performance

- Pipelines; Enhance efficiency of entry loading

## [1.2.1] - 2025-06-14

### 🚀 Features

- Update pipeline for duplicate search
- Instantiate Chain from PBD files; Add residue names mapping
- Module execution

### 🐛 Bug Fixes

- Set invariant to `None` for nan-standard chain

### 📚 Documentation

- Update docs

## [1.2.0] - 2025-05-21

### 🚀 Features

- Use `biotite` to load cif and bcif files
- Update pipeline for bri computation

### 🧪 Testing

- Update structure_test and include dependency

## [1.1.2] - 2025-05-06

### 🚀 Features

- Add AlphaFold database as source
- BID and BIB show only ticks (residue id) being divisible by 10

### 🐛 Bug Fixes

- BIB uses real residue_id, shows last residue

### 📚 Documentation

- Add CHANGELOG
- Update README and pyproject info

## [1.1.0] - 2025-03-04

### 🚀 Features

- Add StructureBase to unify base func
- Add MiniChain, MiniEntry that achieve minimal func for alternative mmCIF format

### 🧪 Testing

- Add tests 

## [1.0.0] - 2024-07-09

### 🚀 Features

- Normalize pdb handling with more check. first handling main loop
- Add methods dealing with duplicated atoms and dna structures
- Process proteins with multiprocessing
- To check missing pdb. align outcome
- Flexible examine code
- Add strength on residue
- Check superposition of atoms in residues
- Change superposition check output format
- Calculate strength with different triples
- Rename class, improve entity framework
- Implement distance test
- Add max test
- Reconstruct protein coordinate demo
- Update construction method
- Pdb basis transformation with stats
- Alpha invariant compare
- Add counting residue function
- Update clean pdb data list
- Column names, recover alpha summary, 20 basic residue
- Add main_chain filter
- Simplify check functions with new decorator
- Add residue_completeness_check. rename some methods
- Add new method for residue_continuity_check
- Add entity info to PDB
- Use KNN to compare alpha-invariant pairwise
- Improve entity check
- Add HETATM param for decorator
- Update pdb entry ids
- Integrated chainwise filter
- Switch auth data to pdb-labelled data
- Entrywise integreated cleaning
- Add entity and seq_length in output
- Add angle check
- Update invariant computation
- Enable count during cleaning
- Add some comments in invariant.py
- Compute invariant for large scale dataset
- Alpha-summary updated
- Adjust invariant_NN_compare to fit new invariant
- Updated alpha-summary and compare
- Chain class for chain-wise computation and visualizaiton
- BIB and BID plot generation
- Get_angle check vector length
- Invariant_ext for dihedral angles
- Enable save of invariant and perturbation
- Improve implementation of perturbing
- Improve perturbation and lipschitz constant compute

### 🐛 Bug Fixes

- Fix iteration count
- Average_duplicated_rows
- Count CA number affected by duplicated
- Exception during superposition check
- Transform basis with chain length < 3
- Calculation branch miss
- Remove wrong pdb id
- Calc coord with orthogonal basis for 1st row only
- On_entry deals with empty dataframe
- Return None if no neighbors
- Seq_compare matches Linf_invariant compare
- Disordered atoms check
- Separate decorated methods; add computation gap for atom-atom+1
- Chain_length not in index
- Handle all HETATM entry; fix gap and clash
- Defect atoms not properly dropped in cleaning
- Pivot table construction and summary computation
- Gap & clash not removed; add sequence output
- Fix sequence extract err; add extra residues
- Length of first residue

### 🚜 Refactor

- Split chain layer calculation
- Format calculation process
- Reformat code, no longer support .pdb
- Move residue definition
- Refine structure for dist

### 📚 Documentation
- Add README.md
- Update .gitignore, add requirements
- Update documents

### ⚡ Performance

- Improve compare function
- Reduce class memory usage, add duplicated symbol
- Less columns for summary; redundant code in cleaning
