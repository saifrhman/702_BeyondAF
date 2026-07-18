4BP8 and 7DSQ epoch-32 inference provenance

Inference jobs:
- 4BP8: Slurm job 9827077
- 7DSQ: Slurm job 9827080

Each inference run produced:
- 100 independently seeded unrelaxed predictions
- 100 corresponding relaxed structures

Seeds:
- 1001 through 1100

Checkpoint:
- OpenFold epoch-32 checkpoint 32-330000.ckpt

Invariant analysis jobs:
- 4BP8-A: Slurm job 9834364
- 7DSQ-B: Slurm job 9834365

BRI implementation:
- backbone-rigid-invariant 1.2.2.2

Important chain-assignment finding:
The inherited FASTA header and prediction filenames identify 7DSQ as
7dsq_A. Sequence analysis showed that the 464-residue sequence instead
matches experimental 7DSQ chain B with 100% sequence identity and 100%
coverage. Experimental comparisons therefore use 7DSQ chain B.

The OpenFold runs used 224 downloaded template mmCIF entries. Their exact
PDB identifiers are preserved in template_ids.txt. The approximately
396 MB externally retrievable RCSB template cache itself is not duplicated
in this Git repository.
