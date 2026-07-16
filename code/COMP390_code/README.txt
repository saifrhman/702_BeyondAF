COMP390 Project Artifact
Project title: AlphaFold and Other Machine Learning Methods for Proteins

Author: Minhao Wu
Module: COMP390 Honours Year Project

Overview
This folder contains the code, selected checkpoints, reduced demo inputs, representative analysis outputs, documentation files and environment information for the COMP390 dissertation project.

The project studies OpenFold/AlphaFold-style protein structure prediction, repeated inference, pLDDT and relaxation behaviour, and geometry-aware evaluation using BRI/LAI/BTI descriptors.

This archive is intended to support the dissertation and to provide a record of the implemented workflow. It is not a complete copy of the full HPC working directory. Complete OpenFold training resources, full alignment databases, raw downloaded datasets, temporary cache directories and all intermediate HPC outputs are not included because they are too large for submission.

Quick start for the runnable demo
The runnable part of this archive is the reduced single-target 2OLO inference demo.

For detailed commands, please read:

docs/INFERENCE_DEMO_README.txt

Alternatively, the notebook version of the same demo is:

COMP390_inference_demo_fixed.ipynb

The demo uses:

demo_single_2olo/

and compares:
1. AlphaFold parameter file
2. OpenFold official checkpoint
3. Project-trained checkpoint

The reduced demo has been tested on the slimmed archive and produced both unrelaxed and relaxed PDB outputs.

Reproducibility scope
The archive supports three main purposes:

1. inspection of the task scripts and analysis workflow;
2. inspection of selected checkpoints, representative outputs and documentation;
3. running a reduced single-target inference demo for 2OLO.

The documentation notebook records the broader project workflow, including setup, data preparation, training, inference, relaxation analysis, BRI/LAI/BTI analysis and visualisation. Some commands in that notebook contain historical Barkla HPC paths or depend on external resources that are not included in this reduced archive. Therefore, not every command in the documentation notebook is expected to run directly from the submitted archive.

For the runnable inference example, use:
- docs/INFERENCE_DEMO_README.txt
- COMP390_inference_demo_fixed.ipynb
- demo_single_2olo/

Main folders
1. task_scripts/
   Task-based scripts used during the project, including data preparation, training, inference, packing, pLDDT analysis, relaxation analysis and display-selection workflows.

2. checkpoints/
   Selected project checkpoint and training-summary files, including epoch-level metrics and the checkpoint used for the reduced inference demo.

3. demo_examples/
   Selected target-level example outputs from AlphaFold/OpenFold-style runs. These are retained as representative examples rather than as a complete archive of all runs.

4. demo_single_2olo/
   A reduced single-target inference demo for 2OLO. This folder contains:
   - fasta_dir/
   - taskcif/
   - alignments_cf/

   This is the main runnable input folder for the minimal inference example.

5. openfold/
   Local OpenFold code snapshot used for training and inference in this project. This is included because the reduced inference workflow depends directly on the local OpenFold installation.

6. openfold_inference/
   Project-level inference input directories, including FASTA files and precomputed alignments used during the project.

7. environment/
   Exported software environment information:
   - openfold_env_full.yml
   - openfold_env_from_history.yml
   - openfold_env_pip_freeze.txt
   - system_info.txt

8. docs/
   Documentation and artifact notes, including:
   - README.txt
   - INFERENCE_DEMO_README.txt
   - OpenFold_Documentation_Notebook.ipynb
   - openfold_metrics_data_dictionary_v4.csv

9. Result / results_examples/
   Selected project outputs used in the dissertation, including pLDDT, relaxation and repeated-inference analysis results.

10. COMP390_inference_demo_fixed.ipynb
   A notebook-based runnable demo for the reduced 2OLO inference example.

Running a minimal inference demo
A minimal demo can be run using the reduced 2OLO example together with the local OpenFold code and the selected checkpoint.

Key required components are:
- openfold/run_pretrained_openfold.py
- openfold/resources/params/params_model_1.npz
- openfold/resources/openfold_params/finetuning_2.pt
- checkpoints/run_subset_gpu32_h100_ep33/csv_logs/version_0/checkpoints/32-330000.ckpt
- demo_single_2olo/fasta_dir/
- demo_single_2olo/taskcif/
- demo_single_2olo/alignments_cf/

The software environment should be recreated using the files in environment/. For detailed commands, see:

docs/INFERENCE_DEMO_README.txt

Notes on scope
- This archive is a reduced dissertation artifact, not a full production-ready redistribution of OpenFold resources.
- Full alignment databases, full OpenFold training-resource downloads, raw datasets, temporary cache directories and all intermediate HPC files are not included.
- Some scripts and notebook cells retain project-specific or historical Barkla HPC paths for documentation purposes.
- These historical commands are included to show how the project was developed, but they are not all expected to run directly from the submitted archive.
- The reduced 2OLO inference demo is the intended runnable example.
- The dissertation should be treated as the main guide to the purpose and interpretation of each workflow component.

Suggested reading order
1. Dissertation PDF
2. docs/INFERENCE_DEMO_README.txt
3. COMP390_inference_demo_fixed.ipynb
4. docs/README.txt
5. environment/
6. task_scripts/
7. checkpoints/
8. demo_single_2olo/
9. Result / results_examples/

End of README