COMP390 Inference Demo Guide
Project: AlphaFold and Other Machine Learning Methods for Proteins

Purpose
This guide explains how to run a minimal inference example from the extracted
COMP390_code folder.

The demo uses a reduced single-target example for 2OLO and compares three inference
settings that differ only in the source of the model parameters or checkpoint:

1. AlphaFold parameter run
   This run uses a DeepMind AlphaFold parameter file in JAX/NPZ format:
   openfold/resources/params/params_model_1.npz.
   OpenFold officially supports inference with AlphaFold's original pretrained
   parameters through the --jax_param_path argument.

2. OpenFold official checkpoint run
   This run uses an official OpenFold-trained parameter checkpoint in PyTorch format:
   openfold/resources/openfold_params/finetuning_2.pt.
   OpenFold documentation describes this as the alternative public monomer-inference
   option to DeepMind's AlphaFold parameters, loaded through
   --openfold_checkpoint_path.

3. Project-trained checkpoint run
   This run uses the project checkpoint retained in the submitted artifact:
   checkpoints/run_subset_gpu32_h100_ep33/csv_logs/version_0/checkpoints/32-330000.ckpt.
   It is loaded through the same OpenFold checkpoint interface as the official
   OpenFold parameter run, but the weights come from this project's resumed
   OpenFold training workflow rather than from the official OpenFold release.

These three runs therefore provide a compact comparison between:
- the original AlphaFold-style pretrained parameter baseline;
- the official OpenFold-trained checkpoint baseline; and
- the checkpoint produced by this dissertation project.

This demo is the runnable component of the submitted archive. The full project workflow
also involved data preparation, checkpoint-resume training, repeated inference,
relaxation analysis, pLDDT extraction, BRI/LAI/BTI analysis and visualisation. However,
the full OpenFold training resources, complete alignment databases, raw downloaded
datasets, temporary cache directories and all intermediate HPC outputs are not included
because they are too large for submission. The broader documentation notebook should be
read as a project workflow record, not as a fully runnable end-to-end reproduction script.

Assumed starting point
After extracting the archive, enter the project root:

cd COMP390_code

All paths below are written relative to this folder.

Required files
The following files and folders must be present inside COMP390_code:

- openfold/run_pretrained_openfold.py
- openfold/resources/params/params_model_1.npz
- openfold/resources/openfold_params/finetuning_2.pt
- checkpoints/run_subset_gpu32_h100_ep33/csv_logs/version_0/checkpoints/32-330000.ckpt
- demo_single_2olo/fasta_dir/2olo.fasta
- demo_single_2olo/taskcif/2olo.cif
- demo_single_2olo/alignments_cf/2olo_A/

Environment setup
This project was run in a conda environment called openfold_env.

A similar environment can be recreated using:

conda env create -f environment/openfold_env_full.yml
conda activate openfold_env

If the full export does not solve cleanly, try:

conda env create -f environment/openfold_env_from_history.yml
conda activate openfold_env
pip install -r environment/openfold_env_pip_freeze.txt

If the environment already exists, activate it with:

conda activate openfold_env

Before running inference
From the COMP390_code directory, set:

export PYTHONPATH=$(pwd)/openfold:$PYTHONPATH

Optional: set a local Triton cache directory to avoid using an NFS-backed default cache:

export TRITON_CACHE_DIR=$(pwd)/triton_cache
mkdir -p "$TRITON_CACHE_DIR"

Check that a GPU is visible:

nvidia-smi

Create output folders:

mkdir -p demo_run_outputs_2olo/alphafold
mkdir -p demo_run_outputs_2olo/openfold_official
mkdir -p demo_run_outputs_2olo/your_model_ep33

Optional quick file check
Before running inference, the following commands may be used to confirm that the main inputs exist:

ls openfold/run_pretrained_openfold.py
ls openfold/resources/params/params_model_1.npz
ls openfold/resources/openfold_params/finetuning_2.pt
ls checkpoints/run_subset_gpu32_h100_ep33/csv_logs/version_0/checkpoints/32-330000.ckpt
ls demo_single_2olo/fasta_dir/2olo.fasta
ls demo_single_2olo/taskcif/2olo.cif
ls demo_single_2olo/alignments_cf/2olo_A

Inference commands

1. AlphaFold parameter run

Purpose:
This command runs OpenFold inference while loading DeepMind AlphaFold model_1
parameters from the retained NPZ file. It serves as the AlphaFold-parameter baseline
in the reduced demo.

python openfold/run_pretrained_openfold.py \
  demo_single_2olo/fasta_dir \
  demo_single_2olo/taskcif \
  --output_dir demo_run_outputs_2olo/alphafold \
  --use_precomputed_alignments demo_single_2olo/alignments_cf \
  --config_preset model_1 \
  --jax_param_path openfold/resources/params/params_model_1.npz \
  --model_device cuda:0 \
  --data_random_seed 101

2. OpenFold official checkpoint run

Purpose:
This command runs inference with the official OpenFold finetuning_2 checkpoint
retained in the archive. It provides the OpenFold-trained public checkpoint baseline
for comparison with the AlphaFold-parameter run and the project checkpoint run.

python openfold/run_pretrained_openfold.py \
  demo_single_2olo/fasta_dir \
  demo_single_2olo/taskcif \
  --output_dir demo_run_outputs_2olo/openfold_official \
  --use_precomputed_alignments demo_single_2olo/alignments_cf \
  --config_preset model_1 \
  --openfold_checkpoint_path openfold/resources/openfold_params/finetuning_2.pt \
  --model_device cuda:0 \
  --data_random_seed 101

3. Project-trained checkpoint run

Purpose:
This command runs inference with the project-trained checkpoint selected from the
resume-training workflow used in the dissertation. It uses the same inference script,
input FASTA file, precomputed alignments and model configuration as the other two
commands, but the checkpoint source is this project's epoch-33 retained checkpoint.

python openfold/run_pretrained_openfold.py \
  demo_single_2olo/fasta_dir \
  demo_single_2olo/taskcif \
  --output_dir demo_run_outputs_2olo/your_model_ep33 \
  --use_precomputed_alignments demo_single_2olo/alignments_cf \
  --config_preset model_1 \
  --openfold_checkpoint_path checkpoints/run_subset_gpu32_h100_ep33/csv_logs/version_0/checkpoints/32-330000.ckpt \
  --model_device cuda:0 \
  --data_random_seed 101

Interpretation of the three runs
All three commands use:
- the same 2OLO FASTA input;
- the same task CIF directory;
- the same precomputed alignment directory;
- the same model configuration preset, model_1;
- the same CUDA device selection and random seed.

The main intended difference is therefore the parameter source:
- AlphaFold run: DeepMind AlphaFold parameters loaded from an NPZ/JAX parameter file;
- OpenFold official run: OpenFold's released trained checkpoint loaded from a PyTorch PT file;
- Project checkpoint run: the dissertation project's selected training checkpoint loaded
  from a PyTorch CKPT file.

This design makes the demo a small controlled comparison of parameter/checkpoint sources
within the same OpenFold inference pipeline.

Expected outputs
The three commands should create:

- demo_run_outputs_2olo/alphafold
- demo_run_outputs_2olo/openfold_official
- demo_run_outputs_2olo/your_model_ep33

Each output directory should contain a predictions/ subdirectory. The expected prediction
files include:

- 2olo_A_model_1_unrelaxed.pdb
- 2olo_A_model_1_relaxed.pdb

The exact auxiliary output files may depend on the OpenFold script behaviour.

Smoke-test note
This reduced 2OLO demo was tested on the slimmed COMP390_code archive. The AlphaFold
parameter run, OpenFold official checkpoint run and project-trained checkpoint run all
completed inference and relaxation successfully, producing both unrelaxed and relaxed
PDB outputs.

GPU note
The commands use:

--model_device cuda:0

So a CUDA-capable GPU is expected. If a different GPU index is needed, change cuda:0
accordingly.

Recommended execution order
1. Extract COMP390_code.zip
2. cd COMP390_code
3. Recreate and activate the conda environment
4. Export PYTHONPATH
5. Optionally set TRITON_CACHE_DIR
6. Check that a GPU is available with nvidia-smi
7. Create the output directories
8. Run one or more of the three inference commands above
9. Check the output predictions/ directory

End of guide