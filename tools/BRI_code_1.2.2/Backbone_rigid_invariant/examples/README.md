# Examples

## Introduction

There are 4 executable examples named as `{task_name}_example.py`, showing how this program cleans PDB data, computes Backbone Rigid Invariants with triangular invariants and Backbone Rigid Average Invariants, and quickly compare BRIs using kd-tree to find near-duplicates.

1. [cleaning_example.py](cleaning_example.py)

   Cleans PDB data with defined regulations and produces qualified protein chains.
   - input: [example_input/cleaning_example_input.txt](example_input/cleaning_example_input.txt), a list of PDB entry IDs.
   - output: 
     - [example_output/cleaned_connective_chains.csv](example_output/cleaned_connective_chains.csv), a csv file including identifier of cleaned chains.
     - [example_output/dropout_chains.csv](example_output/dropout_chains.csv), a csv file including information of chains being dropped out.
2. [invariant_compute_example.py](invariant_compute_example.py)

   Computes Backbone Rigid Invariants (BRI) with given chain identifiers.
   - input: [example_input/invariant_and_summary_compute_example_input.csv](example_input/invariant_and_summary_compute_example_input.csv), a csv file including identifier of cleaned chains.
   - output: [example_output/invariant_compute_output](example_output/invariant_compute_output), a folder containing BRI of chains stored by chain IDs.
4. [summary_compute_example.py](summary_compute_example.py)

   Computes Average Invariants (Brain) with given chain identifiers.
   - input: [example_input/invariant_and_summary_compute_example_input.csv](example_input/invariant_and_summary_compute_example_input.csv), a csv file including identifier of cleaned chains.
   - output: [example_output/summary_example_output.csv](example_output/summary_example_output.csv), a csv file that includes Brains of cleaned chains, one row per chain.
5. [duplicate_search_example.py](duplicate_search_example.py)

   Computes maximal distance of pairs of chains' BRIs and hamming distance of sequences.
   - input: [example_input/duplicate_search_example_input.csv](example_input/duplicate_search_example_input.csv), a csv file including identifier of cleaned chains.
   - output: [example_output/duplicate_search_example_results.csv](example_output/duplicate_search_example_results.csv), a csv file that includes the maximum distances for pairs of chain BRI and the Hamming distance for sequences.

## Usage

> To run a specific example, firstly please make sure the Python Interpreter is properly configured and dependencies in [requirements.txt](../requirements.txt) have been installed.

1. Navigate to the examples directory:

   ```bash
   cd ./examples
   ```

2. Run the example script:

   ```bash
   python {task_name}_example.py
   ```

