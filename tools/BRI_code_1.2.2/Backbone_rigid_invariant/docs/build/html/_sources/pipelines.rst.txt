Use of Pipelines
================

This package so far provides convenient functions to streamline common 
workflows cleaning protein chains, computing Backbone Rigid Invariants, and 
searching near-duplicates in groups.

This section introduces the available pipelines and demonstrates how they can 
be used directly from the command line. 

Overview
--------

The pipelines are designed to:

* Clean protein chain data from mmCIF files 
* Compute Backbone Rigid Invariants (BRIs) in a reproducible manner
* Search for near-duplicate structures within or across datasets

All pipelines are exposed through a command-line interface (CLI), making them 
easy to integrate into shell scripts and larger workflows.

General Command-Line Usage
--------------------------

After installing the package, the main entry point can be accessed via ::

    python -m bri --task <pipeline> [OPTIONS]

Where:

* ``<pipeline>`` specifies the workflow to execute
* ``[OPTIONS]`` control input files, output locations, and algorithm-specific 
  parameters

You can list all available pipelines with ::

    python -m bri --help

.. And obtain help for a specific pipeline with:

.. ```
.. package-name <pipeline> --help
.. ```

Pipeline 1: Cleaning Protein Chains
-----------------------------------

This pipeline applies standard cleaning steps, including selection of chains, 
and removal of chains broken by missing atoms and residues, and atom clashes.

Example
^^^^^^^
::

    python -m bri --task clean \
        --input-data data_list/entry_ids.txt \
        --output results/ \

Options
^^^^^^^

* ``--input-data``: Path to a directory or ``.txt`` file containing accessible protein 
  structure data
* ``--output``: Directory to save cleaned protein chains
* ``--max-samples``: Number of entries to use from the input data. Default is 
  to use all.
* ``--process``: If to show a process bar. Default is not to show.
* ``--extend``: If to use an extended workflow. Defaults to `True`. Use 
  ``--no-extend`` to disable.

    .. note::
       For the cleaning pipeline, the extended version includes additional 
       steps including entity inspection and author data extraction. The 
       output of cleaned chains will be represented with columns:

        * pdb_id: represents the entry name of the mmCIF file
        * entity_id: represents the labelled entity id of the chain
        * model_id: represents the model id (pdbx_PDB_model_num) of the chain
        * chain_id: represents the chain id (label_asym_id)
        * start_residue: represents the first residue number (label_seq_id) of 
          the cleaned chain
        * chain_length: represents the length of this cleaned chain
        * auth_chain_id: represents the author-assigned chain id (auth_asym_id)
        * auth_seq_id_start: represents the first author-assigned residue 
          number (auth_seq_id) of the cleaned chain  
        * auth_seq_id_end: represents the last author-assigned residue number 
          (auth_seq_id) of the cleaned chain  
        * seq: represents the amino acid sequence of the cleaned chain
    
       However, entity information was not always included in protein strcutures 
       (especially computed strcutures). In this case, the non-extended version 
       is recommended to ensure compatibility. Also, the output of it will be 
       more concise as the minimal required information for identifying a 
       protein chain:

        * pdb_id
        * model_id
        * chain_id
        * start_residue
        * chain_length
        * seq

Outputs
^^^^^^^

* File of cleaned protein chains ``chains_cleaned.csv`` written to the 
  specified output directory
* File of removed items ``chains_filtered.csv`` written to the specified 
  output directory

Pipeline 2: Computing Backbone Rigid Invariants
-----------------------------------------------

This pipeline computes Backbone Rigid Invariants (BRIs) for cleaned protein 
chains. 
.. These invariants provide a compact geometric description of backbone conformations.

Example
^^^^^^^
::

    python -m bri --task bri \
        --input results/chains_cleaned.csv \
        --output results/bri/ \

Options
^^^^^^^

* ``--input-data``: A ``.csv`` file containing cleaned protein chains, a ``.txt`` file 
  containing accessable entries, or a directory of mmCIF/bcif files.
* ``--output``: Directory to save BRIs seperately for each chain, or a single file for 
  all chains.
* ``--max-samples``: Number of entries to use from the input data. Default is to use all.
* ``--process``: If to show a process bar. Default is not to show.
* ``--extend``: If to use an extended workflow. Defaults to `True`. Use ``--no-extend`` 
  to disable.

    .. note:: 
       For the BRI computation pipeline, the extended version will include bond angles 
       and tausion angles additional to the standard BRIs.

Outputs
^^^^^^^

* A file containing all BRIs, or files saved in the target directory seperately by chains.

Pipeline 3: Searching for Near-Duplicates
-----------------------------------------

This pipeline searches for near-duplicate protein structures based on their BRIs.

Example 
^^^^^^^
::

    python -m bri --task duplicate \
        --input results/chains_cleaned.csv \
        --output results/duplicates

Options
^^^^^^^

* ``--input``: A ``.csv`` file containing cleaned protein chains.
  `More input options will be added in the future.`
* ``--output``: If not specified, the full results of pairwise distances will be saved 
  into the directory by chain length. Two tables will always be saved even without 
  the specified output option:

    1. ``L_inf_lt1_<time>.csv``: All pairs of structures with L-infinity distances 
       less than 1.0 Angstrom.
    2. ``L_inf_lt001_eq_seq_<time>.csv``: All pairs of structures with the same amino 
       sequences and L-infinity distances less than 0.01 Angstrom.
* ``--max-samples``: Number of entries to use from the input data. Default is to use all.

Outputs
^^^^^^^

* ``.csv`` files listing pairs of near-duplicate structures

Chaining Pipelines
------------------

The pipelines are designed to be composable. For example, cleaning, invariant 
computation, and duplicate search can be executed sequentially::

    python -m bri --task clean --input entry_ids.txt --output results/
    python -m bri --task bri --input results/chains_cleaned.csv --output results/bri/
    python -m bri --task duplicate --input results/chains_cleaned.csv --output results/

Notes
-----

* Logs will be printed and saved in the current working directory by default for 
  troubleshooting.

.. Further details on each pipeline and their underlying algorithms are provided in the following sections.
