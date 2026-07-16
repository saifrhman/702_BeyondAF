1. **inv (Compute Invariants)**

    Calculates invariants from raw PDB files and saves them as CSVs, including Backbone Rigid Invariants (BRI), Length Angle Invariants (LAI) and Backbone Torsion Invariants (BTI).

    Syntax: 
    ```bash
    pdb_assess inv [INPUT_DIR] [OUTPUT_DIR] [OPTIONS]
    ```

    - Input Expected: A directory containing .pdb files.

    - Options:

        -n, --num-processor <int>: Defines the number of CPU cores to utilize for parallel processing. Defaults to half of your system's available cores.

2. compare (Distance Matrices)

    Measures structural distances between all processed protein chains in the input directory using Root Mean Square (RMS) and Chebyshev (L-infinity) metrics, outputting symmetric distance matrices and a comprehensive CSV table.

    Syntax: 
    ```bash
    pdb_assess compare [INPUT_DIR] [OUTPUT_DIR]
    ```

    - Input Expected: A directory containing computed invariant CSVs (typically the output from the inv command).

3. proj (Scatter Projections)

    Calculates statistical summaries (mean and standard deviation) for the invariants and generates 2D scatter plots to project and visualize these distributions across the dataset.

    Syntax: 
    ```bash
    pdb_assess proj [INPUT_DIR] [OUTPUT_DIR]
    ```

    - Input Expected: A directory containing computed invariant CSVs.

4. **plot (Structure Comparison)**

    Generates high-resolution overlay plots comparing the invariant sequences of your computed models against a known experimental reference structure.

    Syntax: 
    ```bash
    pdb_assess plot [INPUT_DIR] [OUTPUT_DIR] [OPTIONS]
    ```

    - Input Expected: A directory containing computed invariant CSVs.

    - Options:

        -s, --structure <string>: (Required) The identifier for the experimental structure to serve as the baseline comparison (e.g., 1HHO-1-A-1-141).

5. **pipe (Full Pipeline)**

    A unified command that sequentially executes inv, compare, and proj. Ideal for processing a raw batch of PDBs through to final statistical projections in a single execution.

    Syntax: 
    ```bash
    pdb_assess pipe [INPUT_DIR] [OUTPUT_DIR] [OPTIONS]
    ```

    - Input Expected: A directory containing raw .pdb files.

    - Options:

        -n, --num-processor <int>: Defines the number of CPU cores to utilize for the invariant computation step. Defaults to half of your system's available cores.