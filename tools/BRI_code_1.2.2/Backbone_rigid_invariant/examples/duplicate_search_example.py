import os
import sys

sys.path.append(
    os.path.join(os.path.dirname(__file__), "..")
)  # in case the package is not properly included
import pandas as pd
from bri import MiniChain
from bri.invariant_compare import invariant_nearest_neighbor_compare


if __name__ == "__main__":
    # read chain identifiers
    input_file = "example_input/duplicate_search_example_input.csv"
    output_folder = "example_output"
    cleaned_input = pd.read_csv(
        input_file, usecols=list(range(6)), na_values="", keep_default_na=False
    )
    cleaned_input = cleaned_input.drop(columns=["seq"])
    cleaned_input = cleaned_input.rename(columns={"start_residue": "start_index"})

    compare_results = []

    # extract invariants of the same length
    chain_groups = cleaned_input.groupby("chain_length")
    for length, chain_group in chain_groups:
        chain_ids = chain_group.to_dict("records")
        invariants_same_length = []
        for chain_id in chain_ids:
            # get invariants of chains
            invariant = MiniChain(**chain_id).invariant
            invariant["pdb_id"] = chain_id["pdb_id"]
            invariants_same_length.append(invariant)

        invariants_same_length = pd.concat(
            invariants_same_length, ignore_index=True
        )
        fmt_bri_group = invariants_same_length[
            [
                "pdb_id",
                "model_id",
                "chain_id",
                "residue_id",
                "residue_label",
                "x(N)",
                "y(N)",
                "z(N)",
                "x(C)",
                "y(C)",
                "z(C)",
                "x(A)",
                "y(A)",
                "z(A)",
                "chain_length",
            ]
        ]  # choose columns

        # compare invariants
        compare_results.append(
            invariant_nearest_neighbor_compare(fmt_bri_group, seq_compare=True)
        )

    compare_results = pd.concat(compare_results, ignore_index=True)
    compare_results.to_csv(
        f"{output_folder}/duplicate_search_example_results.csv", index=False
    )
    print(len(compare_results))
