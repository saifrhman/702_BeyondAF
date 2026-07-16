# -*- coding = utf-8 -*-
# @File: invariant_compute.PY
import os
import sys

sys.path.append(
    os.path.join(os.path.dirname(__file__), "..")
)  # in case the package is not properly included
import pandas as pd
from bri import Chain


if __name__ == "__main__":
    # read chain identifiers
    input_file = "example_input/invariant_and_summary_compute_example_input.csv"
    output_folder = "example_output/invariant_compute_output"
    cleaned_input = pd.read_csv(
        input_file, usecols=list(range(6)), na_values="", keep_default_na=False
    )

    # map entry ID with chain identifiers
    entry_config_dict = cleaned_input.to_dict("records")

    # compute and save results
    for chains in entry_config_dict:
        chain = Chain(**chains)
        chain_invariant = (
            chain.get_chain_invariant()
        )  # compute invariant of cleaned chains

        # save by their chain names
        chain_invariant.to_csv(f"{output_folder}/{chain}.csv", index=False)
