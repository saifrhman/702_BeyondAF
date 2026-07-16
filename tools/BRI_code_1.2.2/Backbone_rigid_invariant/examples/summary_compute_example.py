import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))  # in case the package is not properly included
import pandas as pd
from util.pdbx2df import PDB
from util.invariant import get_invariant_summary


if __name__ == '__main__':
    # read chain identifiers
    input_file = 'example_input/invariant_and_summary_compute_example_input.csv'
    output_folder = 'example_output'
    cleaned_input = pd.read_csv(input_file, usecols=list(range(6)), na_values='', keep_default_na=False)

    # map entry ID with chain identifiers
    entry_config_groups = cleaned_input.groupby('pdb_id')
    entry_config_dict = {entry_id: chains.to_dict('records') for entry_id, chains in entry_config_groups}

    # compute and save results
    entry_summary = []
    for entry_id, chains in entry_config_dict.items():
        entry_invariant = PDB(entry_id).get_entry_invariant(chains)  # compute invariant of an entry with cleaned chains
        output_split_group = entry_invariant.groupby(['model_id', 'chain_id'])  # split results by chain name

        for chain_name, chain_invariant in output_split_group:
            if len(chain_invariant) > 1:
                summary = get_invariant_summary(chain_invariant)    # invariant-summary
                summary['pdb_id'] = entry_id
                summary[['model_id', 'chain_id']] = chain_name
                summary = summary.iloc[:, list(range(-3, len(summary.columns) - 3))]
                entry_summary.append(summary)
    # save summary table
    pd.concat(entry_summary, ignore_index=True).to_csv('example_output/summary_example_output.csv', index=False)

