# -*- coding = utf-8 -*-
# @File: cleaning.PY
import os
import sys

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))  # in case the package is not properly included

import pandas as pd
import numpy as np
from bri.filter import integrated_chainwise_filter, field_check
from bri.base.base_util import on_entry
from bri.pdbx2df import Entry


def entry_integrated_cleaning(pdb_id):
    entry = field_check(pdb_id)
    if not isinstance(entry, Entry):
        return None, pd.DataFrame(entry, index=[0])

    entry_integrated_filter = on_entry()(integrated_chainwise_filter)
    cleaning_res = entry_integrated_filter(entry)
    if cleaning_res is None:
        return None

    dirty_set = cleaning_res[cleaning_res['type'] != 'clean']
    if 'chain_length' in dirty_set.columns:
        dirty_set = dirty_set[['pdb_id', 'model_id', 'chain_id', 'residue_id', 'type', 'residue_label', 'chain_length']]
    else:
        dirty_set = dirty_set[['pdb_id', 'model_id', 'chain_id', 'residue_id', 'residue_label', 'type']]
        dirty_set['chain_length'] = np.nan

    clean_set = cleaning_res[cleaning_res['type'] == 'clean']
    clean_set = clean_set.groupby(['pdb_id', 'model_id', 'chain_id']).agg(
        start_residue=('residue_id', 'min'), seq_max=('residue_id', 'max'),
        residue_label=('residue_label', lambda x: ''.join(x[::3]))).reset_index()
    clean_set['chain_length'] = clean_set.apply(lambda row: row.seq_max - row.start_residue + 1, axis=1)
    clean_set.start_residue = clean_set.start_residue.astype('int')
    clean_set.chain_length = clean_set.chain_length.astype('int')

    chain_entity_dict = {i['id']: int(i['entity_id']) for i in entry.label_entity_info}
    clean_set['entity_id'] = clean_set['chain_id'].map(chain_entity_dict)
    clean_set = clean_set.drop(columns=['seq_max'])
    clean_set = clean_set.rename(columns={'residue_label': 'seq'})
    return clean_set, dirty_set



if __name__ == '__main__':
    with open('example_input/cleaning_example_input.txt', 'r', encoding='utf-8') as f:
        entry_list = f.read().split(',')

    # obtain online protein data from RCSB by entry IDs, and start cleaning
    cleaning_result = [entry_integrated_cleaning(id) for id in entry_list]
    # remove None type
    cleaning_result = [i for i in cleaning_result if i is not None]
    # separate cleaned chains and noisy chains
    clean_chains = [result[0] for result in cleaning_result if result[0] is not None]
    noisy_chains = [result[1] for result in cleaning_result if result[1] is not None]
    # concatenate all results into one DataFrame
    clean_chains = pd.concat(clean_chains, ignore_index=True)
    noisy_chains = pd.concat(noisy_chains, ignore_index=True)
    # save results
    clean_chains = clean_chains[['pdb_id', 'entity_id', 'model_id', 'chain_id', 'start_residue', 'chain_length', 'seq']]
    clean_chains.to_csv('example_output/cleaned_connective_chains.csv', index=False)
    noisy_chains.to_csv('example_output/dropout_chains.csv', index=False)

    print('number of cleaned chains:', len(clean_chains))
    print('number of noisy chains:', len(noisy_chains))
