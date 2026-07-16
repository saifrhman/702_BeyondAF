import multiprocessing as mp
import tqdm
from matplotlib import pyplot as plt
import pandas as pd

from bri.main import *
from bri import MiniEntry
from bri.filter import entry_integrated_cleaning, minientry_integrated_cleaning, integrated_chainwise_filter

# clean, dirty = entry_integrated_cleaning(r"..\tests\test_data\AF-A0A0G3QJE1-F1-model_v4_sel80.0_extract.cif")
e = MiniEntry(r"..\tests\test_data\AF-A0A0G3QJE1-F1-model_v4_sel80.0_extract.cif")
dirty = integrated_chainwise_filter(e.chains[0].get_feature())
# clean, dirty = entry_integrated_cleaning(r"101M")


def task(pid):
    results = []
    try:
        entry = MiniEntry(pid, extra_keys={'metric': 'ma_qa_metric_global'})
        for chain in entry.chains:
            results.append({'pdb_id': pid, 'chain_id': chain._chain_id, 'num_atoms': entry.coordinates.shape[0], 'pLDDt': entry.metric[0]['metric_value']})
        return pd.DataFrame(results)
    except Exception as e:
        print(f"Error processing {pid}: {e}")
        return pd.DataFrame()

if __name__ == '__main__':

    dataset = input_process('csm_clash_ids.txt', 0, None)

    # Task start
    start_time = time.time()
    pool = mp.Pool(mp.cpu_count())
    result = tqdm.tqdm(pool.imap_unordered(task, dataset), total=len(dataset))

    result = list(result)
    pool.close()

    # Post-process based on task type
    result = [df for df in result if not None]
    result_df = pd.concat(result)
    result_df.to_csv('atom_counts.csv', index=False)
