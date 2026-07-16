from bri import Chain, Entry, MiniChain, MiniEntry
from matplotlib import pyplot as plt
from pathlib import Path


def run():
    from bri import MiniChain
    dir_path = Path('../../holomon_backbone_15us_dt1ns')
    files = list(dir_path.glob('*.pdb'))
    for f in files[:50]:
        a = MiniChain.from_pdb(f)
        l = a.invariant.residue_id.to_list()
        s = set(range(1, 314))
        print(s - set(l))


import atomium
atomium.fetch('9a9q')