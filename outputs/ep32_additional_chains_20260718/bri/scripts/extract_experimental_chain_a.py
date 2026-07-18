#!/usr/bin/env python3

from pathlib import Path
import argparse

from Bio.PDB import MMCIFParser, PDBIO, Select
from Bio.PDB.Polypeptide import is_aa


class ChainASelect(Select):
    def accept_model(self, model):
        return model.id == 0

    def accept_chain(self, chain):
        return chain.id == "A"

    def accept_residue(self, residue):
        return residue.id[0] == " " and is_aa(residue, standard=False)

    def accept_atom(self, atom):
        altloc = atom.get_altloc()
        return altloc in (" ", "A", "1")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cif", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    cif_path = Path(args.cif)
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    parser_obj = MMCIFParser(QUIET=True, auth_chains=True)
    structure = parser_obj.get_structure(cif_path.stem, str(cif_path))

    model = next(structure.get_models())
    chain_ids = [chain.id for chain in model]

    print(f"Available chains in {cif_path.name}: {chain_ids}")

    if "A" not in chain_ids:
        raise RuntimeError(
            f"Chain A was not found in {cif_path}. Available chains: {chain_ids}"
        )

    io = PDBIO()
    io.set_structure(structure)
    io.save(str(output_path), ChainASelect())

    atom_count = 0
    with output_path.open() as handle:
        for line in handle:
            if line.startswith("ATOM"):
                atom_count += 1

    if atom_count == 0:
        raise RuntimeError(f"No ATOM records were written to {output_path}")

    print(f"Wrote {output_path}")
    print(f"ATOM records: {atom_count}")


if __name__ == "__main__":
    main()
