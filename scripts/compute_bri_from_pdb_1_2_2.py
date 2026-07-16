from pathlib import Path
import argparse

from bri import MiniChain


def main():
    parser = argparse.ArgumentParser(
        description="Compute BRI/LAI/BTI invariants from a single PDB file."
    )
    parser.add_argument("--pdb", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    if not args.pdb.is_file():
        raise FileNotFoundError(args.pdb)

    args.output.parent.mkdir(parents=True, exist_ok=True)

    chain = MiniChain.from_pdb(args.pdb)
    invariant = chain.get_chain_invariant(angles=True)

    if invariant is None or invariant.empty:
        raise RuntimeError(f"No invariants produced for {args.pdb}")

    required_columns = {
        "residue_id",
        "x(N)", "y(N)", "z(N)",
        "x(A)", "y(A)", "z(A)",
        "x(C)", "y(C)", "z(C)",
        "length(N)", "length(A)", "length(C)",
        "angle(N)", "angle(A)", "angle(C)",
        "tau(NA)", "tau(AC)", "tau(CN)",
    }

    missing = required_columns.difference(invariant.columns)
    if missing:
        raise RuntimeError(
            f"Missing expected invariant columns: {sorted(missing)}"
        )

    invariant.to_csv(args.output, index=False)

    print(f"Input: {args.pdb}")
    print(f"Residue rows: {len(invariant)}")
    print(f"Chains: {sorted(invariant['chain_id'].astype(str).unique())}")
    print(f"Saved: {args.output}")


if __name__ == "__main__":
    main()
