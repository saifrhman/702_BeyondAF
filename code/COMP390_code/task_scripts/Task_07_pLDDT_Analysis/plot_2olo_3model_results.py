from pathlib import Path
import argparse
import pandas as pd
import matplotlib.pyplot as plt


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--accuracy-dir", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    prediction_df = pd.read_csv(
        args.accuracy_dir / "2olo_accuracy_per_prediction.csv"
    )
    residue_df = pd.read_csv(
        args.accuracy_dir / "2olo_accuracy_per_residue.csv"
    )

    order = ["AlphaFold", "OpenFold", "Retrained_epoch_32"]

    summary = (
        prediction_df.groupby("model")
        .agg(
            mean_plddt=("mean_plddt", "mean"),
            std_plddt=("mean_plddt", "std"),
            mean_rmsd=("ca_rmsd_angstrom", "mean"),
            std_rmsd=("ca_rmsd_angstrom", "std"),
        )
        .reindex(order)
    )

    plt.figure(figsize=(7, 5))
    plt.bar(
        summary.index,
        summary["mean_plddt"],
        yerr=summary["std_plddt"],
        capsize=5,
    )
    plt.ylabel("Mean pLDDT")
    plt.title("2OLO prediction confidence")
    plt.xticks(rotation=15)
    plt.tight_layout()
    plt.savefig(args.out_dir / "2olo_mean_plddt_by_model.png", dpi=200)
    plt.close()

    plt.figure(figsize=(7, 5))
    plt.bar(
        summary.index,
        summary["mean_rmsd"],
        yerr=summary["std_rmsd"],
        capsize=5,
    )
    plt.ylabel("CA RMSD to experimental 2OLO (Å)")
    plt.title("2OLO structural accuracy")
    plt.xticks(rotation=15)
    plt.tight_layout()
    plt.savefig(args.out_dir / "2olo_ca_rmsd_by_model.png", dpi=200)
    plt.close()

    residue_summary = (
        residue_df.groupby(["model", "matched_pair_index"])
        .agg(mean_error=("ca_error_angstrom", "mean"))
        .reset_index()
    )

    plt.figure(figsize=(10, 5))
    for model in order:
        data = residue_summary[residue_summary["model"] == model]
        plt.plot(
            data["matched_pair_index"],
            data["mean_error"],
            label=model,
        )

    plt.xlabel("Aligned residue position")
    plt.ylabel("Mean CA error (Å)")
    plt.title("Residue-level structural error for 2OLO")
    plt.legend()
    plt.tight_layout()
    plt.savefig(
        args.out_dir / "2olo_residue_error_profile.png",
        dpi=200,
    )
    plt.close()

    plt.figure(figsize=(8, 5))
    for model in order:
        data = residue_df[residue_df["model"] == model]
        plt.scatter(
            data["plddt"],
            data["ca_error_angstrom"],
            s=8,
            alpha=0.25,
            label=model,
        )

    plt.xlabel("pLDDT")
    plt.ylabel("CA error after alignment (Å)")
    plt.title("pLDDT versus structural error")
    plt.legend()
    plt.tight_layout()
    plt.savefig(
        args.out_dir / "2olo_plddt_vs_structural_error.png",
        dpi=200,
    )
    plt.close()

    print(summary)
    print(f"\nSaved figures to: {args.out_dir}")


if __name__ == "__main__":
    main()
