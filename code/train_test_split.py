import numpy as np
import constants
import split_dataset as sd
import argparse

import sys, os
sys.path.append("./")

import pandas as pd
import split_dataset as sd

BASE_PATH = "../data"
SEEDS = [10, 24, 36, 42, 56]


# ----------------------------
# REDUCED TRAINING SET SIZES
# ----------------------------
TRAIN_SIZES = {
    "avgfp": [50, 500, 1000, 5000, 10000, 25000, 35000],
    "bgl3": [50, 500, 1000, 5000, 10000, 15000, 20000],
    "gb1": [100, 500, 1000, 5000, 10000, 25000, 50000, int(1e6)],
    "pab1": [50, 500, 1000, 5000, 10000, 25000, 30000],
    "ube4b": [50, 500, 1000, 5000, 10000, 25000, 50000, 70000],
}

def generate_reduced_splits(protein_name, seeds, ratios):
    """Generate splits with different reduced training set sizes."""
    train_ns = TRAIN_SIZES.get(protein_name)
    if train_ns is None:
        raise ValueError(f"Protein '{protein_name}' not supported for reduced splits. If this is a new protein, please add it to TRAIN_SIZES.")

    data_path = os.path.join(BASE_PATH, protein_name, f"{protein_name}.tsv")
    if not os.path.exists(data_path):
        raise FileNotFoundError(f"Dataset file not found: {data_path}")

    print(f"Generating reduced splits for {protein_name} with seeds {seeds}...")

    for seed in seeds:
        for train_n in train_ns:
            ds = pd.read_csv(data_path, sep="\t")
            split, split_dir = sd.train_tune_test(
                ds,
                train_size=ratios[0],
                tune_size=ratios[1],
                test_size=ratios[2],
                train_n=int(train_n),
                rseed=seed,
                out_dir=os.path.join(BASE_PATH, protein_name, "splits"),
                overwrite=True
            )
            print(f"train_n={train_n}, seed={seed}: {split_dir}")

    print("Reduced splits complete.")


def generate_standard_splits(protein_name, seeds, ratios):
    """Generate standard random and supertest splits."""
    print(f"Generating standard splits for {protein_name}...")

    ds_fn = constants.DATASETS[protein_name]["ds_fn"]
    ds = pd.read_csv(ds_fn, sep="\t")
    out_dir = os.path.join("data", protein_name, "splits")

    for seed in seeds:
        # Standard split
        split, split_dir = sd.train_tune_test(
            ds,
            train_size=ratios[0],
            tune_size=ratios[1],
            test_size=ratios[2],
            rseed=seed,
            out_dir=out_dir,
            overwrite=True
        )
        print(f"Random split (seed={seed}): {split_dir}")

def generate_supertest_split(protein_name, seeds):
    """Generate standard random and supertest splits."""
    print(f"Generating standard splits for {protein_name}...")

    ds_fn = constants.DATASETS[protein_name]["ds_fn"]
    ds = pd.read_csv(ds_fn, sep="\t")
    out_dir = os.path.join("data", protein_name, "splits")

    supertest_idxs, supertest_fn = sd.supertest(ds, size=0.1, rseed=seeds, out_dir=out_dir, overwrite=True)
    print(f"Supertest split saved: {supertest_fn}")


def generate_mutation_split(protein_name, seed):
    """Generate mutation-based splits."""
    print(f"Generating mutation-based splits for {protein_name}...")

    ds_fn = constants.DATASETS[protein_name]["ds_fn"]
    ds = pd.read_csv(ds_fn, sep="\t")
    out_dir = os.path.join("data", protein_name, "mutation_splits")

    mut_split, mut_split_dir, _ = sd.mutation_split(
        ds,
        train_muts_size=0.7,
        tune_size=0.5,
        out_dir=out_dir,
        rseed=seed,
        overwrite=True
    )
    print(f"Mutation split saved: {mut_split_dir}")

def generate_position_split(protein_name, seed):
    """Generate position-based splits."""
    print(f"Generating position-based splits for {protein_name}...")

    ds_fn = constants.DATASETS[protein_name]["ds_fn"]
    ds = pd.read_csv(ds_fn, sep="\t")
    out_dir = os.path.join("data", protein_name, "position_splits")

    seq_len = len(constants.DATASETS[protein_name]["wt_aa"])
    wt_ofs = constants.DATASETS[protein_name]["wt_ofs"]

    pos_split, pos_split_dir, _ = sd.position_split(
        ds,
        seq_len=seq_len,
        wt_ofs=wt_ofs,
        train_pos_size=0.7,
        tune_size=0.5,
        out_dir=out_dir,
        rseed=seed,
        overwrite=True
    )
    print(f"Position split saved: {pos_split_dir}")

def main():
    parser = argparse.ArgumentParser(
        description="Generate dataset splits for protein datasets."
    )

    parser.add_argument(
        "--protein", "-p", required=True, type=str,
        help="Protein name (e.g., avgfp, pab1, gb1)"
    )

    parser.add_argument(
        "--split_type", "-t", required=True, type=str, choices=["reduced", "standard", "supertest", "mutation", "position"],
        help="Type of dataset split to generate"
    )

    parser.add_argument(
        "--seeds", "-s", nargs="+", type=int, default=[42],
        help="List of random seeds to use for splits"
    )

    parser.add_argument(
        "--ratios", "-r", nargs=3, type=float, default=[0.8, 0.1, 0.1],
        metavar=("TRAIN", "TUNE", "TEST"),
        help="Train, tune, test split ratios (must sum to 1.0)"
    )

    args = parser.parse_args()

    protein = args.protein.lower()
    split_type = args.split_type.lower()
    seeds = args.seeds
    ratios = args.ratios

    if not np.isclose(sum(ratios), 1.0, atol=1e-3):
        parser.error(f"Ratios must sum to 1.0 (got {sum(ratios):.3f})")

    print(f"Starting dataset preparation for protein: {protein}")
    print(f"Split type: {split_type}")
    print(f"Seeds: {seeds}")
    print(f"Ratios (train/tune/test): {ratios}")

    if split_type == "reduced":
        generate_reduced_splits(protein, seeds, ratios)
    elif split_type == "standard":
        generate_standard_splits(protein, seeds, ratios)
    elif split_type == "supertest":
        generate_supertest_split(protein, seeds)
    elif split_type == "mutation":
        generate_mutation_split(protein, seeds)
    elif split_type == "position":
        generate_position_split(protein, seeds)
    else:
        parser.error(f"Unknown split type: {split_type}")

    print(f"Finished generating '{split_type}' splits for {protein}!")

if __name__ == "__main__":
    main()
