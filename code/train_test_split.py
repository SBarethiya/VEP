import sys
sys.path.append("/home/shrishti/Documents/Projects/SequenceToStructure_Gitter/nn4dms-final/code")

import pandas as pd
import split_dataset as sd
import utils

protein_name = sys.argv[1]

base_path = "/home/shrishti/Documents/Projects/SequenceToStructure_Gitter/nn4dms-final"

if protein_name == "avgfp":
    train_ns = [50, 500, 1000, 5000, 10000, 25000, 35000 ]
elif protein_name == "bgl3":
    train_ns = [50, 500, 1000, 5000, 10000, 15000, 20000 ]
elif protein_name == "gb1":
    train_ns = [100, 500, 1000, 5000, 10000, 25000, 50000, 10e5 ]
elif protein_name == "pab1":
    train_ns = [50, 500, 1000, 5000, 10000, 25000, 30000 ]
elif protein_name == "ube4b":
    train_ns = [50, 500, 1000, 5000, 10000, 25000, 50000, 70000 ]

# load the full dataset (really, we only need the # of variants in the dataset, but this is easier)
for seed in [10, 24, 36, 42, 56]:
    for train_n in train_ns:
        ds = pd.read_csv(f"{base_path}/data/{protein_name}/{protein_name}.tsv", sep="\t")
        split, split_dir = sd.train_tune_test(ds, train_size=.81, tune_size=.09, test_size=.1, train_n=int(train_n), rseed=seed, out_dir=f"{base_path}/data/{protein_name}/splits", overwrite=True)
ds = pd.read_csv(f"{base_path}/data/{protein_name}/{protein_name}.tsv", sep="\t")
split, split_dir = sd.train_tune_test(ds, train_size=.81, tune_size=.09, test_size=.1, train_n=None, rseed=seed, out_dir=f"{base_path}/data/{protein_name}/splits", overwrite=True)
    
print(split_dir)
