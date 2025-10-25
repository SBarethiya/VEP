import os, sys
import pandas as pd
import re

# Directory containing variant folders
DIR = sys.argv[1]
output_file = sys.argv[2]

# Columns to include
COLUMNS = [
    "total_score","dslf_fa13","fa_atr","fa_dun","fa_elec","fa_intra_rep",
    "fa_rep","fa_sol","hbond_bb_sc","hbond_lr_bb","hbond_sc","hbond_sr_bb",
    "omega","p_aa_pp","pro_close","rama","ref","yhh_planarity"
]

# Step 1: Read first numeric SCORE line from each folder
scores_dict = {}   # key: folder, value: list of floats

for folder in os.listdir(DIR):
    folder_path = os.path.join(DIR, folder)
    if not os.path.isdir(folder_path):
        continue

    score_file = os.path.join(folder_path, "score.sc")
    if not os.path.isfile(score_file):
        continue

    with open(score_file, "r") as f:
        lines = f.readlines()

    # Find first numeric SCORE line
    for line in lines:
        if line.startswith("SCORE:") and "total_score" not in line:
            # Remove SCORE: and split into floats
            values = line.strip().split()[1:19]  # first 18 numeric values
            values = [float(v) for v in values]
            scores_dict[folder] = values
            break

# Step 2: Build reference map: key = first letter + number, value = folder of wild-type
ref_map = {}
for folder in scores_dict.keys():
    match = re.match(r"([A-Z])(\d+)([A-Z])", folder)
    if match:
        wt, res_num, mut = match.groups()
        key = f"{wt}{res_num}"        # first letter + residue number
        # If mutation = wild-type, this folder is the reference
        if mut == wt:
            ref_map[key] = folder

# Step 3: Subtract reference scores
output_rows = []
for folder, values in scores_dict.items():
    match = re.match(r"([A-Z])(\d+)([A-Z])", folder)
    if not match:
        continue
    wt, res_num, mut = match.groups()
    key = f"{wt}{res_num}"
    if key not in ref_map:
        # No reference found, skip or use zeros
        print(f"No reference found for {folder}, skipping.")
        continue

    ref_folder = ref_map[key]
    ref_values = scores_dict[ref_folder]

    # Subtract reference column-wise
    subtracted = [round(v - rv, 6) for v, rv in zip(values, ref_values)]
    subtracted.append(folder)  # add mutation column
    output_rows.append(subtracted)

# Step 4: Write to CSV
df = pd.DataFrame(output_rows, columns=COLUMNS + ["mutation"])
df.to_csv(output_file, index=False)
print(f"CSV with subtracted values created: {output_file}")

