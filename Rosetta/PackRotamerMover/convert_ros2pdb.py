import pandas as pd
import re

# input the variant file obtained from ddG_calculation.py
data = pd.read_csv("variant.csv")
# The offset between the Rosetta and Original residue index
offset = 0
# The residue number from which the offset will start
offset_residue_num = 1
# New file name
new_file = "new_file.csv"

# Get all the variants from the dataframe
variant = data["variant"]

# Funtion to change the numbering from Rosetta index to the origianl PDB index
def add_to_number_in_variant(variant, offset, offset_residue_num):
    """
    Input Parameters:
    ----------------
    variants (str): currnet variant in the list
    offset (int): The offset between the Rosetta and Original residue index    
    offset_residue_num (int): The residue number from which the offset will start

    Output Parameters:
    -----------------
    new_variant (str): Original residue index in the current varaint
    """

    # It divides the variant into wild-type residue alphabet, Rosetta residue index, mutated residue alphabet
    match = re.match(r"([A-Z]+)(\d+)([A-Z]+)", variant)
    prefix, number, suffix = match.groups()
    new_number = int(number)
    # The offset between the Rosetta and Original residue index
    if int(number) >= offset_residue_num:
        new_number = int(number) + offset
    new_variant = f"{prefix}{new_number}{suffix}"
    return new_variant

new_variants = [add_to_number_in_variant(var, offset, offset_residue_num) for var in variant]
data["variant"] = new_variants
print(data["variant"]) 

# new datafram with correct residue index
data.to_csv(new_file, index=False)
