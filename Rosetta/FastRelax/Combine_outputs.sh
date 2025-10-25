#!/bin/bash

# Output CSV file
#OUTPUT="scores_with_mutation.csv"
#Folder="test"

Folder="$1"
OUTPUT="$2"

# Define the column names you want
COLS="total_score,dslf_fa13,fa_atr,fa_dun,fa_elec,fa_intra_rep,fa_rep,fa_sol,hbond_bb_sc,hbond_lr_bb,hbond_sc,hbond_sr_bb,omega,p_aa_pp,pro_close,rama,ref,yhh_planarity,mutation"

# Write header to CSV
echo "$COLS" > "$OUTPUT"

# Loop over all subfolders in test/
for folder in $Folder/*/; do
    mutation=$(basename "$folder")   # folder name
    score_file="$folder/score.sc"

    # Skip if file does not exist
    [ -f "$score_file" ] || continue

    # Extract numeric SCORE line (skip header)
    line=$(grep "^SCORE:" "$score_file" | tail -n 1 | sed 's/^SCORE:[[:space:]]*//')

    # Convert spaces to commas
    line=$(echo "$line" | tr -s ' ' ',')

    # Take only the first 18 values
    line=$(echo "$line" | cut -d',' -f1-18)

    # Append folder name as mutation
    echo "$line,$mutation" >> "$OUTPUT"
done

echo "CSV file generated: $OUTPUT"

