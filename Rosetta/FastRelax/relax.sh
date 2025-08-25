#!/usr/bin/bash
export ROSETTA=/path_to_rosetta/rosetta_src_2018.33.60351_bundle/main
export ROSETTA_BIN=$ROSETTA/source/bin
export ROSETTA_DB=$ROSETTA/database

# Do the mutation and get the energetics
echo "Mutation in progress..."
$ROSETTA_BIN/rosetta_scripts.linuxgccrelease @flags_mutate -out:level 200 

echo "Energy breakdown in progress..."
$ROSETTA_BIN/residue_energy_breakdown.linuxgccrelease @flags_score -out:level 200

# Additional information
# Read more here: https://docs.rosettacommons.org/docs/latest/full-options-list
# Can use different levels: 0 - fatal, 100 - error, 200 - warning, 300 - info, 400 - debug, 500 - trace.

