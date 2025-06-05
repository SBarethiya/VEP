#!/usr/bin/bash

export ROSETTA=/home/shrishti/Documents/Softwares/Rosetta/rosetta_src_2018.33.60351_bundle/main
export ROSETTA_BIN=$ROSETTA/source/bin
export ROSETTA_DB=$ROSETTA/database
# this script might be called from the root directory, but it does all its work in the rosetta working directory
echo "do mutation in ROSETTA..."
$ROSETTA_BIN/rosetta_scripts.linuxgccrelease @flags_mutate -out:level 200
echo "do energy breakdown..."
$ROSETTA_BIN/residue_energy_breakdown.linuxgccrelease @flags_score -out:level 200
