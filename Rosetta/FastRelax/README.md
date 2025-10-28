# Rosetta Variant List Creation and  FastRelax Procedure

Tutourial describes how to generate a variant list using Rosetta, relax a protein structure with FastRelax, and run the `ddG_calculation.py` script on the lowest energy pose.

## Steps
###  1. Relax the Structure (FastRelax)
Relax the structure preferably 5 times to obtain energetically favorable conformations:
``` Run in command line
$ROSETTA_BIN/relax.default.linuxgccrelease \
    -s test.pdb \
    -restore_talaris_behavior \
    -nstruct 1 \
    -relax:default_repeats 5 \
    -out:path:pdb .
```

### 2. Generate the Variant File
Generate a variant list for a given structure using the following PyRosetta script `variant_file.py`:
``` 
python variant_file.py --pdb_file test.pdb --output_file variant.txt --target_chain A --offset 0
```

### 3. Calculate Rosetta ddG score
Before running, ensure that the paths in the following files are correctly set:
- `relax.py`
- `templates/mutate_template.xml`
- `templates/mutation_template.resfile`

Then, calculate ddG scores using:
```
python ddG_calculation.py \
    --pdb_fn test.pdb \
    --variants_fn variant.txt \
    --save_raw /absolute/path/to/save \
    --raw_files /path/to/rosetta_working_dir/
```
**Notes**: Use the lowest energy relaxed structure for ddG calculations.

### 4. Combine All Scores into a CSV File
- To generate a CSV file containing the raw scores from all variants, run:
    ```
    ./Combine_outputs.sh FolderName OutputCSVFile
    ```
- To generate a CSV file where the scores are subtracted relative to the wild-type, run:
    ```
    python combine_substract.py FolderName OutputCSVFile
    ```

## Acknowledgment
The code is based on the [nn4dms](https://github.com/gitter-lab/nn4dms).

