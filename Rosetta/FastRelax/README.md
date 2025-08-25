# PyRosetta Variant List Creation and  FastRelax Procedure

This document outlines the steps to generate a variant list using PyRosetta, relax the structure, and run the `ddG_calculation.py` script with the lowest energy pose.

## Steps

### 1. Create Variant List
To create a variant list for a given structure, use the following PyRosetta script:

```bash
python variant_file.py --pdb_file test.pdb --output_file varaint.txt --target_chain A --offset 0
```

###  2. Relax the Structure N (preferably 5) times

``` Run in command line
$ROSETTA_BIN/relax.default.linuxgccrelease -s test.pdb -restore_talaris_behavior -nstruct 2 -relax:default_repeats 5 -out:path:pdb .
```

### 3. Create variant file using variant_file.py
``` 
python variant_file.py --pdb_file test.pdb --output_file varaint.txt --target_chain A --offset 0
```

### 4. Calculate the Rosetta ddG score
```
Before that please change the directories in the the path of Rosetta dir in relax.py, templates/mutate_template.xml and templates/mutation_template.resfile files.

python ddG_calculation.py --pdb_fn structure.pdb --variants_fn variant.txt --save_raw dir_to_save --files_raw rosetta_working_dir/

structure.pdb: lowest energy structure
variant.txt: list of varaints for the ddG calculation
dir_to_save: absolute path of the directory where you want to save your data
--files_raw: rosetta_working_dir/
```
