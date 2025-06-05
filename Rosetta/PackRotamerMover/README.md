# Structure Relaxation and ddG Calculation Pipeline

This pipeline performs structure relaxation and calculates the ddG for all possible variants of each position. Follow the steps below to run the necessary scripts.

## Steps

1. **Relax the Structure:**
   - Use the `relax.py` script to relax the initial structure.
   - Provide the initial PDB file as input.
   - Choose an optimized `max_iter` value for the relaxation process.
   - Run `relax.py` with the optimized `max_iter` 4 times to obtain the best-optimized structure.

2. **Feed the Lowest Energy Pose:**
   - After performing the relaxation steps, select the lowest energy pose.
   - Feed this lowest energy pose into the `run.sh` script.

3. **Run ddG Calculations:**
   - The `run.sh` script will automatically call `ddG_calculation.py` for all possible variants at each position.
   - Provide the lowest ebery input pdb and the total residue number

## Notes

- **Residue Index Mismatch:**
  - The residue index in the original PDB file may differ from the index used by Rosetta.
  - If you encounter this issue, use the `convert_ros2pdb.py` script to convert between the Rosetta and PDB residue numbering.

## Example Usage

```bash (give 0 max iteration if want to find the max_iteration)
python relax.py initial_structure.pdb 0
./run.sh

