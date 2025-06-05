from pyrosetta import *
init()  # Initializes the PyRosetta library

from rosetta.protocols import *
sfxn = get_fa_scorefxn()  # Create a score function object using the default "fa_scorefxn" (full atom scoring function).

pdb = sys.argv[1]
max_iter = sys.argv[2]
max_iter = int(max_iter)
print(max_iter)
lowest_energy_pdb = pdb
pose = pose_from_pdb(pdb)
lowest_energy = sfxn.score(pose) 

if max_iter == 0:
    # To find the maximum number of iteration
    for i in [100,200,300,400,500]:
        pose = pose_from_pdb(pdb)
        # Initialize the FastRelax protocol object, which performs energy minimization to relax the protein structure.
        fr = rosetta.protocols.relax.FastRelax()
        # Set the scoring function for the relaxation protocol.
        fr.set_scorefxn(sfxn)
        # Set the maximum number of iterations for the relaxation process.
        fr.max_iter(i)       
        # Apply the FastRelax protocol to the pose (structure), modifying the pose in place.
        fr.apply(pose)
        # Save the relaxed pose as a new PDB file with a unique name based on the iteration number.
        pose.dump_pdb(f'output_{i}.pdb')
        if sfxn.score(pose) <= lowest_energy:
            max_iter = i
            lowest_energy_pdb = f'output_{i}.pdb'
            lowest_energy = sfxn.score(pose)

# After finding the max_iter, Loop to perform relaxation 4 times
if max_iter != 0:
    for i in range(0, 4):
        # Load the protein structure from the PDB file for each iteration.
        pose = pose_from_pdb(pdb)
        # Initialize the FastRelax protocol object, which performs energy minimization to relax the protein structure.
        fr = rosetta.protocols.relax.FastRelax()
        # Set the scoring function for the relaxation protocol.
        fr.set_scorefxn(sfxn)
        # Set the maximum number of iterations for the relaxation process.
        fr.max_iter(max_iter)
        # Apply the FastRelax protocol to the pose (structure), modifying the pose in place.
        fr.apply(pose)
        # Save the relaxed pose as a new PDB file with a unique name based on the iteration number.
        pose.dump_pdb(f'output_{max_iter}_{i}.pdb')
        pose.dump_pdb(f'output_{max_iter}_{i}.pdb')
        if sfxn.score(pose) <= lowest_energy:
            lowest_energy_pdb = f'output_{max_iter}_{i}.pdb'
            lowest_energy = sfxn.score(pose)

print(f"The lowest energy pose is: {lowest_energy_pdb} with energy: {lowest_energy} and max iteration: {max_iter}")

