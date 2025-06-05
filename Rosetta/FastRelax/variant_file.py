from pyrosetta import *
import numpy as np
import sys
init()

'''
This file creates the every individual the variant file at every individual residue position.
'''

amino_acids = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]

def main():
    pdb_file = sys.argv[1]
    file_name = sys.argv[2]
    pose = pose_from_pdb(pdb_file)
    pose_sequence_len = len(pose.sequence())
    variants = [(f"{i}_{aa}", f"{pose.residue(i).name1()}{i}{aa}") for i in range(1, pose_sequence_len+1, 1) for aa in amino_acids]
    np.savetxt(file_name, np.array(variants), fmt="%s")
     
if __name__ == "__main__":
    main()

