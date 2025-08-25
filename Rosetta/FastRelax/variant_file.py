import argparse
from pyrosetta import *
import numpy as np
import sys
init()

'''
This file creates the every individual the variant file at every individual residue position.
'''

amino_acids = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]

def main():
    pdb_file = args.pdb_file
    output_file = args.output_file
    target_chain = args.target_chain
    offset = args.offset

    pose = pose_from_pdb(pdb_file)
    pdb_info = pose.pdb_info()
    variants = []

    for i in range(1, pose.total_residue() + 1):
        chain = pdb_info.chain(i)
        if chain != target_chain:
            continue

        orig_aa = pose.residue(i).name1()
        pdb_resnum = pdb_info.number(i)
        
        for new_aa in amino_acids:
            if new_aa == orig_aa:
                continue  # Skip self-mutations

            # e.g., "24_D" (custom ID), "A24D" (mutation string)
            id_str = f"{u}"
            mut_str = f"{orig_aa}{pdb_resnum}{new_aa}"
            chain_str = chain
            variants.append((id_str, mut_str, chain_str, offset))

    np.savetxt(output_file, np.array(variants), fmt="%s")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate variant file",
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("--pdb_file",
                        help="PDB file in which mutation will be done",
                        type=str)

    parser.add_argument("--output_file",
                        help="Output file name to save the variants",
                        type=str,
                        default="./varaint.txt")

    parser.add_argument("--target_chain",
                        help="Chain in which mutation will be done",
                        type=str)

    parser.add_argument("--offset",
                        help="If there is any offset in the residue number",
                        type=int,
                        default=0)

    main(parser.parse_args())
