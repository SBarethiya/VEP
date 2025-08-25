from pyrosetta import *
from pyrosetta.toolbox import *
from pyrosetta.teaching import *
import pandas as pd
import csv
import os, sys
import numpy as np
#import  aspose.cells
from pyrosetta.rosetta import *
#from aspose.cells import Workbook
import re
from joblib import Parallel, delayed
import multiprocessing as mp
import time
init()

import predict_ddG


# Ref2015 score function
sfxn = get_fa_scorefxn()

# 20 canonical amino acids
amino_acids = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]

# Function for prediction of ddG, energy components and labels of components
def get_energy_components(native_pose, mutated_pose, sfxn):
    """
    Calculate the energy components and ddG of a protein structure before and after a mutation.

    Input Parameters:
    -----------------
    native_pose : Pose
        The native (unmutated) protein structure as a Pose object.
    mutated_pose : Pose
        The mutated protein structure as a Pose object.
    sfxn : ScoreFunction
        A ScoreFunction object used to evaluate and score the energy of the protein poses.

    Output Parameters:
    ------------------
    labels : List[str]
        A list of strings representing the names of the energy terms (e.g., "fa_atr", "fa_rep", "fa_sol").
    ddGs : List[float]
        A list of floating-point numbers representing the energ component values.
       These values are the difference in energy between the mutated and native poses for each term.
    score : float
        The overall ddG score (total energy difference) for the mutated pose versus the native pose.
    """  
    # Get the total energies for the native pose and mutated pose,
    # applying the score function (sfxn) weights to the energies.
    tmp_native = native_pose.energies().total_energies().weighted_string_of(sfxn.weights())
    tmp_mutant = mutated_pose.energies().total_energies().weighted_string_of(sfxn.weights())

    # Split the weighted energy string into a list and remove any None values.
    array_native = list(filter(None, tmp_native.split(' ')))
    array_mutant = list(filter(None, tmp_mutant.split(' ')))

    # Extract only the energy values for native pose (skip labels and take every other value in the array).
    native_scores = []
    for i in range(len(array_native)):
        if (i % 2 != 0):  # Every second element is an energy score.
            native_scores.append(float(array_native[i]))

    # Similarly, extract energy values for the mutated pose.
    mutant_scores = []
    for i in range(len(array_mutant)):
        if (i % 2 != 0):  # Every second element is an energy score.
            mutant_scores.append(float(array_mutant[i]))

    # Calculate the each energy component by subtracting
    # native pose scores from the corresponding mutant scores.
    ddGs = []
    for i in range(len(mutant_scores)):
        ddG_component = mutant_scores[i] - native_scores[i]  # ddG is the difference in energy components.
        ddGs.append(round(ddG_component, 3))  # Round each ddG component to 3 decimal places.

    # Extract the labels (energy terms) from the native pose energy string (every other element is a label).
    labels = []
    for i in range(len(array_native)):
        if (i % 2 == 0):  # Every other element is a label.
            labels.append(array_native[i].translate(':').strip(":"))  # Clean up the label text.

    # Calculate the total score difference (ddG) for the entire structure.
    native_score = sfxn.score(native_pose)  # Get the total score for the native pose.
    mutant_score = sfxn.score(mutated_pose)  # Get the total score for the mutated pose.
    score = mutant_score - native_score  # ddG for the entire structure is the difference in total scores.
    return labels, ddGs, score

def com_split(org_aa, resid_cur, mut_aa, pose):
    """
    Input Parameters:
    -----------------
    org_aa : str
        The original amino acid at the residue position (before mutation).
    resid_cur : int
        The residue index where the mutation is to occur.
    mut_aa : str
        The mutated amino acid that will replace the original one.
    pose : Pose
        The original protein pose (structure) before mutation.

    Output Parameters:
    ------------------
    pose_ref : Pose
        A new Pose object representing the structure with the original amino acid at the specified position.
    pose_mut : Pose
        A new Pose object representing the structure with the mutated amino acid at the specified position.
    """
    # Create a copy of the input pose to retain the original pose for reference
    pose_ref = Pose()
    pose_ref.assign(pose)
    
    # Create another copy of the input pose to apply the mutation
    pose_mut = Pose()
    pose_mut.assign(pose)
    
    # Mutate the residue at position 'resid_cur' to the original amino acid 'org_aa' (this serves as the reference pose)
    mutate_residue(pose_ref, resid_cur, org_aa, 12.0, sfxn)
    
    # Mutate the residue at position 'resid_cur' to the mutated amino acid 'mut_aa' (this serves as the mutated pose)
    mutate_residue(pose_mut, resid_cur, mut_aa, 12.0, sfxn)
    
    # Return the two poses: one with the original residue and one with the mutated residue
    return pose_ref, pose_mut

def worker(i, aa, pose, sfxn, queue):
    """
    Input Parameters:
    -----------------
    i : int
        The index of the residue to mutate in the protein structure.
    aa : str
        The amino acid to mutate the residue to.
    pose : Pose
        The protein pose to perform the mutation on.
    sfxn : ScoreFunction
        The ScoreFunction object used to evaluate and score the energy of the poses.
    queue : Queue
        A multiprocessing Queue to store the results from the worker.

    Output Parameters:
    ------------------
    result : tuple
        A tuple containing the calculated (ddGs), the individual component difference (score_s).
    """
    # Retrieve the wild-type at the specified residue index 'i'
    wt = pose.residue(i).name1()
    
    # Record the start time for performance tracking
    start = time.time()
    
    # Call the 'com_split' function to generate two poses: reference and mutated
    pose_ref, pose_mut = com_split(wt, i, aa, pose)
    
    # Get the energy components and ddG values by comparing the reference pose and mutated pose
    labels, dGs, score_s = get_energy_components(pose_ref, pose_mut, sfxn)
    
    # Create a result tuple containing ddG values, overall score difference, and mutation label
    result = (dGs, score_s, wt + str(i) + aa)
    
    # Calculate and print the runtime for this worker function
    run_time = time.time() - start
    print(f"{i} {aa} {dGs} {score_s} {sfxn(pose_mut)}")
    print("run_time", i, run_time)
    
    # Put the result into the queue for later collection
    queue.put(result)

def main():
    """
    The function reads mutation range and target PDB file, then performs mutations for each residue in the specified range, calculates energy components, and stores the results in a CSV file.
    """
    # Read the start and stop indices for the residue mutation range from bash file arguments
    start = sys.argv[1]
    stop = sys.argv[2]
    file_name = sys.argv[3]
    print(start, stop)  # Print start and stop values to the console for debugging
    start = int(start)  # Convert the start index to an integer
    stop = int(stop)    # Convert the stop index to an integer
    
    # Initialize lists to store results
    dG = []  # List to store ddG values (score difference between mutant and native)
    dG_comp_All = []  # List to store individual energy components
    variant_done = []  # List to store the identifiers of the variants processed

    # Load the PDB file and create a Pose object
    pose = pose_from_pdb(file_name)

    # Load the Ref2015 scoring function
    sfxn = get_fa_scorefxn()

    # Calculate the length of the protein sequence for the specified residue range
    pose_sequence_len = len(pose.sequence()[start:stop])

    # Generate all possible variants (residue index and amino acid) within the specified range
    variants = [(i, aa) for i in range(start, stop+1, 1) for aa in amino_acids]

    # Assume pose_info and sfxn_info are placeholders for data required to reconstruct the pose and scoring function
    pose_info = "pose_info_data"  # Placeholder for pose information
    sfxn_info = "sfxn_info_data"  # Placeholder for scoring function information

    # Initialize a queue to collect results from the worker processes
    queue = mp.Queue()
    processes = []  # List to keep track of all the processes

    # Loop through the variants and start a new process for each mutation
    for i, aa in variants:
        # Create a new process that will run the worker function for the current variant
        p = mp.Process(target=worker, args=(i, aa, pose, sfxn, queue))
        processes.append(p)  # Add the process to the list
        p.start()  # Start the process

    # Wait for all processes to finish
    for p in processes:
        p.join()

    # After all worker processes are done, retrieve results from the queue
    while not queue.empty():
        results = queue.get()  # Get the result from the queue
        dGs, score_s, variant = results  # Unpack the results into individual variables
        dG.append(score_s)  # Append the overall ddG value (score) to the list
        dG_comp_All.append(dGs)  # Append the individual ddG components to the list
        variant_done.append(variant)  # Append the variant identifier to the list

    # Prepare the results for saving in CSV format
    scores = []  # List to store the final result data to be saved

    # Loop through the dG list and create the final output format
    for i in range(0, len(dG), 1):
        v1 = dG_comp_All[i]  # Individual ddG components for the i-th variant
        v2 = []  # Temporary list to store data for the i-th variant
        for j in range(0, len(dG_comp_All[1]), 1):
            v2.append(v1[j])  # Add each ddG component to the list
        v2.append(dG[i])  # Add the overall ddG value (score) for this variant
        v2.append(variant_done[i]) 
        scores.append(v2)  # Append the result for this variant to the scores list

    # Save the results as a CSV file without headers and index, appending to an existing file
    pd.DataFrame(scores).to_csv('ddG_scores.csv', mode='a', header=False, index=False)
 
if __name__ == "__main__":
    main()

