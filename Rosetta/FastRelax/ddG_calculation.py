import argparse
import subprocess
import shutil
import os
import time
import multiprocessing

from gen_rosetta_args import gen_rosetta_args

def copyComplete(source, target):
    # copy content, stat-info (mode too), timestamps...
    shutil.copy2(source, target)
    # copy owner and group
    st = os.stat(source)
    os.chown(target, st.st_uid, st.st_gid)


def calculate_DDG(varaint_index, variant, chain, offset, args):
    """
    - Funtion where to calculate ddG and energy terms for each variant
    Input Parameters:
        - varaint_index (str): name of the varaint
        - variant (str): variant to do the calculation
        - offset (int): offset in the residue number.
        - chain (str): do the mutation on which chain.
    """

    # Create the output directory for each variant
    path = os.path.join(args.save_raw, varaint_index)
    os.mkdir(path)

    # Copies the PDB and other important file in the output directory
    structure_file =  str(args.save_raw) + str(varaint_index) + "/structure.pdb"
    shutil.copyfile(args.pdb_fn, structure_file)
    flgm_file = str(args.save_raw) + str(varaint_index) + "/flags_mutate"
    shutil.copyfile(f"{args.raw_files}flags_mutate", flgm_file)
    flgs_file = str(args.save_raw) + str(varaint_index) + "/flags_score"
    shutil.copyfile(f"{args.raw_files}flags_score", flgs_file)
    list_file = str(args.save_raw) + str(varaint_index) + "/list.txt"
    shutil.copyfile(f"{args.raw_files}list.txt", list_file)

    # generate the rosetta arguments for this variant
    variant_dir = str(args.save_raw) + str(varaint_index)
    gen_rosetta_args(variant, variant_dir, args.raw_files, offset, chain)

    # copy relax file and run relax script
    relax_save = str(args.save_raw) + str(varaint_index) + "/relax.sh"
    copyComplete("relax.sh", relax_save)
    os.chdir(str(args.save_raw) + str(varaint_index))
    relax_cmd = "./relax.sh" # str(args.save_raw) + str(varaint_index) + "/relax.sh"
    os.chmod("./relax.sh", 0o755) 

    # Do the process in parallel
    process = subprocess.Popen(relax_cmd,  shell=True)
    process.wait()

def main(args):
    """
    - Main function to process a list of variants in parallel using multiprocessing.
    """

    # Read the list of variants from the specified file
    with open(args.variants_fn, "r") as file:
        ids_variants = file.readlines()

    # Initialize indices for batching variants based on available CPU cores
    start_idx = 0
    end_idx = 0
    processes = []

    # Process variants in batches parallely
    for batch_start in range(0, len(ids_variants), args.cpu_core):
       end_idx += args.cpu_core
       start_time = time.time()
       for j in range(start_idx, end_idx):
            if j >= (len(ids_variants)):
                continue

            idx, variant, chain, offset = ids_variants[j].split()
            print(f"Processing: {idx}, {variant} on chain {chain} with offset {offset}")

            process = multiprocessing.Process(target=calculate_DDG, args=(idx, variant, chain, offset, args))
            processes.append(process)
            
       # Start all processes in the batch
       for i in range(start_idx, end_idx):
           if i < len(processes):
               processes[i].start()

       # Wait for all processes in the batch to complete
       for i in range(start_idx, end_idx):
           if i < len(processes):
               processes[i].join()

       end_time = time.time()
       print(f"Batch completed in {end_time - start_time:.2f} seconds")
       start_idx = end_idx

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("--variants_fn", help="the file containing the variants",type=str, default="./variant.txt")

    parser.add_argument("--job_id", help="job id is used to save a diagnostic file", type=str, default="no_job_id")

    parser.add_argument("--pdb_fn", help="path to pdb file", type=str, default="./test.pdb")

    parser.add_argument("--save_raw", help="set this to save the raw score.sc and energy.txt files in addition to the parsed ones", type=str, default="./rosetta_working_dir/")

    parser.add_argument("--raw_files", help="where are the needed files for the calculations are saved", type=str, default="./rosetta_working_dir/")

    parser.add_argument("--cpu_core", help="Number of cpu want to run", type=int, default=3)

    main(parser.parse_args())
