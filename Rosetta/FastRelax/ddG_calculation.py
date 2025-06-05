import argparse
import subprocess
import shutil
import os

from gen_rosetta_args import gen_rosetta_args
import multiprocessing


def copyComplete(source, target):
    # copy content, stat-info (mode too), timestamps...
    shutil.copy2(source, target)
    # copy owner and group
    st = os.stat(source)
    os.chown(target, st.st_uid, st.st_gid)


def cal_ddg(varaint_name, variant, args):
    """
    - Funtion where to calculate ddG and energy terms for each variant
    Input Parameters:
        - vid (str): name of the varaint
        - variant (str): variant to do the calculation
    """
    # Create the output directory for each variant
    directory = varaint_name
    files_dir = args.files_raw
    parent_dir = args.save_raw
    path = os.path.join(parent_dir, directory)
    os.mkdir(path)
    # Copies the PDB and other important file in the output directory
    struct =  str(parent_dir) + str(varaint_name) + "/structure.pdb"
    shutil.copyfile(args.pdb_fn, struct)
    flgm = str(parent_dir) + str(varaint_name) + "/flags_mutate"
    shutil.copyfile(f"{files_dir}flags_mutate", flgm)
    flgs = str(parent_dir) + str(varaint_name) + "/flags_score"
    shutil.copyfile(f"{files_dir}flags_score", flgs)
    lis = str(parent_dir) + str(varaint_name) + "/list.txt"
    shutil.copyfile(f"{files_dir}list.txt", lis)

    osd = str(parent_dir) + str(varaint_name)
    gen_rosetta_args(variant, osd, files_dir)
    shs = str(parent_dir) + str(varaint_name) + "/relax.sh"
    copyComplete("relax.sh", shs)

    os.chdir(str(parent_dir) + str(varaint_name))
    prd = str(parent_dir) + str(varaint_name) + "/relax.sh"
    # Do the process in parallel
    process = subprocess.Popen(prd,  shell=True)
    process.wait()

def main(args):
    # Open the variant file where the list of variant
    with open(args.variants_fn, "r") as file:
        ids_variants = file.readlines()
        #ids_variants = [line.strip() for line in file]
    # Number of varaints you want to at once
    x = args.cpu_core
    y = 0
    z = 0
    p_var = []
    # Run the selected number of varaints in parallel
    for k in range(0,len(ids_variants),x):
       z = z + x
       for j in range(y,z,1):
            vid, variant = ids_variants[j].split()
            p_var.append(f"p{j}")
            p_var[j] = multiprocessing.Process(target=cal_ddg, args=(vid, variant, args))
       for j in range(y,z,1):
            p_var[j].start()
       for j in range(y,z,1):
            p_var[j].join()
       y = z

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--variants_fn", help="the file containing the variants",type=str,default="variant.txt")
    parser.add_argument("--job_id", help="job id is used to save a diagnostic file", type=str, default="no_job_id")
    parser.add_argument("--pdb_fn", help="path to pdb file", type=str, default="Bgl3.pdb")
    parser.add_argument("--save_raw", help="set this to save the raw score.sc and energy.txt files in addition to the parsed ones", type=str, default="./rosetta_working_dir/")
    parser.add_argument("--files_raw", help="where needed files for the calculations are saved", type=str, default="./rosetta_working_dir/")
    parser.add_argument("--cpu_core", help="Number of cpu want to run", type=int, default=10)
    main(parser.parse_args())
