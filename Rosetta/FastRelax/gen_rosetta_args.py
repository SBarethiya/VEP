""" Generates arguments for Rosetta run """

import argparse
import os


def gen_rosetta_script_str(variant, template_fn):

    # for the rosetta script xml, we just need the mutated residue numbers and chain
    # Note: ROSETTA USES 1-BASED INDEXING
    resnums = []
    for mutation in variant.split(","):
        resnum_0_idx_raw = int(mutation[1:-1])
        resnum_0_idx_offset = resnum_0_idx_raw 
        resnum_1_index = resnum_0_idx_offset + 1
        resnums.append("{}A".format(resnum_1_index))
        # resnums.append("{}".format(resnum_1_index))

    resnum_str = ",".join(resnums)

    # load the template
    with open(f"{template_fn}templates/mutate_template.xml", "r") as f:
        template_str = f.read()

    # fill in the template
    formatted_template = template_str.format(resnum_str)

    return formatted_template


def gen_resfile_str(variant, template_fn):
    "residue_number chain PIKAA replacement_AA"

    mutation_strs = []
    for mutation in variant.split(","):
        resnum_0_idx_raw = int(mutation[1:-1])
        resnum_0_idx_offset = resnum_0_idx_raw
        resnum_1_index = resnum_0_idx_offset + 1
        new_aa = mutation[-1]

        mutation_strs.append("{} A PIKAA {}".format(resnum_1_index, new_aa))

    # add new lines between mutation strs
    mutation_strs = "\n".join(mutation_strs)

    # load the template
    with open(f"{template_fn}templates/mutation_template.resfile", "r") as f:
        template_str = f.read()

    formatted_template = template_str.format(mutation_strs)

    return formatted_template


def gen_rosetta_args(variant, out_dir, files_rose):

    rosetta_script_str = gen_rosetta_script_str(variant, files_rose)
    resfile_str = gen_resfile_str(variant, files_rose)

    with open(os.path.join(out_dir, "mutate.xml"), "w") as f:
        f.write(rosetta_script_str)
    
    with open(os.path.join(out_dir, "mutation.resfile"), "w") as f:
        f.write(resfile_str)


def main(args):
    gen_rosetta_args(args.vid,args.variant, 0, args.out_dir)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("--vid", 
                        help="the variant for which to generate rosetta args_gb1", 
                        type=str)

    parser.add_argument("--variant",
                        help="the variant for which to generate rosetta args_gb1",
                        type=str)

    parser.add_argument("--out_dir",
                        help="the directory in which to save mutation.resfile and mutate.xml",
                        type=str)

    parser.add_argument("--files_rose",
                        help="the directory in which to save mutation.resfile and mutate.xml",
                        type=str)
    main(parser.parse_args())
