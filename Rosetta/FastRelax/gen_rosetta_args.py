import argparse
import os

def gen_rosetta_script_xml_str(variant, template_fn, offset, chain):
    """
    Generate RosettaScript XML string from mutation variant(s).    
    Parameters:
        variant (str): Mutation string, e.g., "A23D,C45K"
        template_fn (str): Path to the template directory (without filename)
        offset (int): offset in the residue number.
        chain (str): do the mutation on which chain.
    Returns:
        str: Filled RosettaScript XML string
    """

    # ROSETTA USES 1-BASED INDEXING
    resnums = []
    for mutation in variant.split(","):
        resnum = int(mutation[1:-1])  # middle is the residue number
        rosetta_idx_offset = resnum + offset  # add offset
        resnums.append(f"{rosetta_idx_offset}{chain}")  # e.g., 24A

    resnum_str = ",".join(resnums)

    # load the XML template
    with open(f"{template_fn}templates/mutate_template.xml", "r") as f:
        template_str = f.read()

    # fill in the template
    formatted_template = template_str.format(resnum_str)
    return formatted_template


def gen_resfile_str(variant, template_fn, offset, chain):
    """
    Generate a Rosetta resfile string from mutation data: residue_number chain PIKAA replacement_AA.
    Parameters:
        variant (str): Mutation string like "A23D,C45K"
        template_fn (str): Path to the folder containing mutation_template.resfile
        offset (int): offset in the residue number.
        chain (str): do the mutation on which chain.
    Returns:
        str: Filled resfile content
    """

    mutation_strs = []

    for mutation in variant.split(","):
        resnum = int(mutation[1:-1])
        resnum_idx_offset = resnum + offset
        new_aa = mutation[-1]

        mutation_strs.append("{} {} PIKAA {}".format(resnum_idx_offset, chain, new_aa))

    # add new lines between mutation strs
    mutation_strs = "\n".join(mutation_strs)

    # load the template
    with open(f"{template_fn}templates/mutation_template.resfile", "r") as f:
        template_str = f.read()

    formatted_template = template_str.format(mutation_strs)
    return formatted_template


def gen_rosetta_args(variant, output_dir, rosetta_files, offset, chain):
    """
    Generate and write RosettaScript XML and resfile based on variant.

    Parameters:
        variant (str): Mutation string like "A23D, C45K"
        output_dir (str): Directory to write output files
        rosetta_files (str): Path prefix to template files (no trailing slash)
        offset (int): offset in the residue number.
        chain (str): do the mutation on which chain.
    Returns:
        dict: Paths to generated XML and resfile
    """

    rosetta_script_str = gen_rosetta_script_xml_str(variant, rosetta_files, int(offset), chain)
    resfile_str = gen_resfile_str(variant, rosetta_files, int(offset), chain)

    with open(os.path.join(output_dir, "mutate.xml"), "w") as f:
        f.write(rosetta_script_str)
    
    with open(os.path.join(output_dir, "mutation.resfile"), "w") as f:
        f.write(resfile_str)


def main(args):
    gen_rosetta_args(args.variant, args.output_dir, args.rosetta_files, args.offset, args.chain_id)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate Rosetta resfile and XML for a mutation",
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("--variant",
                        help="the variant for which to generate rosetta args",
                        type=str)

    parser.add_argument("--output_dir",
                        help="the directory in which to save mutation.resfile and mutate.xml",
                        type=str)

    parser.add_argument("--rosetta_files",
                        help="the directory in which mutation.resfile and mutate.xml exist",
                        type=str)

    parser.add_argument("--offset",
                        help="If there is any offset in the residue number",
                        type=int)

    parser.add_argument("--chain_id",
                        help="Chain id on which mutation will be done",
                        type=str)
    main(parser.parse_args())
