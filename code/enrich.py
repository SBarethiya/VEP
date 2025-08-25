""" Common functions for creating enrich2 datasets (for GB1 and Bgl3) """
import os
from os.path import join
from subprocess import call
import pandas as pd
import utils


def create_e2_config_file(inp_fn, sel_fn, e2_output_dir, config_file_save_dir):
    """ create a config file specifying parameters for enrich2 """
    # text = """{"libraries": [
    #                {
    #                  "counts file": "INP_COUNTS_FN", 
    #                  "identifiers": {
    #                   "type": "protein",
    #                   "min count": 100, 
    #                   "use aligner": false, 
    #                   "wild type": {
    #                     "coding": "protein", 
    #                     "sequence": "MQYKLILNGKTLKGETTTEAVDAATAEKVFKQYANADNGVDGEWTYDDAATKTFTVTE"
    #                   }
    #                 }, 
    #                  "name": "T0", 
    #                  "report filtered reads": false, 
    #                  "timepoint": 0
    #                }, 
    #                {
    #                  "counts file": "SEL_COUNTS_FN", 
    #                  "identifiers": {
    #                   "type": "protein",
    #                   "min count": 100, 
    #                   "use aligner": false, 
    #                   "wild type": {
    #                     "coding": "protein", 
    #                     "sequence": "MQYKLILNGKTLKGETTTEAVDAATAEKVFKQYANADNGVDGEWTYDDAATKTFTVTE"
    #                   }
    #                 }, 
    #                  "name": "T1", 
    #                  "report filtered reads": false, 
    #                  "timepoint": 1
    #                }
    #              ], 
    #              "name": "C1", 
    #              "output directory": "OUTPUT_DIR"}"""

    text = """{"libraries": [
                   {
                     "counts file": "INP_COUNTS_FN", 
                     "identifiers": { }, 
                     "name": "T0", 
                     "report filtered reads": false, 
                     "timepoint": 0
                   }, 
                   {
                     "counts file": "SEL_COUNTS_FN", 
                     "identifiers": { }, 
                     "name": "T1", 
                     "report filtered reads": false, 
                     "timepoint": 1
                   }
                 ], 
                 "name": "C1", 
                 "output directory": "OUTPUT_DIR"}"""

    text = text.replace("INP_COUNTS_FN", inp_fn)
    text = text.replace("SEL_COUNTS_FN", sel_fn)
    text = text.replace("OUTPUT_DIR", e2_output_dir)

    with open(join(config_file_save_dir, "e2_config"), "w") as f:
        f.write(text)

    return join(config_file_save_dir, "e2_config")


def create_tsv_dataset(library_size, e2_scores_fn, save_fn):
    """ create a simple tsv dataset file using the output from enrich2 """
    e2_data = pd.HDFStore(e2_scores_fn)
    # get the e2 scores, removing the wild-type and moving the variant index into a column
    e2_scores = e2_data.select("/main/identifiers/scores")
    e2_scores = e2_scores.loc[e2_scores.index != "_wt"]
    e2_scores.reset_index(inplace=True)
    # get the input and selected counts
    e2_counts = e2_data.select("/main/identifiers/counts")
    if "_wt" in e2_counts.index:
        e2_counts = e2_counts.drop("_wt")
    variants = e2_scores["index"].values
    num_mutations = e2_scores["index"].apply(lambda x: len(x.split(","))).values
    scores = e2_scores["score"].values
    inp = e2_counts["c_0"].values
    sel = e2_counts["c_1"].values

    cols = ["variant", "num_mutations", "inp", "sel", "score"]
    data = {"variant": variants, "num_mutations": num_mutations, "inp": inp, "sel": sel, "score": scores}

    df = pd.DataFrame(data, columns=cols)
    df = df[:library_size]
    if save_fn is not None:
        df.to_csv(save_fn, sep="\t", index=False)

    return df

def ensure_dir_exists(path):
    if not os.path.exists(path):
        os.makedirs(path)

def create_e2_dataset(count_df, output_dir, library_size, output_fn=None):
    """ creates an enrich2 dataset (saves input files for enrich2, runs enrich2, then converts to my format) """
    # create a special output directory for enrich2
    e2_output_dir = join(output_dir, "e2_output")
    ensure_dir_exists(e2_output_dir)
    count_df.loc["_wt"] = {
        "variant": "_wt",
        "inp": 1759616,
        "sel": 3041819
    }
    # create enrich2 input files
    inp = count_df[["variant", "inp"]]
    sel = count_df[["variant", "sel"]]
    # inp.index = inp.index.astype(str)
    # sel.index = sel.index.astype(str)
    inp_fn = join(e2_output_dir, "idents_inp.tsv")
    sel_fn = join(e2_output_dir, "idents_sel.tsv")
    
    inp.to_csv(inp_fn, sep="\t", header=["index", "count"], index=False)
    sel.to_csv(sel_fn, sep="\t", header=["index","count"], index=False)

    # create enrich2 config file for this ds
    e2_config_fn = create_e2_config_file(inp_fn, sel_fn, e2_output_dir, e2_output_dir)

    # run e2
    call(['enrich_cmd', '--no-plots', '--no-tsv', e2_config_fn, 'ratios', 'wt'])
    if output_fn is None:
        output_fn = "dataset.tsv"

    create_tsv_dataset(library_size, join(e2_output_dir, "C1_sel.h5"), save_fn=join(output_dir, output_fn))
