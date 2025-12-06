""" Encodes data in different formats """

from os.path import join, isfile
import argparse

import numpy as np
import pandas as pd
from sklearn.preprocessing import OneHotEncoder, MinMaxScaler

import constants

def enc_aa_index(int_seqs):
    """
    encodes data in aa index properties format
    param: 
        int_seqs: integer encoded sequences
    return: 
        aa index encoded sequences 
    """
    aa_features = np.load("data/aaindex/pca-19_norm.npy")
    # add all zero features for stop codon
    aa_features = np.insert(aa_features, 0, np.zeros(aa_features.shape[1]), axis=0)
    aa_features_enc = aa_features[int_seqs]
    return aa_features_enc

def enc_rmsf(rmsf_file, int_seqs):
    """ 
    encodes data in RMSF format
    param: 
        rmsf_file: path to RMSF file
        int_seqs: integer encoded sequences
    return:
        rmsf encoded sequences
    """
    rmsf = pd.read_csv(rmsf_file,sep="\t")
    rmsf = np.array(rmsf["rmsf"])
    # add all zero features for stop codon
    rmsf = np.insert(rmsf, 0, np.zeros(rmsf.shape[0]), axis=0)
    rmsf_enc = rmsf[int_seqs]
    rmsf_enc = np.expand_dims(rmsf_enc, axis=-1)    
    return rmsf_enc

def enc_ss(ss_file, int_seqs):
    """
    encodes data in secondary structure format
    param: 
        ss_file: path to secondary structure file
        int_seqs: integer encoded sequences
    return:
        secondary structure encoded sequences
    """
    ss = np.load(ss_file)
    # add all zero features for stop codon
    ss = np.insert(ss, 0, np.zeros(ss.shape[0]), axis=0)
    ss_enc = ss[int_seqs]
    ss_enc = np.expand_dims(ss_enc, axis=-1)    
    return ss_enc

def enc_rosetta(rosetta_file, variant, encoded,wt_offset):
    """
    encodes data in rosetta score format
    param: 
        rosetta_file: path to rosetta score file
        variant: list of variants
        encoded: previously encoded sequences
        wt_offset: wild type offset
    return:
        rosetta score encoded sequences
    """
    rost = pd.read_csv(rosetta_file, sep="\t")
    scores = rost.loc[:, rost.columns != 'mutation']
    a = np.zeros_like(encoded[:,:,0])
    c = np.zeros_like(scores.iloc[0:1])
    b = c[:, np.newaxis]
    new_score = a[:, :, np.newaxis] + b
    for j, v in enumerate(variant):
        for m in v.split(","):
            for i, r in enumerate(rost['mutation']):
                if m == r:
                    position = int(m[1:-1])-wt_offset
                    new_score[j,position,:] = scores.iloc[i].values
    new_encoded = np.concatenate((encoded,new_score),axis=2)
    return new_encoded

def enc_one_hot(int_seqs):
    """
    encodes data in one-hot format
    param: 
        int_seqs: integer encoded sequences
    return:
        one-hot encoded sequences
    """
    enc = OneHotEncoder(categories=[range(constants.NUM_CHARS)] * int_seqs.shape[1], dtype=np.bool, sparse=False)
    one_hot = enc.fit_transform(int_seqs).reshape((int_seqs.shape[0], int_seqs.shape[1], constants.NUM_CHARS))
    return one_hot


def enc_int_seqs_from_char_seqs(char_seqs):
    """
    converts character sequences to integer sequences
    param: 
        char_seqs: list of character sequences
    return:
        seq_ints: numpy array of integer encoded sequences
    """
    seq_ints = []
    for char_seq in char_seqs:
        int_seq = [constants.C2I_MAPPING[c] for c in char_seq]
        seq_ints.append(int_seq)
    seq_ints = np.array(seq_ints)
    return seq_ints


def enc_int_seqs_from_variants(variants, wild_type_seq, wt_offset=0):
    """
    converts variants to integer sequences based on the wild-type sequence
    param: 
        variants: list of variants
        wild_type_seq: wild-type character sequence
        wt_offset: wild-type offset
    return:
        seq_ints: numpy array of integer encoded sequences
        wt_mut_seq: numpy array of integer encoded wild-type minus mutant sequences
    """
    wild_type_int = np.zeros(len(wild_type_seq), dtype=np.uint8)
    wild_type_int_wt = np.zeros(len(wild_type_seq), dtype=np.uint8)
    for i, c in enumerate(wild_type_seq):
        wild_type_int[i] = constants.C2I_MAPPING[c]

    seq_ints = np.tile(wild_type_int, (len(variants), 1))
    wt_mut_seq = np.tile(wild_type_int_wt, (len(variants), 1))
    for i, variant in enumerate(variants):
        if variant == "_wt":
            continue
        variant = variant.split(",")
        for mutation in variant:
            position = int(mutation[1:-1])
            replacement = constants.C2I_MAPPING[mutation[-1]]
            seq_ints[i, position-wt_offset] = replacement
            wt_mut_seq[i, position-wt_offset] = replacement

    return seq_ints, wt_mut_seq


def encode_int_seqs(char_seqs=None, variants=None, wild_type_aa=None, wild_type_offset=None):
    """
    encodes either character sequences or variants to integer sequences
    param: 
        char_seqs: list of character sequences
        variants: list of variants
        wild_type_aa: wild-type character sequence
        wild_type_offset: wild-type offset
    return: 
        int_seqs: numpy array of integer encoded sequences
        wt_minus_mt_seqs: numpy array of integer encoded wild-type minus mutant sequences
        single: boolean indicating if a single sequence was provided
    """
    single = False
    if variants is not None:
        if not isinstance(variants, list):
            single = True
            variants = [variants]
        int_seqs, wt_minus_mt_seqs = enc_int_seqs_from_variants(variants, wild_type_aa, wild_type_offset)

    elif char_seqs is not None:
        if not isinstance(char_seqs, list):
            single = True
            char_seqs = [char_seqs]

        int_seqs = enc_int_seqs_from_char_seqs(char_seqs)
    return int_seqs, wt_minus_mt_seqs, single


def encode(args, encoding, char_seqs=None, variants=None, ds_name=None, wt_aa=None, wt_offset=None):
    """ the main encoding function that will encode the given sequences or variants and return the encoded data """

    if variants is None and char_seqs is None:
        raise ValueError("must provide either variants or full sequences to encode")
    if variants is not None and ((ds_name is None) and ((wt_aa is None) or (wt_offset is None))):
        raise ValueError("if providing variants, must also provide (wt_aa and wt_offset) or "
                         "ds_name so I can look up the WT sequence myself")

    if ds_name is not None:
        wt_aa = constants.DATASETS[ds_name]["wt_aa"]
        wt_offset = constants.DATASETS[ds_name]["wt_ofs"]

    # get integer encoded sequences from either char seqs or variants
    int_seqs, wt_minus_mt_seqs, single = encode_int_seqs(char_seqs=char_seqs, variants=variants,
                                       wild_type_aa=wt_aa, wild_type_offset=wt_offset)

    # encode variants using int seqs
    encodings = encoding.split(",")
    encoded_data = []
    for enc in encodings:
        if enc == "one_hot":
            encoded_data.append(enc_one_hot(int_seqs))
        elif enc == "aa_index":
            encoded_data.append(enc_aa_index(wt_minus_mt_seqs))
        elif enc == "rmsf":
            encoded_data = np.concatenate(encoded_data, axis=-1)
            encoded_rmsf = enc_rmsf(args.rmsf_file, wt_minus_mt_seqs)
            encoded_data = np.concatenate((encoded_data,encoded_rmsf), axis=2)
        elif enc == "rosetta":
            encoded_data = enc_rosetta(args.rosetta_file, variants,encoded_data,wt_offset)
        else:
            raise ValueError("err: encountered unknown encoding: {}".format(enc))

    if single:
        encoded_data = encoded_data[0]

    return encoded_data


def encode_full_dataset(ds_name, encoding, args):
    # load the dataset
    ds_fn = constants.DATASETS[ds_name]["ds_fn"]
    ds =  pd.read_csv(ds_fn, sep="\t")
    # encode the data
    encoded_data = encode(args, encoding=encoding, variants=ds["variant"].tolist(), ds_name=ds_name)
    return encoded_data


def encode_full_dataset_and_save(ds_name, encoding, args):
    """ encoding a full dataset """
    out_fn = join(constants.DATASETS[ds_name]["ds_dir"], "enc_{}_{}.npy".format(ds_name, encoding))
    if isfile(out_fn):
        print("err: encoded data already exists: {}".format(out_fn))
        return
    encoded_data = encode_full_dataset(ds_name, encoding, args)
    np.save(out_fn, encoded_data)
    return encoded_data


def main(args):

    if args.ds_name == "all":
        ds_names = constants.DATASETS.keys()
    else:
        ds_names = [args.ds_name]

    if args.encoding == "all":
        encodings = ["one_hot", "aa_index", "rmsf", "rosetta"]
    else:
        encodings = [args.encoding]

    for ds_name in ds_names:
        for encoding in encodings:
            encode_full_dataset_and_save(ds_name, encoding, args)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("ds_name",
                        help="name of the dataset",
                        type=str)
    parser.add_argument("encoding",
                        help="what encoding to use",
                        type=str,
                        choices=["one_hot", "aa_index", "rmsf", "rosetta", "all"])    
    main(parser.parse_args())
