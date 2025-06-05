#!/bin/bash

# Create the list of datasets and corresponding net_files
my_list=("avgfp" "bgl3" "gb1" "pab1" "ube4b")
net_files=("cnn-5xk3f128.yml" "cnn-1xk17f32.yml" "cnn-3xk17f128.yml" "cnn-3xk17f128.yml" "cnn-5xk3f128.yml")
batch_sizes=(64 64 64 128 128)
learning_rates=(0.001 0.001 0.001 0.0001 0.0001)

# Loop through each dataset and corresponding net_file using index
for i in "${!my_list[@]}"; do
    protein="${my_list[$i]}"
    netfile="${net_files[$i]}"
    batch_size="${batch_sizes[$i]}"
    learning_rate="${learning_rates[$i]}"
    for split_dir in data/$protein/splits/*/; do
        split_dir="${split_dir%/}"
        out="${split_dir##*/}"
        python code/regression.py \
            --dataset_name "$protein" \
            --dataset_file "data/$protein/$protein.txt" \
            --net_file "network_specs/cnns/${netfile}" \
            --encoding "one_hot,aa_index" \
            --batch_size "$batch_size" \
            --learning_rate "$learning_rate" \
            --epochs 500 \
            --early_stopping \
            --split_dir "$split_dir" \
            --log_dir_base "output/${protein}/cnn/${out}" \
            --rmsf_file "data/$protein/rmsf.tsv" \
            --ss_file "data/$protein/ss.npy" \
            --rosetta_file "data/$protein/ddG.tsv"
    done
done

my_list=("avgfp" "bgl3" "gb1" "pab1" "ube4b")
net_files=("gc-2xf128.yml" "gc-1xf32.yml" "gc-5xf128.yml" "gc-1xf32.yml" "gc-3xf128.yml")
graphs=("dist_thresh_7.graph" "dist_thresh_6.graph" "dist_thresh_7.graph" "dist_thresh_7.graph" "dist_thresh_7.graph")
batch_sizes=(32 128 32 128 128)

# Loop through each dataset and corresponding net_file using index
for i in "${!my_list[@]}"; do
    protein="${my_list[$i]}"
    netfile="${net_files[$i]}"
    graph="${graphs[$i]}"
    batch_size="${batch_sizes[$i]}"
    for split_dir in data/$protein/splits/*/; do
        split_dir="${split_dir%/}"
        out="${split_dir##*/}"
        python code/regression.py \
            --dataset_name "$protein" \
            --dataset_file "data/${protein}/${protein}.txt" \
            --net_file "network_specs/gcns/${netfile}" \
            --encoding "one_hot,aa_index" \
            --graph_fn "data/${protein}/graphs/${graph}" \
            --batch_size "$batch_size" \
            --learning_rate 0.0001 \
            --epochs 500 \
            --early_stopping \
            --split_dir "$split_dir" \
            --log_dir_base "output/${protein}/gcn${out}" \
            --rmsf_file "data/$protein/rmsf.tsv" \
            --ss_file "data/$protein/ss.npy" \
            --rosetta_file "data/$protein/ddG.tsv"
    done
done

