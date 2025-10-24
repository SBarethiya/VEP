# Variant effect predictor (VEP)
This project integrates biophysics-based features from Rosetta into existing variant effect prediction models to address data scarcity and improve the robustness of mutational impact predictions for mutational and positional extrapolation.
By combining computational protein modeling with machine learning, VEP enhances predictive performance for protein variant effects.

## Rosetta data preparation

There are two procedures to obtain Rosetta-derived features used in this project:
1. **[FastRelax](https://github.com/SBarethiya/VEP/tree/main/Rosetta/FastRelax)**: 
    - Used for structure relaxation and energy minimization.
    - Require [Rosetta installation](https://docs.rosettacommons.org/demos/latest/tutorials/install_build/install_build)

2. **[PackRotamer](https://github.com/SBarethiya/VEP/tree/main/Rosetta/PackRotamerMover)**: 
    - Used for sidechain repacking and biophysics-based features.
    - Requires [PyRosetta installation](https://www.pyrosetta.org/downloads#h.ydwwhv85t3cc) 

Please refer to the [Rosetta directory](https://github.com/SBarethiya/VEP/tree/main/Rosetta) for scripts and details for both procedures.

## Prerequisites
Before running the code, ensure the following dependencies are installed:

```bash
conda create -n vep-env python=3.6
conda activate vep-env

pip install tensorflow==1.14

# if GPU (avaiable)
pip install tensorflow-gpu==1.14

conda install cudatoolkit=10.0.130
pip install -r requirements.txt
```
## Usage
```
# Train the model
python code/regression --dataset_name bgl3 --dataset_file data/bgl3/bgl3.tsv --encoding one_hot,aa_index,rmsf,rosetta --net_file network_specs/lr.yml --epochs 2 --rmsf_file data/bgl3/rmsf.tsv  --rosetta_file data/bgl3/ddG.tsv
```
Modify configuration files under network_specs/ to adjust hyperparameters or Python code/regression -h

## Example Workflow
1. Prepare the new dataset: create new directory under [data/](https://github.com/SBarethiya/VEP/tree/main/data), put a relaxed PDB, add information about the protein in [code/constants.py](https://github.com/SBarethiya/VEP/tree/main/code/constants.py)
2. Prepare Rosetta features (FastRelax, PackRotamerMover) and RMSF features and put under the directory data/your_protein.
3. Create splits: Either you can create new splits by providing the arguments during training (for example: --train_size 0.7 --tune_size 0.15 --test_size 0.15) or you can create new splits by using: [code/train_test_split.py](https://github.com/SBarethiya/VEP/tree/main/code/train_test_split.py). 
    ```
    python code/train_test_split.py --protein bgl3 --split_type standard --seeds 42 --ratios 0.7 0.15 0.15
    ```
    If you are using GCN then also create: 
    ```
    python code/gen_structure_graph.py --datasets bgl3 --thresholds 6 7 8 --graph-types DIST_THRESH
    ```
4. Train or fine-tune the neural network model

## Acknowledgment
The code and networks are adapted from [nn4dms](https://github.com/gitter-lab/nn4dms) developed by the Gitter Lab. 
We appreciate their contributions for the open-source tools.
