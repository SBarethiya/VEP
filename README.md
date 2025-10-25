# Variant effect predictor (VEP)
This project integrates biophysics-based features from Rosetta into existing variant effect prediction models to address data scarcity and improve the robustness of mutational impact predictions for mutational and positional extrapolation.

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

### Conda Environment
```bash
conda create -n vep-env python=3.6
conda activate vep-env
```
### Python Packages
```
# CPU-only
pip install tensorflow==1.14

# GPU-enabled (optional)
pip install tensorflow-gpu==1.14

# Other packages
pip install -r requirements.txt
```
## Usage
### Train the model
```
python code/regression.py \
    --dataset_name bgl3 \
    --dataset_file data/bgl3/bgl3.tsv \
    --encoding one_hot,aa_index,rmsf,rosetta \
    --net_file network_specs/lr.yml \
    --epochs 2 \
    --rmsf_file data/bgl3/rmsf.tsv \
    --rosetta_file data/bgl3/ddG.tsv
```
- Modify hyperparameters in network_specs/ or see hyperparameters or Python code/regression -h
- Input files should be tab-delimited .tsv files.

## Example Workflow
1. Prepare the new dataset: 
    - create new directory under [data/](https://github.com/SBarethiya/VEP/tree/main/data)
    - Place relaxed PDBs and RMSF/ddG data as .tsv files.
    - Update protein information in [code/constants.py](https://github.com/SBarethiya/VEP/tree/main/code/constants.py).
2. Prepare Rosetta Features 
    - Run FastRelax and PackRotamerMover to generate features.
    - Include RMSF from MD simulations.
3. Create splits: 
    - Either you can create new splits by providing the arguments during training (for example: --train_size 0.7 --tune_size 0.15 --test_size 0.15) 
    - Or you can create new splits by using: [code/train_test_split.py](https://github.com/SBarethiya/VEP/tree/main/code/train_test_split.py). 
        ```
        python code/train_test_split.py --protein bgl3 --split_type standard --seeds 42 --ratios 0.7 0.15 0.15
        ```
    - For GCN models, also generate structure graphs: 
        ```
        python code/gen_structure_graph.py --datasets bgl3 --thresholds 6 7 8 --graph-types DIST_THRESH
        ```
4. Train or Fine-tune Models
    - Use the training script with desired dataset and encoding.
    - A few example network configurations are provided in the `network_specs/` directory. For more examples, refer to [nn4dms](https://github.com/gitter-lab/nn4dms).

## Acknowledgment
The code and networks are adapted from [nn4dms](https://github.com/gitter-lab/nn4dms) developed by the Gitter Lab. 
We appreciate their contributions for the open-source tools.
