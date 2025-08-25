# Variant effect predictor (VEP)

This project focuses on adding the physics-based features in the existing models to improve the data scarity issues.

## Prerequisites
Before running the code, ensure the following dependencies are installed:
```bash
conda create -n vep-env python=3.6
conda activate vep-env

conda install cudatoolkit=10.0.130
pip install -r requirements.txt
```

## Rosetta data preparation

There are two procedure to get the Rosetta scores:
1. **[FastRelax](https://github.com/SBarethiya/VEP/tree/main/Rosetta/FastRelax)**
2. **[PackRotamer](https://github.com/SBarethiya/VEP/tree/main/Rosetta/PackRotamerMover)**

[PyRosetta installation](https://www.pyrosetta.org/downloads#h.ydwwhv85t3cc) 
[Rosetta installation](https://docs.rosettacommons.org/demos/latest/tutorials/install_build/install_build)

This repository contains script for these two procedure see in [Rosetta](https://github.com/SBarethiya/VEP/tree/main/Rosetta)

## Acknowledgment
The code is adapted from the [nn4dms](https://github.com/gitter-lab/nn4dms) developed by the Gitter Lab. 
We appreciate their contributions for the open-source tools.
