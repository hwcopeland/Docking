#!/bin/bash

# Create a new conda environment
conda create -n docking python=3.7

# Activate the environment
source activate docking

apt install -y autodock-vina
# Install conda packages
conda install -c conda-forge fpocket
conda install -c bioconda vina
conda install -c schrodinger pymol
conda install -c conda-forge biopython
conda install -c conda-forge openbabel

pip install rdkit

export PATH="$HOME/miniconda/envs/docking/bin:$PATH"
export PATH="$HOME/miniconda/bin:$PATH"0

mkdir Ligands Protein Results
cd Ligands
mkdir PDB PDBQT
cd ..