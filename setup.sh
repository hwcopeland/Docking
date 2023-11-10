#!/bin/bash

# Install conda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/miniconda

# Add conda to PATH
export PATH="$HOME/miniconda/bin:$PATH"

# Create a new conda environment
conda create -n docking python=3.7

# Activate the environment
source activate docking

# Install rdkit and fpocket with conda
conda install -c rdkit rdkit
conda install -c bioconda fpocket

# Install autodock vina
conda install -c bioconda autodock-vina

# Install biopython and openbabel with pip
pip install biopython openbabel

cd $HOME

mkdir Docking
mkdir Ligands
mkdir Protein
mkdir Results
cd Docking
wget https://github.com/hwcopeland/Docking/blob/main/docking.py > docking.py