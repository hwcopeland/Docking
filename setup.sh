#!/bin/bash

# Check the operating system
OS="$(uname)"
MACHINE="$(uname -m)"

if [ "$OS" == "Linux" ]; then
    # Install conda
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/miniconda

    # Add conda to PATH
    export PATH="$HOME/miniconda/bin:$PATH"
elif [ "$OS" == "Darwin" ]; then
    if [ "$MACHINE" == "arm64" ]; then
        # Install conda for ARM
        wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-arm64.sh
        bash Miniforge3-MacOSX-arm64.sh -b -p $HOME/miniconda
    else
        # Install conda for Intel
        wget https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
        bash Miniconda3-latest-MacOSX-x86_64.sh -b -p $HOME/miniconda
    fi

    # Add conda to PATH
    export PATH="$HOME/miniconda/bin:$PATH"
fi

# Create a new conda environment
conda create -n docking python=3.7

# Activate the environment
source activate docking

# Install conda packages
conda install -c conda-forge rdkit
conda install -c conda-forge fpocket
conda install -c conda-forge vina
conda install -c schrodinger pymol
conda install -c conda-forge biopython
conda install -c conda-forge openbabel

export PATH="$HOME/miniconda/envs/docking/bin:$PATH"
export PATH="$HOME/miniconda/bin:$PATH"0

mkdir Ligands Protein Results
cd Ligands
mkdir PDB PDBQT
cd ..