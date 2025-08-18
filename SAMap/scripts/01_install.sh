#!/bin/bash
# ------------------------------------------------------------------------------
## 01_install.sh
# SAMap requires BLAST results between all pairs of species. This script
# estimates the BLAST results using BLAST.
#
# usage: source scripts/01_install.sh
# ------------------------------------------------------------------------------

# load module and print information 
module load conda               # load OMA (standalone version)
module list                     # list the loaded modules
lscpu | grep "Model name"       # CPU Architecture
date                            # Time of the submitted jobmodule load conda

# Install SAMap dependencies availabe in conda
conda create -n SAMap -c conda-forge python=3.9 numpy=1.23.5 pip pybind11 h5py=3.8.0 leidenalg python-igraph texttable
conda activate SAMap

# Install GitRepository (version 1.0.15)
git clone https://github.com/atarashansky/SAMap.git samap_directory
cd samap_directory
pip install .

# Define NCBI BLAST version.
ncbi_blast_version='2.9.0'

# Download NCBI BLAST tarball.
wget "ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/${ncbi_blast_version}/ncbi-blast-${ncbi_blast_version}+-x64-linux.tar.gz"

# Extract NCBI BLAST binaries in current conda environment bin directory.
tar -xzvf "ncbi-blast-${ncbi_blast_version}+-x64-linux.tar.gz" \
    -C "${CONDA_PREFIX}/bin/" \
    --strip-components=2 \
    "ncbi-blast-${ncbi_blast_version}+/bin/"

cd ..
