#!/bin/bash
# ------------------------------------------------------------------------------
## 02b_measure_mmseqs2_hog.sh
# This script load the MMseqs2 HOG configuration and then estimate the OrthoMAP
# clusters for the random state: 1 to 20 (Leiden Algorithm use random seed!). 
# The clusters are stored with the cell ids as .csv
# 
#   - 02a: OMA analysis on OMA groups (also called Orthologous Groups)
#   - 02b: OrthoFinder MMseqs2 on HOG (Hierarchical Orthologous Groups)
#   - 02c: OrthoFinder MMseqs2 on OOG (1:1 Orthologs)
#   - 02d: OrthoFinder Diamond on HOG
#   - 02e: OrthoFinder Diamond on OOG
# ------------------------------------------------------------------------------
#SBATCH --job-name=measure            # Job name    (default: sbatch)
#SBATCH --output=measure-%j.out       # Output file (default: slurm-%j.out)
#SBATCH --error=measure-%j.err        # Error file  (default: slurm-%j.out)
#SBATCH --constraint=amd                # Target AMD CPU (7452 or 9554)
#SBATCH --partition=basic               # Unlimited time (no hpc project)
#SBATCH --ntasks=1                      # Number of tasks
#SBATCH --cpus-per-task=1               # Number of CPUs per task
#SBATCH --mem-per-cpu=200GB              # Memory per CPU (in MB)
#SBATCH --time=12:00:00                 # Wall clock time limit (H:M:S)
#SBATCH --constraint=amd                # Target AMD CPU (7452 or 9554)
#SBATCH --partition=basic               # Unlimited time (no hpc project)

# load module and print information 
module load R-bundle-CRAN/2024.11-foss-2024a    # Rcpp
module load seurat5             # load seurat5 and dependency
module list                     # list the loaded modules
lscpu | grep "Model name"       # CPU Architecture
date                            # Time of submitted job

# Run the pipeline
Rscript src/benchmark/measure_orthomap_cluster.R \
        data/configuration/05b_config_mmseqs2_hog.yaml \
        mmseqs2_hog.csv