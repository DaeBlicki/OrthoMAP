#!/bin/bash
# ------------------------------------------------------------------------------
## 05_run_orthomap.sh
# Run OrthoMAP workflow. This script generates sparse cell vs. cell distance 
# matrices along all species used in the orthology interference part. The
# results are UMAP plots and finalized Seurat Objects (v5.0.1).
#
# When this script is used to replicate the results, make sure that the right
# configuration is loaded in advanced! 
# 
# Use: source("scripts/05a_load_oma_config.sh")
#   - 05a: OMA analysis on OMA groups (also called Orthologous Groups)
#   - 05b: OrthoFinder MMseqs2 on HOG (Hierarchical Orthologous Groups)
#   - 05c: OrthoFinder MMseqs2 on OOG (1:1 Orthologs)
#   - 05d: OrthoFinder Diamond on HOG
#   - 05e: OrthoFinder Diamond on OOG
# ------------------------------------------------------------------------------
#SBATCH --job-name=orthomap             # Job name    (default: sbatch)
#SBATCH --output=orthomap-%j.out        # Output file (default: slurm-%j.out)
#SBATCH --error=orthomap-%j.err         # Error file  (default: slurm-%j.out)
#SBATCH --constraint=amd                # Target AMD CPU (7452 or 9554)
#SBATCH --partition=basic               # Unlimited time (no hpc project)
#SBATCH --ntasks=1                      # Number of tasks
#SBATCH --cpus-per-task=1               # Number of CPUs per task
#SBATCH --mem-per-cpu=50GB              # Memory per CPU (in MB)
#SBATCH --time=03:00:00                 # Wall clock time limit (H:M:S)
#SBATCH --constraint=amd                # Target AMD CPU (7452 or 9554)
#SBATCH --partition=basic               # Unlimited time (no hpc project)

# load module and print information 
module load R-bundle-CRAN/2024.11-foss-2024a    # Rcpp
module load seurat5             # load seurat5 and dependency
module list                     # list the loaded modules
lscpu | grep "Model name"       # CPU Architecture
date                            # Time of submitted job

# Run the pipeline
Rscript src/OrthoMAP_main.R x