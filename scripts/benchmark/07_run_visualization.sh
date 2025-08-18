#!/bin/bash
# ------------------------------------------------------------------------------
## 07_run_visualization.sh
# This script visualize the results for ARI, NMI, ASW, LISI and harmonized
# Shannon Entropy. In addition, it runs and creates the UMAP plots
# ------------------------------------------------------------------------------
#SBATCH --job-name=visualize            # Job name    (default: sbatch)
#SBATCH --output=visualize-%j.out       # Output file (default: slurm-%j.out)
#SBATCH --error=visualize-%j.err        # Error file  (default: slurm-%j.out)
#SBATCH --constraint=amd                # Target AMD CPU (7452 or 9554)
#SBATCH --partition=basic               # Unlimited time (no hpc project)
#SBATCH --ntasks=1                      # Number of tasks
#SBATCH --cpus-per-task=1               # Number of CPUs per task
#SBATCH --mem-per-cpu=64GB              # Memory per CPU (in MB)
#SBATCH --time=04:00:00                 # Wall clock time limit (H:M:S)
#SBATCH --constraint=amd                # Target AMD CPU (7452 or 9554)
#SBATCH --partition=basic               # Unlimited time (no hpc project)

# load module and print information 
module load R-bundle-CRAN/2024.11-foss-2024a    # Rcpp
module load seurat5             # load seurat5 and dependency
module list                     # list the loaded modules
lscpu | grep "Model name"       # CPU Architecture
date                            # Time of submitted job

# ------------------------------------------------------------------------------
# Shannon Entropy plots in form of donut charts
echo "Run Donut Charts, started at `date`"
Rscript src/benchmark/create_donut_charts.R

# ------------------------------------------------------------------------------
# UMAP plots
echo "Run UMAP, started at `date`"
Rscript src/benchmark/create_umap_results.R

echo "Finished at `date`"