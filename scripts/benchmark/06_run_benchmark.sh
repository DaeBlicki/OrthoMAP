#!/bin/bash
# ------------------------------------------------------------------------------
## 06_run_benchmark.sh
# This script run the benchmark using the ARI, NMI, ASW, LISI and harmonized
# Shannon Entropy. The batch and bio conservations are plotted as Lolipop plots.
# The entropy is presented as boxplot.
# ------------------------------------------------------------------------------
#SBATCH --job-name=benchmark            # Job name    (default: sbatch)
#SBATCH --output=benchmark-%j.out       # Output file (default: slurm-%j.out)
#SBATCH --error=benchmark-%j.err        # Error file  (default: slurm-%j.out)
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

Rscript src/benchmark/compare_configuration.R
Rscript src/benchmark/compare_resolution.R
Rscript src/benchmark/benchmark_measurements.R
# Rscript src/benchmark/benchmark_best_results.R

echo "Finished at `date`"