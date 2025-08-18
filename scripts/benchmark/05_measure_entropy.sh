#!/bin/bash
# ------------------------------------------------------------------------------
## 05_measure_ari_nmi.sh
# This script estimates the normalized Shannon entropy for each estimated
# clusters in dependency of the celltype ("IDs").
# ------------------------------------------------------------------------------
#SBATCH --job-name=entropy          # Job name    (default: sbatch)
#SBATCH --output=entropy-%j.out     # Output file (default: slurm-%j.out)
#SBATCH --error=entropy-%j.err      # Error file  (default: slurm-%j.out)
#SBATCH --constraint=amd            # Target AMD CPU (7452 or 9554)
#SBATCH --partition=basic           # Unlimited time (no hpc project)
#SBATCH --ntasks=1                  # Number of tasks
#SBATCH --cpus-per-task=1           # Number of CPUs per task
#SBATCH --mem-per-cpu=8GB           # Memory per CPU (in MB)
#SBATCH --time=00:30:00             # Wall clock time limit (H:M:S)
#SBATCH --constraint=amd            # Target AMD CPU (7452 or 9554)
#SBATCH --partition=basic           # Unlimited time (no hpc project)
#SBATCH --partition=basic           # Unlimited time (no hpc project)

# load module and print information 
module load R-bundle-CRAN/2024.11-foss-2024a    # Rcpp
module load seurat5             # load seurat5 and dependency
module list                     # list the loaded modules
lscpu | grep "Model name"       # CPU Architecture
date                            # Time of submitted job

# Run the measurements
echo "Measurement for OMA `date`"
Rscript src/benchmark/measure_entropy.R \
        results/benchmark/leiden_clusters/oma.csv \
        oma_entropy.csv

# Run the measurements
echo "Measurement for MMseqs2 OOG `date`"
Rscript src/benchmark/measure_entropy.R \
        results/benchmark/leiden_clusters/mmseqs2_oog.csv \
        mmseqs2_oog_entropy.csv

# Run the measurements
echo "Measurement for Diamond OOG `date`"
Rscript src/benchmark/measure_entropy.R \
        results/benchmark/leiden_clusters/diamond_oog.csv \
        diamond_oog_entropy.csv

# Run the measurements
echo "Measurement for SAMAp `date`"
Rscript src/benchmark/measure_entropy.R \
        results/benchmark/leiden_clusters/samap_best.csv \
        samap_entropy.csv

