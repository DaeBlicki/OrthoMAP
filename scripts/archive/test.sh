#!/bin/bash
#SBATCH --job-name=silhouette        # Job name    (default: sbatch)
#SBATCH --output=silhouette-%j.out   # Output file (default: slurm-%j.out)
#SBATCH --error=silhouette-%j.err    # Error file  (default: slurm-%j.out)
#SBATCH --constraint=amd              # Target AMD CPU (7452 or 9554)
#SBATCH --partition=basic             # Unlimited time (no hpc project)
#SBATCH --ntasks=1                    # Number of tasks
#SBATCH --cpus-per-task=1            # Number of CPUs per task
#SBATCH --mem-per-cpu=200GB             # Memory per CPU (in MB)
#SBATCH --time=1-00:00:00             # Wall clock time limit (H:M:S)
#SBATCH --constraint=amd              # Target AMD CPU (7452 or 9554)
#SBATCH --partition=basic             # Unlimited time (no hpc project)

# load module and print information 
module load R-bundle-CRAN/2024.11-foss-2024a    # Rcpp
module load seurat5             # load seurat5 and dependency
module list                     # list the loaded modules
lscpu | grep "Model name"       # CPU Architecture
date                            # Time of submitted job

# Run the pipeline
Rscript src/archive/silhouette_test.R