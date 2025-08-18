#!/bin/bash
# ------------------------------------------------------------------------------
## 03_measure_ari_nmi.sh
# This script estimates the batch effect removal and bio conservations ARI and
# NMI. It uses the estimated Leiden clusters to calculate the overlap of the
# estimated clusters with the celltype and batch.
# ------------------------------------------------------------------------------
#SBATCH --job-name=ari_nmi          # Job name    (default: sbatch)
#SBATCH --output=ari_nmi-%j.out     # Output file (default: slurm-%j.out)
#SBATCH --error=ari_nmi-%j.err      # Error file  (default: slurm-%j.out)
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

# --------------------------------------------------------------------------
# Measure ARI, NMI, and harmonized Shannon entropy
echo "Measurement for ARI and NMI, started at `date`"
Rscript src/benchmark/measure_ari_nmi.R \
        results/benchmark/leiden_clusters/oma.csv \
        oma_ari_nmi.csv
Rscript src/benchmark/measure_ari_nmi.R \
        results/benchmark/leiden_clusters/mmseqs2_hog.csv \
        mmseqs2_hog_ari_nmi.csv
Rscript src/benchmark/measure_ari_nmi.R \
        results/benchmark/leiden_clusters/mmseqs2_oog.csv \
        mmseqs2_oog_ari_nmi.csv
Rscript src/benchmark/measure_ari_nmi.R \
        results/benchmark/leiden_clusters/diamond_hog.csv \
        diamond_hog_ari_nmi.csv
Rscript src/benchmark/measure_ari_nmi.R \
        results/benchmark/leiden_clusters/diamond_oog.csv \
        diamond_oog_ari_nmi.csv
Rscript src/benchmark/measure_ari_nmi.R \
        results/benchmark/leiden_clusters/samap_2.0.csv \
        samap_ari_nmi.csv

# --------------------------------------------------------------------------
# Measure harmonized Shannon Entropy
echo "Measurement for harmonized Shannon Entropy, started at `date`"
Rscript src/benchmark/measure_entropy.R \
        results/benchmark/leiden_clusters/oma.csv \
        oma_entropy.csv
Rscript src/benchmark/measure_entropy.R \
        results/benchmark/leiden_clusters/mmseqs2_hog.csv \
        mmseqs2_hog_entropy.csv
Rscript src/benchmark/measure_entropy.R \
        results/benchmark/leiden_clusters/mmseqs2_oog.csv \
        mmseqs2_oog_entropy.csv
Rscript src/benchmark/measure_entropy.R \
        results/benchmark/leiden_clusters/diamond_hog.csv \
        diamond_hog_entropy.csv
Rscript src/benchmark/measure_entropy.R \
        results/benchmark/leiden_clusters/diamond_oog.csv \
        diamond_oog_entropy.csv
Rscript src/benchmark/measure_entropy.R \
        results/benchmark/leiden_clusters/samap_2.0.csv \
        samap_entropy.csv

echo "Measurement finished at `date`"