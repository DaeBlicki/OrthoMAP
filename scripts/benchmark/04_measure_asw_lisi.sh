#!/bin/bash
# ------------------------------------------------------------------------------
## 04_measure_asw_lisi.sh
# This script estimates the batch effect removal and bio conservations ASW and
# LISI. It uses the 80% of the cells to calculate intrinsic conservation of the
# biological variability and batch score. This is done 20 times.
# ------------------------------------------------------------------------------
#SBATCH --job-name=asw_lisi         # Job name    (default: sbatch)
#SBATCH --output=asw_lisi-%j.out    # Output file (default: slurm-%j.out)
#SBATCH --error=asw_lisi-%j.err     # Error file  (default: slurm-%j.out)
#SBATCH --constraint=amd            # Target AMD CPU (7452 or 9554)
#SBATCH --partition=basic           # Unlimited time (no hpc project)
#SBATCH --ntasks=1                  # Number of tasks
#SBATCH --cpus-per-task=1           # Number of CPUs per task
#SBATCH --mem-per-cpu=80GB           # Memory per CPU (in MB)
#SBATCH --time=24:00:00             # Wall clock time limit (H:M:S)
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
Rscript src/benchmark/measure_asw_lisi.R \
        results/OrthoMAP_Seurat_Objects/oma/result.Robj \
        oma_asw_lisi.csv

# Run the measurements
echo "Measurement for MMseqs2 HOG `date`"
Rscript src/benchmark/measure_asw_lisi.R \
        results/OrthoMAP_Seurat_Objects/mmseqs2/HOG/result.Robj \
        mmseqs2_hog_asw_lisi.csv

# Run the measurements
echo "Measurement for MMseqs2 OOG `date`"
Rscript src/benchmark/measure_asw_lisi.R \
        results/OrthoMAP_Seurat_Objects/mmseqs2/OOG/result.Robj \
        mmseqs2_oog_asw_lisi.csv

# Run the measurements
echo "Measurement for Diamond HOG `date`"
Rscript src/benchmark/measure_asw_lisi.R \
        results/OrthoMAP_Seurat_Objects/diamond/HOG/result.Robj \
        diamond_hog_asw_lisi.csv

# Run the measurements
echo "Measurement for Diamond OOG `date`"
Rscript src/benchmark/measure_asw_lisi.R \
        results/OrthoMAP_Seurat_Objects/diamond/OOG/result.Robj \
        diamond_oog_asw_lisi.csv

echo "Measurement finished at `date`"
