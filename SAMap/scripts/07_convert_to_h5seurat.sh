#!/bin/bash
# ------------------------------------------------------------------------------
## 05_convert_to_h5seurat.sh
# Convert the single-cell Atlas from Python AnnData to R Seurat Object.
# ------------------------------------------------------------------------------
#SBATCH --job-name=convert2             # Job name    (default: sbatch)
#SBATCH --output=convert2-%j.out        # Output file (default: slurm-%j.out)
#SBATCH --error=convert2-%j.err         # Error file  (default: slurm-%j.out)
#SBATCH --constraint=amd                # Target AMD CPU (7452 or 9554)
#SBATCH --partition=basic               # Unlimited time (no hpc project)
#SBATCH --ntasks=1                      # Number of tasks
#SBATCH --cpus-per-task=1               # Number of CPUs per task
#SBATCH --mem-per-cpu=32GB               # Memory per CPU (in MB)
#SBATCH --time=04:00:00                 # Wall clock time limit (H:M:S)
#SBATCH --constraint=amd                # Target AMD CPU (7452 or 9554)
#SBATCH --partition=basic               # Unlimited time (no hpc project)

# load module and print information 
module load python3 conda seurat5   # list the loaded modules
module list                         # list the loaded modules
lscpu | grep "Model name"           # CPU Architecture
date                                # Time of the submitted job

conda activate SAMap

# ------------------------------------------------------------------------------
## Create Metadata for the SAMap object
echo "[R] Create SAMap Object Metadata, started at `date`"
python3 src/create_samap_metadata.R

# ------------------------------------------------------------------------------
## Convert SAMap to .h5ad file
echo "[Python3] SAMap -> h5ad, started at `date`"
python3 src/create_samap_h5ad_object.py

# ------------------------------------------------------------------------------
## Convert .h5ad file to .h5seurat
echo "[R] h5ad -> h5seurat, started at `date`"
Rscript src/create_samap_h5seurat_object.R

echo "Conversion of SAMap finished at `date`"