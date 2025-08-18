#!/bin/bash
# ------------------------------------------------------------------------------
## 02_convert_to_h5ad.sh
# Convert the single-cell Atlas from R Seurat to Object Python AnnData Object.
# ------------------------------------------------------------------------------
#SBATCH --job-name=convert1             # Job name    (default: sbatch)
#SBATCH --output=convert1-%j.out        # Output file (default: slurm-%j.out)
#SBATCH --error=convert1-%j.err         # Error file  (default: slurm-%j.out)
#SBATCH --constraint=amd                # Target AMD CPU (7452 or 9554)
#SBATCH --partition=basic               # Unlimited time (no hpc project)
#SBATCH --ntasks=1                      # Number of tasks
#SBATCH --cpus-per-task=1               # Number of CPUs per task
#SBATCH --mem-per-cpu=8GB               # Memory per CPU (in MB)
#SBATCH --time=00:15:00                 # Wall clock time limit (H:M:S)
#SBATCH --constraint=amd                # Target AMD CPU (7452 or 9554)
#SBATCH --partition=basic               # Unlimited time (no hpc project)

# load module and print information
module load seurat5 python3 conda
module list                     # list the loaded modules
lscpu | grep "Model name"       # CPU Architecture
date                            # Time of the submitted jobmodule load conda

conda activate SAMap

# ------------------------------------------------------------
# Hydra vulgaris
echo "Hydra vulgaris single-cell atlas, started at `date`"
Rscript src/convert_to_h5ad_object.R \
        data/single_cell_atlas/Hv.Robj \
        data/single_cell_atlas/Hv.h5Seurat
rm -r data/single_cell_atlas/Hv.h5Seurat

# ------------------------------------------------------------
# Aurelia coerulea
echo "Aurelia coerulea single-cell atlas, started `date`"
Rscript src/convert_to_h5ad_object.R \
        data/single_cell_atlas/Ac.Robj \
        data/single_cell_atlas/Ac.h5Seurat
rm -r data/single_cell_atlas/Ac.h5Seurat

# ------------------------------------------------------------
# Nematostella vectensis
echo "Nematostella vectensis single-cell atlas, started `date`"
Rscript src/convert_to_h5ad_object.R \
        data/single_cell_atlas/Nv.Robj \
        data/single_cell_atlas/Nv.h5Seurat
rm -r data/single_cell_atlas/Nv.h5Seurat

# ------------------------------------------------------------
echo "Conversion finished at `date`"