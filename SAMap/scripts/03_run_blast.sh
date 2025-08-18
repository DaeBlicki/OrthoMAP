#!/bin/bash
# ------------------------------------------------------------------------------
## 03_run_blast.sh
# SAMap requires BLAST results between all pairs of species. This script
# estimates the BLAST results using BLAST.
# ------------------------------------------------------------------------------
#SBATCH --job-name=blast                # Job name    (default: sbatch)
#SBATCH --output=blast-%j.out           # Output file (default: slurm-%j.out)
#SBATCH --error=blast-%j.err            # Error file  (default: slurm-%j.out)
#SBATCH --constraint=amd                # Target AMD CPU (7452 or 9554)
#SBATCH --partition=basic               # Unlimited time (no hpc project)
#SBATCH --ntasks=1                      # Number of tasks
#SBATCH --cpus-per-task=64              # Number of CPUs per task
#SBATCH --mem-per-cpu=1GB               # Memory per CPU (in MB)
#SBATCH --time=12:00:00                 # Wall clock time limit (H:M:S)
#SBATCH --constraint=amd                # Target AMD CPU (7452 or 9554)
#SBATCH --partition=basic               # Unlimited time (no hpc project)

# load module and print information 
module list                     # list the loaded modules
lscpu | grep "Model name"       # CPU Architecture
date                            # Time of the submitted jobmodule load conda

# move to samap_directory
cd samap_directory

## Run BLAST results
# Ac vs. Hv
echo "[BLAST] Ac vs. Hv started at `date`"
bash map_genes.sh --tr1 ../data/DB/Ac.fa \
                  --t1 prot \
                  --n1 Ac \
                  --tr2 ../data/DB/Hv.fa \
                  --t2 prot \
                  --n2 Hv \
# Ac vs. Nv
echo "[BLAST] Ac vs. Nv started at `date`"
bash map_genes.sh --tr1 ../data/DB/Ac.fa \
                  --t1 prot \
                  --n1 Ac \
                  --tr2 ../data/DB/Nv.fa \
                  --t2 prot \
                  --n2 Nv \
# Hv vs. Nv
echo "[BLAST] Hv vs. Nv started at `date`"
bash map_genes.sh --tr1 ../data/DB/Hv.fa \
                  --t1 prot \
                  --n1 Hv \
                  --tr2 ../data/DB/Nv.fa \
                  --t2 prot \
                  --n2 Nv \

# move directories


# Job finished
echo "Run NCBI BLAST finished at `date`"
