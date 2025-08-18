#!/bin/bash
# oma_part3.sh
# This is the orthology interference part

#SBATCH --job-name=oma3               # Job name    (default: sbatch)
#SBATCH --output=oma3-%j.out          # Output file (default: slurm-%j.out)
#SBATCH --error=oma3-%j.err           # Error file  (default: slurm-%j.out)
#SBATCH --constraint=amd              # Target AMD CPU (7452 or 9554)
#SBATCH --partition=basic             # Unlimited time (no hpc project)
#SBATCH --ntasks=1                    # Number of tasks
#SBATCH --cpus-per-task=1             # Number of CPUs per task
#SBATCH --mem-per-cpu=30GB             # Memory per CPU
#SBATCH --time=7-00:00:00             # Wall clock time limit (H:M:S)
#SBATCH --constraint=amd              # Target AMD CPU (7452 or 9554)
#SBATCH --partition=basic             # Unlimited time (no hpc project)

# Inform the user
echo "[STEP 3]: Orthology Interference Part"

# load module and print information 
module load omastandalone/2.6.0 # load OMA (standalone version)
module list                     # list the loaded modules
lscpu | grep "Model name"       # CPU Architecture
date                            # Time of the submitted job

# run oma
cd data
OMA

# Inform the user with time
echo "[STEP 3]: Process finished at `date`"

# move directory to results
mv Cache ../results/OMA
