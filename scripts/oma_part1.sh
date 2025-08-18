#!/bin/bash
## oma_part1.sh
# This is the database conversion part

#SBATCH --job-name=oma1               # Job name    (default: sbatch)
#SBATCH --output=oma1-%j.out          # Output file (default: slurm-%j.out)
#SBATCH --error=oma1-%j.err           # Error file  (default: slurm-%j.out)
#SBATCH --constraint=amd              # Target AMD CPU (7452 or 9554)
#SBATCH --partition=basic             # Unlimited time (no hpc project)
#SBATCH --ntasks=1                    # Number of tasks
#SBATCH --cpus-per-task=1             # Number of CPUs per task
#SBATCH --mem-per-cpu=1GB             # Memory per CPU (in MB)
#SBATCH --time=02:00:00               # Wall clock time limit (H:M:S)
#SBATCH --constraint=amd              # Target AMD CPU (7452 or 9554)
#SBATCH --partition=basic             # Unlimited time (no hpc project)

# Inform the user
echo "[STEP 1]: Database Conversion Part"

# load module and print information 
module load omastandalone/2.6.0 # load OMA (standalone version)
module list                     # list the loaded modules
lscpu | grep "Model name"       # CPU Architecture
date                            # Time of the submitted job


# run OMA part 1
cd data 
OMA -c

# Inform the user with time
echo "[STEP 1]: Process finished at `date`"
