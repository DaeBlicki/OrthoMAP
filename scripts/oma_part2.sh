#!/bin/bash
# oma_part2.sh
# This is the all-against-all part, it is split into 64 parallel jobs

#SBATCH --job-name=oma2               # Job name    (default: sbatch)
#SBATCH --output=oma2-%j.out          # Output file (default: slurm-%j.out)
#SBATCH --error=oma2-%j.err           # Error file  (default: slurm-%j.out)
#SBATCH --constraint=amd              # Target AMD CPU (7452 or 9554)
#SBATCH --partition=basic             # Unlimited time (no hpc project)
#SBATCH --ntasks=1                    # Number of tasks
#SBATCH --cpus-per-task=32            # Number of CPU cores per task
#SBATCH --mem-per-cpu=2GB             # Memory per CPU
#SBATCH --time=7-00:00:00             # Wall clock time limit (H:M:S)
#SBATCH --constraint=amd              # Target AMD CPU (7452 or 9554)
#SBATCH --partition=basic             # Unlimited time (no hpc project)

# Inform the user
echo "[STEP 2]: All-Against-All Part"

# load module and print information 
module load omastandalone/2.6.0 # load OMA (standalone version)
module list                     # list the loaded modules
lscpu | grep "Model name"       # CPU Architecture
date                            # Time of the submitted job

# run OMA part 2
cd data
export NR_PROCESSES=32
OMA -s -n $NR_PROCESSES

# Inform the user with time
echo "[STEP 2]: Process finished at `date`"
