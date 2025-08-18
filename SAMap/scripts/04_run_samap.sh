#!/bin/bash
# ------------------------------------------------------------------------------
## 04_run_samap.sh
# Run SAMap workflow
# ------------------------------------------------------------------------------
#SBATCH --job-name=samap                # Job name    (default: sbatch)
#SBATCH --output=samap-%j.out           # Output file (default: slurm-%j.out)
#SBATCH --error=samap-%j.err            # Error file  (default: slurm-%j.out)
#SBATCH --constraint=amd                # Target AMD CPU (7452 or 9554)
#SBATCH --partition=basic               # Unlimited time (no hpc project)
#SBATCH --ntasks=1                      # Number of tasks
#SBATCH --cpus-per-task=1               # Number of CPUs per task
#SBATCH --mem-per-cpu=130GB             # Memory per CPU (in MB)
#SBATCH --time=24:00:00                 # Wall clock time limit (H:M:S)
#SBATCH --constraint=amd                # Target AMD CPU (7452 or 9554)
#SBATCH --partition=basic               # Unlimited time (no hpc project)

# load module and print information 
module load python3 conda       # list the loaded modules
module list                     # list the loaded modules
lscpu | grep "Model name"       # CPU Architecture
date                            # Time of the submitted job

# Run SAMap
conda activate SAMap
python3 -u src/SAMap_main.py

echo "SAMap finished at `date`"