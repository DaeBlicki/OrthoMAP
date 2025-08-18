#!/bin/bash
## orthofinder_update.sh
# This is script update the previous calculated data base with new species.
# Submit this script in the project directory!

#SBATCH --job-name=orthofinder        # Job name    (default: sbatch)
#SBATCH --output=orthofinder-%j.out   # Output file (default: slurm-%j.out)
#SBATCH --error=orthofinder-%j.err    # Error file  (default: slurm-%j.out)
#SBATCH --constraint=amd              # Target AMD CPU (7452 or 9554)
#SBATCH --partition=basic             # Unlimited time (no hpc project)
#SBATCH --ntasks=1                    # Number of tasks
#SBATCH --cpus-per-task=64            # Number of CPUs per task
#SBATCH --mem-per-cpu=1GB             # Memory per CPU (in MB)
#SBATCH --time=2-00:00:00             # Wall clock time limit (H:M:S)
#SBATCH --constraint=amd              # Target AMD CPU (7452 or 9554)
#SBATCH --partition=basic             # Unlimited time (no hpc project)

# Inform the user
echo "OrthoFinder update"

# load module and print information 
module load orthofinder         # load orthofinder
module list                     # list the loaded modules
lscpu | grep "Model name"       # CPU Architecture
date                            # Time of submitted job

# update orthofinder with 64 CPU cores
orthofinder.py -b results/OrthoFinder/diamond -f data/DB -a 64

# Inform the user with time
echo "Orthofinder update finished at `date`"
