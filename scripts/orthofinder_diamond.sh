#!/bin/bash
## orthofinder.sh
# This is script to analyse proteomes with orthofinder using 64 AMD CPU-cores
# Submit this script in the project directory!

#SBATCH --job-name=orthofinder_diamond          # Job name    (default: sbatch)
#SBATCH --output=orthofinder_diamond-%j.out     # Output file (default: slurm-%j.out)
#SBATCH --error=orthofinder_diamond-%j.err       # Error file  (default: slurm-%j.out)
#SBATCH --constraint=amd                        # Target AMD CPU (7452 or 9554)
#SBATCH --partition=basic                       # Unlimited time (no hpc project)
#SBATCH --ntasks=1                              # Number of tasks
#SBATCH --cpus-per-task=64                      # Number of CPUs per task
#SBATCH --mem-per-cpu=1GB                       # Memory per CPU (in MB)
#SBATCH --time=4-00:00:00                       # Wall clock time limit (H:M:S)
#SBATCH --constraint=amd                        # Target AMD CPU (7452 or 9554)
#SBATCH --partition=basic                       # Unlimited time (no hpc project)

# Inform the user
echo "OrthoFinder (Diamond) pipeline"

# load module and print information 
module load orthofinder         # load orthofinder
module list                     # list the loaded modules
lscpu | grep "Model name"       # CPU Architecture
date                            # Time of submitted job

## run orthofinder with 64 CPU cores
orthofinder.py -f 'data/DB' -o 'results/OrthoFinder/diamond' -t 64 -a 64 -S diamond
res_dir=$(find results/OrthoFinder/diamond/ -mindepth 1 -maxdepth 1 -type d)

# Sanity check
if [ -d "$res_dir" ]; then
    echo "Moving results from: $res_dir"
    mv "$res_dir"/* results/OrthoFinder/diamond/
    rmdir "$res_dir"
else
    echo "❌ OrthoFinder run failed or no result directory found."
fi

# Inform the user with time
echo "Orthofinder pipeline finished at `date`"
