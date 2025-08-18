#!/bin/bash
## mmseqs2_clusterupdate.sh
# This script update new protemoes in DB for MMseqs2 clustering technique
# using greedy lincluster. Add new proteomes FASTA files in DB/
# It uses 8 AMD CPU-cores. Can only mbe used for greedy clustering
# Submit this script in the project directory!

#SBATCH --job-name=mmseqs2_update       # Job name    (default: sbatch)
#SBATCH --output=mmseqs2_update-%j.out  # Output file (default: slurm-%j.out)
#SBATCH --error=mmseqs2_update-%j.err   # Error file  (default: slurm-%j.out)
#SBATCH --constraint=amd                # Target AMD CPU (7452 or 9554)
#SBATCH --partition=basic               # Unlimited time (no hpc project)
#SBATCH --ntasks=1                      # Number of tasks
#SBATCH --cpus-per-task=8               # Number of CPUs per task
#SBATCH --mem-per-cpu=1GB               # Memory per CPU (in MB)
#SBATCH --time=3-00:00:00               # Wall clock time limit (H:M:S)
#SBATCH --constraint=amd                # Target AMD CPU (7452 or 9554)
#SBATCH --partition=basic               # Unlimited time (no hpc project)

# Inform the user
echo "MMseqs2 cluster update started"

# load module and print information 
module load mmseqs2/17-b804f    # load mmseqs2 and run completion
source /lisc/app/mmseqs2/17-b804f/util/bash-completion.sh
module list                     # list the loaded modules
lscpu | grep "Model name"       # CPU Architecture
date                            # Time of submitted job

#-----------------------------------------------------#
# PROJECT DEPENDANT VARIABLES AND CONSTANT (ADJUST!)  #
#-----------------------------------------------------#
species=("Ac" "Hv" "Nv")    # .fa files in DB

#-----------------------------------------------------#
# Execute MMseqs2 pipeline                            #
#-----------------------------------------------------#
## Step 1: Create MMseqs2 data base
# Create database in as one single .fa file (temporary)
# Convert the single .fa file into mmseq2 database structure
cd data/DB
# Create empty DB.fa file (temporary .fa file)
true > DB.fa
output="DB.fa"

# Loop over each species
for sp in "${species[@]}"; do
  infile="${sp}.fa"
  if [[ -f "$infile" ]]; then
    cat "$infile" >> "$output"
  else
    echo "Warning: File $infile not found, skipping." >&2
  fi
done

# convert into mmseq2 data structure
mmseqs createdb DB.fa mmseqs2/Greedy_clust/DB_new
rm -r DB.fa

## Step 2: Cluster update
# Update DB (old) with DB_new and get DB_new_update (new). The DB_clu gets updated as
# DB_update_clu. The older version get deleted and the new version get renamed.
cd mmseqs2/Greedy_clust
mkdir tmp
mmseqs clusterupdate DB DB_new DB_clu DB_new_updated DB_update_clu tmp
# rmdir tmp DB DB_new DB_clu
# mv DB_new_updated DB
# mv DB_update_clu DB_clu
mmseqs createtsv DB_new DB_new DB_update_clu DB_update_clu.tsv

# Inform the user with time
echo "MMseqs2 linear cluster finished at `date`"

