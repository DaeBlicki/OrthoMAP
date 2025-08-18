#!/bin/bash
## mmseqs2_linclust_BLAST.sh
# This is script to analyse protemoes with MMseqs2 clustering technique
# in O(n). This enhance the runtime in cost with sensitivity. The output
# is a .tsv file and cluster folder. It uses 12 AMD CPU-cores.
# Submit this script in the project directory!

#SBATCH --job-name=mmseqs2_BLAST      # Job name    (default: sbatch)
#SBATCH --output=mmseqs2_BLAST-%j.out # Output file (default: slurm-%j.out)
#SBATCH --error=mmseqs2_BLAST-%j.err  # Error file  (default: slurm-%j.out)
#SBATCH --constraint=amd              # Target AMD CPU (7452 or 9554)
#SBATCH --partition=basic             # Unlimited time (no hpc project)
#SBATCH --ntasks=1                    # Number of tasks
#SBATCH --cpus-per-task=64            # Number of CPUs per task
#SBATCH --mem-per-cpu=1GB             # Memory per CPU (in MB)
#SBATCH --time=72:00:00               # Wall clock time limit (H:M:S)
#SBATCH --constraint=amd              # Target AMD CPU (7452 or 9554)
#SBATCH --partition=basic             # Unlimited time (no hpc project)

# Inform the user
echo "MMseqs2 linear cluster started"
echo "[MMseqs2 clustering]: BLAST"
echo "[PARAMETER]: --min-seq-id 0.9"
echo "[PARAMETER]: --cov-mode 1"
echo "[PARAMETER]: -c 0.8"
echo "[PARAMETER]: --cluster-mode 1"

# load module and print information 
module load mmseqs2/17-b804f    # load mmseqs2 and run completion
source /lisc/app/mmseqs2/17-b804f/util/bash-completion.sh
module list                     # list the loaded modules
lscpu | grep "Model name"       # CPU Architecture
date                            # Time of submitted job

#-----------------------------------------------------#
# PROJECT DEPENDANT VARIABLES AND CONSTANT (ADJUST!)  #
#-----------------------------------------------------#
species=("ta" "ed" "aq" "em" "sd" "sl" "xs" "ad" "am" "ep" "tripc" "ch"\
         "hm" "hv" "re" "ac" "nv2" "sc" "sm" "ce" "cr" "dm" "tc" \
         "py" "sb" "bf" "lo" "hs" "mm" "sa" "re2" "kostya")

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
mkdir -p ../../results/mmseqs2/BLAST_clust/
mmseqs createdb DB.fa ../../results/mmseqs2/BLAST_clust/DB
rm -r DB.fa

## Step 2: Cluster with MPI and result output
# Linclust is a clustering in linear time. It is magnitudes faster 
# but a bit less sensitive than clustering. Output as .tsv file.
#     parameter list:
#         --min-seq-id FLOAT  : 
#         --cov-mode 1        : coverage of target
#         -c FLOAT            : aligned (covered) residues
#         --cluster-mode 1    : BLASTclust
#          
cd ../../results/mmseqs2/BLAST_clust
mkdir tmp
mmseqs linclust DB DB_clu tmp \
              --min-seq-id 0.9\
              --cov-mode 1 -c 0.8\
              --cluster-mode 1\
              --threads 64
rm -r tmp
mmseqs createtsv DB DB DB_clu DB_clu.tsv

# Inform the user with time
echo "MMseqs2 linear cluster finished at `date`"

