#!/bin/bash
## oma_cheatsheet.sh
# This cheatsheet explains the progress of the pipeline. The pipeline contains
# three part: database conversion part, all-against-all part, and orthology 
# interference part.
# The second part can be parallelized and hence, this pipeline execute the 
# second part in parallel. However, parallel computing can cause some processors
# to crash and be killed. Therefore, below is some troubleshooting.

## Run Orthologuos MAtrix (OMA) pipeline
# Run this code in project repository, ensure the job is finished before submit
# the next part. Status can be checked at any point of the pipeline.
sbatch < scripts/oma_part1.sh # Database Conversion (takes seconds to minutes)
sbatch < scripts/oma_part2.sh # All-Against All (takes hours to days)
sbatch < scripts/oma_part3.sh # Orthology Interference

## Troubleshooting: oma-cleanup and oma-status
# Ensure to be in the data repository (must contain DB/ and Cashe/ subdirectory).
module load omastandalone/2.6.0

## oma-cleanup: Removes unfinished all-against-all alignements in oma part 2
# Goal: Clean-up OMA part 2 in case of interuption (Cluster time out)
# Require: OMA analysis stopped and go data/ subdirectory, need `parameter.drw`
#       1) Ensure no OMA job is running -> [y]
#       2) Ensure no OMA process is running -> [y]
/lisc/app/omastandalone/2.6.0/OMA/OMA.2.6.0/bin/oma-cleanup

## oma-status: Checks the status of the OMA pipeline
# Goal: Debugging the process and estimate required tiem
# Require: Go in data/ subdirectory, need `parameter.drw`
# example output:
#   Summary of OMA standalone All-vs-All computations:
#   --------------------------------------------------
#   Nr chunks started: 50 (4.42%)
#   Nr chunks finished: 325 (28.71%)
#   Nr chunks finished w/o exported genomes: 325 (28.71%)
/lisc/app/omastandalone/2.6.0/OMA/OMA.2.6.0/bin/oma-status  # status
