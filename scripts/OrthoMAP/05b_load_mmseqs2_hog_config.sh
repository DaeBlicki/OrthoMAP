#!/bin/bash
# ------------------------------------------------------------------------------
## 05c_load_diamond_hog_config.sh
# Script to create symbolic link to the configurations of the project. This
# script loads (forced!) the configuration with the OrthoFinder using MMseqs2.
# It maps to Hierarchical Orthologs Groups. 
# 
# Load the configuration before run `05_run_orthomap.sh` 
# ------------------------------------------------------------------------------
cd data
rm config.yaml
ln -s configuration/05b_config_mmseqs2_hog.yaml config.yaml
cd ..