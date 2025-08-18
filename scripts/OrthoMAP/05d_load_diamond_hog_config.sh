#!/bin/bash
# ------------------------------------------------------------------------------
## 05c_load_diamond_hog_config.sh
# Script to create symbolic link to the configurations of the project. This
# script loads (forced!) the configuration with the OrthoFinder using Diamond.
# It maps to Hierarchical Orthologs Groups. 
# 
# Load the configuration before run `05_run_orthomap.sh` 
# ------------------------------------------------------------------------------
cd data
rm config.yaml
ln -s configuration/05d_config_diamond_hog.yaml config.yaml
cd ..