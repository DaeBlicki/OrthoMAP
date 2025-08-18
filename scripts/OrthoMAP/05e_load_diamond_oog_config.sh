#!/bin/bash
# ------------------------------------------------------------------------------
## 05b_load_diamond_oog_config.sh
# Script to create symbolic link to the configurations of the project. This
# script loads (forced!) the configuration with the OrthoFinder using Diamond.
# It maps to 1:1 orthologs. 
# 
# Load the configuration before run `05_run_orthomap.sh` 
# ------------------------------------------------------------------------------
cd data
rm config.yaml
ln -s configuration/05e_config_diamond_oog.yaml config.yaml
cd ..