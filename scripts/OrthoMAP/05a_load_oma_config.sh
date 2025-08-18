#!/bin/bash
# ------------------------------------------------------------------------------
## 05a_load_oma_config.sh
# Script to create symbolic link to the configurations of the project. This
# script loads (forced!) the configuration with the OMA standalone inputs mapped
# to OMA groups.
#
# Load the configuration before run `05_run_orthomap.sh` 
# ------------------------------------------------------------------------------
cd data
rm config.yaml
ln -s configuration/05a_config_oma.yaml config.yaml
cd ..