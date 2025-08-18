#!/bin/bash
# ------------------------------------------------------------------------------
## 04_preprocess_orthology.sh
# Preprocess OMA and OrthoFinder results. Create consistent input for the fully
# automatized OrthoMAP workflow.
# ------------------------------------------------------------------------------

# Load modules
module load R

# Run pre-processing
Rscript src/data_preprocessing/create_oma_table.R results/OMA
Rscript src/data_preprocessing/create_orthofinder_table.R results/OrthoFinder/diamond/Phylogenetic_Hierarchical_Orthogroups
Rscript src/data_preprocessing/create_orthofinder_table.R results/OrthoFinder/mmseqs2/Phylogenetic_Hierarchical_Orthogroups
# Rscript src/data_preprocessing/create_orthofinder_table.R results/OrthoFinder/blast/Phylogenetic_Hierarchical_Orthogroups