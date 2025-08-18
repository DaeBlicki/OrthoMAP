#!/bin/bash
# ------------------------------------------------------------------------------
## 01_download.sh
# Download the protein models `.fa` and single-cell data `.Robj` for the 
# published Hydra vulgaris [1]. For data access for Nematostella vectensis and
# Aurelia coerulea, please message my supervisors below. Ask for the FASTA,
# .Robj, and Annotation table as .csv.
#
# Dr. Cole A.G (alison.cole@univie.ac.at) 
# Dr. Montenegro J. (juan.montenegro@univie.ac.at)
#
# [1] J. F. Cazet et al. “A chromosome-scale epigenetic map of the Hydragenome 
# reveals conserved regulators of cell state”. In: Genome Research 33 (2023), 
# pp. 283–298. DOI: https://doi.org/10.1101/gr.277772.122
# ------------------------------------------------------------------------------
mkdir download
cd download

# Download Hydra vulgaris [1]
wget -O HVAEP.Robj https://research.nhgri.nih.gov/HydraAEP/download/scriptsdata/aepAtlasDubInclude.rds
wget -O HVAEP.csv https://raw.githubusercontent.com/cejuliano/brown_hydra_genomes/refs/heads/main/ID_Conversion/HVAEP.master.conversion.table.csv
wget -O HVAEP_nondubs.tsv https://raw.githubusercontent.com/cejuliano/brown_hydra_genomes/refs/heads/main/05_hydraAtlasReMap/04_finalize/nondub.tsv
wget -O Hv.fa.gz https://research.nhgri.nih.gov/HydraAEP/download/sequences/hv_aep/HVAEP.prot.fa.gz
gunzip Hv.fa.gz
mv Hv.fa ../data/DB
cd ../