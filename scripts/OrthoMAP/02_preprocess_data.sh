#!/bin/bash
# ------------------------------------------------------------------------------
## 02_preprocess_data.sh
# Preprocess the protein models `.fa` and single-cell data `.Robj` for the 
# published Hydra vulgaris [1]. It removes doublets using `HVAEP_nondubs.tsv`
# and remove cluster 37 and 41 as described in their paper.
#
# [1] J. F. Cazet et al. “A chromosome-scale epigenetic map of the Hydragenome 
# reveals conserved regulators of cell state”. In: Genome Research 33 (2023), 
# pp. 283–298. DOI: https://doi.org/10.1101/gr.277772.122
# ------------------------------------------------------------------------------

# Load modules
module load seurat5
module load bioawk

## Step 1: Run preprocessing hydra data
# 1) Removes doublets using the protocol from their paper
# 2) Create another column `IDs` in the meta.data by removing prefix 
#    in `CuratedIdent`.
# 3) Create the annotation table for hydra
Rscript src/data_preprocessing/create_Hv_seurat_object.R
Rscript src/data_preprocessing/create_Hv_csv.R
Rscript src/data_preprocessing/create_tissue_coloring.R

## Step 2: Run preprocess Aurelia fasta file
# 1) generate table: geneID | isoform | sequence length
bioawk -c fastx '{ split($name, a, ".t"); print a[1], $name, length($seq) }' data/DB/Ac.fa  > download/Ac_isoform_table.txt
# 2) select the longest isoform per gene
awk '
{
  if ($1 in max) {
    if ($3 > max[$1]) {
      max[$1] = $3;
      best[$1] = $2;
    }
  } else {
    max[$1] = $3;
    best[$1] = $2;
  }
}
END {
  for (g in best) print best[g];
}
' download/Ac_isoform_table.txt > download/Ac_selected_isoform.txt
# 3) create fasta file by keeping the selected isoforms
bioawk -c fastx 'BEGIN {
  while ((getline < "download/Ac_selected_isoform.txt") > 0) keep[$1] = 1
}
keep[$name] {
  gene = $name
  sub(/\.t[0-9]+$/, "", gene)  # remove .t%d from end
  print ">" gene
  for (i = 1; i <= length($seq); i += 100)
    print substr($seq, i, 100)
}' data/DB/Ac.fa > data/DB/tmp.fa
rm -r data/DB/Ac.fa
mv data/DB/tmp.fa data/DB/Ac.fa