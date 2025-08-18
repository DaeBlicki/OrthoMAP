#!/bin/bash
# ------------------------------------------------------------------------------
## 03_run_orthology.sh
# Run orthology interference tools OMA [1] and OrthoFinder [2]. OMA is used
# to identify OMA groups (also called Orthogroups), OrthoFinder is used to 
# determine Hierarchical Orthologous Groups and 1:1 orthologs. Furthermore, the
# OrthoFinder can apply Diamond [3], MMseqs2 [4] or Blast [5] algorithm.
#
# [1] Adrian Altenhoff et al. “OMA standalone: Orthology inference among
# public and custom genomes and transcriptomes”. In: Genome Research 29 
# (June 2019), pp. 1152–1163. DOI: https://doi.org/10.1101/gr.243212.118.
#
# [2] David Emms and Steven Kelly. “OrthoFinder: Phylogenetic orthology
# inference for comparative genomics”. In: Genome Biology 20 (Nov. 2019),
# p. 238. DOI: https://doi.org/10.1186/s13059-019-1832-y.
#
# [3] Benjamin Buchfink, Chao Xie, and Daniel H. Huson. “Fast and sensitive
# protein alignment using DIAMOND”. In: Nature Methods 12.1 (2015),
# pp. 59–60. DOI: https://doi.org/10.1038/nmeth.3176.
#
# [4] Martin Steinegger and Johannes S ¨oding. “MMseqs2 enables sensitive
# protein sequence searching for the analysis of massive data sets”. In:
# Nature Biotechnology 35 (Oct. 2017), pp. 1026–1028. 
# DOI: https://doi.org/10.1038/nbt.3988
#
# [5] Camacho, C., Coulouris, G., Avagyan, V. et al. BLAST+: architecture and 
# applications. BMC Bioinformatics 10, 421 (2009). 
# DOI: https://doi.org/10.1186/1471-2105-10-421
# ------------------------------------------------------------------------------

# Run OMA standalone
jid1=$(sbatch --parsable scripts/oma_part1.sh) &&
jid2=$(sbatch --parsable --dependency=afterok:$jid1 scripts/oma_part2.sh) &&
sbatch --dependency=afterok:$jid2 scripts/oma_part3.sh

## Run OrthoFinder (diamond, mmseqs2 or blast)
sbatch < scripts/orthofinder_diamond.sh
sbatch < scripts/orthofinder_mmseqs2.sh
# sbatch < scripts/orthofinder_blast.sh