mkdir tmp
mv Citation.txt Comparative_Genomics_Statistics/ Gene_Duplication_Events/ Gene_Trees/ Log.txt MultipleSequenceAlignments/ Orthogroups Orthogroup_Sequences/ Orthologues/ Phylogenetically_Misplaced_Genes/ Putative_Xenologs/ Resolved_Gene_Trees/ Single_Copy_Orthologue_Sequences/ Species_Tree/ WorkingDirectory/ tmp
tar -czvf tmp.tar.gz tmp
rm -r tmp