# Project Journal

This **Markdown** have the purpose to replicate the results of the Thesis and provide a step-by-step guide.

## 1. Download and pre-process single-cell RNA sequencing atlas

### Acces single-cell atlas and FASTA files

Please contact Dr. Cole A.G (alison.cole@univie.ac.at) or Dr. Montenegro J. (juan.montenegro@univie.ac.at) for data access for *Nematostella vectensis* and *Aurelia coerulea*. Ask for the Annotation table `.csv`, Fasta file `.fa`, and single-cell atlas `.Robj` or `.rds`. Place the annotiation table in `data/annotation_table`, the `.fa` in `data/DB/` and the `.Robj` in `data/single_cell_atlas`. After the data access was sucessful, run following in the project root:

```
source scripts/project/01_download.sh
source scripts/project/02_preprocess_data.sh
```

1) It downloads and pre-processs the scRNA-seq data from *Hydra vulgaris*. The script will download [`aepAtlasDubInclude.rds`](https://pubmed.ncbi.nlm.nih.gov/36639202/) and the list of the non-doublets. In the next step, possible doublets are getting removed (see [05_hydraAtlasReMap.md](https://github.com/cejuliano/brown_hydra_genomes/blob/main/05_hydraAtlasReMap.md)) as described in their paper. Additionally, the annotation table is created to match `.fa` and `.Robj`. Last but not least, it creates a new meta.data `IDs` by removing the prefix: `I_`, `Ec_`, and `En_` in `CuratedIdent` and update the Seurat Object to the newest version. 

2) It removes the isoforms in Aurelias protein model by choosing the longest sequence per geneID. 

## 2. Run orthology interference analysis
The orthology interference was analyzed with OMA standalone and OrthoFinder (Diamond and MMseqs2 algorithm). Run ```source scripts/project/03_run_orthology.sh``` in the project repository. This can take a while (2-10h, for each analyis using 64 AMD cores). **After** the results are produced, run ```source scripts/project/04_preprocess_orthology.sh``` to evaluate the orthologs table used in the OrthoMAP workflow.

## 3. Run OrthoMAP
In `data/configuration/` are all `config.yaml` stored used in the project. Below is an overview:

- `05a` : OMA analysis on OMA groups
- `05b` : OrthoFinder MMseqs2 on HOG (Hierarchical Orthologous Groups)
- `05c` : OrthoFinder MMseqs2 on OOG (1:1 Orthologs)
- `05d` : OrthoFinder Diamond on HOG
- `05e` : OrthoFinder Diamond on OOG

For example, if you want to load OMA configuration, run following command in project root ```source scripts/project/05a_load_oma_config.sh```. After the configuration is loaded, run ```sbatch < scripts/project/05_run_orthomap.sh``` in SLURM or ````source scripts/project/05_run_orthomap.sh``` on the local device.