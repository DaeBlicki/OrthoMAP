# Guideline

## 1 Environemnt and Requirement

Computational results of this work have been achieved using the *Life Science Compute Cluster* (LiSC) of the University of Vienna. Make sure to have following softwares installed on your device:

## 2 Run Orthology Interference Analysis

Place protein models in questions as `.fa` in the `data/DB` subdirectory. Next, make sure to adapt `parameters.drw` when using OMA standalone. Finally, run following command in console:

```
## Run OMA analysis
# Make sure to adapt 'parameters.drw' in the results/ subdirectory
jid1=$(sbatch --parsable scripts/oma_part1.sh) &&
jid2=$(sbatch --parsable --dependency=afterok:$jid1 scripts/oma_part2.sh) &&
sbatch --dependency=afterok:$jid2 scripts/oma_part3.sh

## Run OrthoFinder (diamond, mmseqs2 or blast) analysis
sbatch < scripts/orthofinder_diamond.sh
sbatch < scripts/orthofinder_mmseqs2.sh
sbatch < scripts/orthofinder_blast.sh
```

## 3 Run OrthoMAP

## Getting started

The `data/` subdirectory is structured as follows:

- In `data/scsRNA_metadata.csv`: Metadata used in the pipeline
- In `data/DB`: Transcriptome objects as `.fa`  
- In `data/single_cell_atlas`: Single-cell data as `.Robj`
- In `data/annotation_table` : Maps for gene.id (`.fa`) to gene.name (`.Robj`)
- In `data/config.yaml`: Contains configuration for the resulting path in the pipeline and which input the user wants to use

### scsRNA_metadata.csv
This file contains **mandatory** and **important** information used in the pipeline. Here is a description for the colnames with an example.

- @**file_id**: Name of the file for the species, example: "Nv", "nv"
- @**scientific_name**: Full species name, example: "Nematostella vectensis"
- @**orig.ident**: Source sample or batch, example:"orig.ident"
- @**ID.separate**: Cell type identifier, example: "ID.separate"
- @**sc_version**: Single-cell (seurat) data version, example: "5.0.1"
    - seurat objects can be loaded using `get(load(<path_to_seurat_object>))`
    - get version using `Version(<your_seurat_object>)`
    - pipeline automatically calls `UpdateSeuratObject`

⚠️ **Trouble-shooting: Error loading seurat objects**

- Try using `readRDS(<path_to_seurat_object>)`
- The single-cell data may be outdated and thus, these files are labelled as "old" in sc_version
- pipeline automatically handels this case, but it needs **user information stored in scsRNA_metadata.csv**.
- However, the pipeline load the seurat file multiple time and when it uses an "old" seurat object, the seurat file will updated everytime when is getting loaded. To improve the runtime, please use updated versions of seurat (5.0.1).

### annotation_table.csv
This file must need **geneID** (gene name in `.fa`) and **gene.name** (gene name in `.Robj`). In case that the proteome transcript in FASTA files contains multiple transcript for one gene (usually marked at the end with `.t1`,`.t2`, ...) and the annotation table does not distinct the isoforms, the pipeline invalidate the projection from gene to orthogroup when multiple isoforms map to different orthogroups. Read the README.md for the discussion.

### config.yaml
This file contains information which input is used and how the results are stored. Make sure to remain the configuration structure.

## Run OrthoMAP

OrthoMAP is a Rscript to create a Seurat Object with an orthologous gene cluster (orthogroups) against cells expression. In the project repository, run following command: ```Rscript src/main.R (opt)[f]```

**[f] Pipeline flag**
- **x** : execute everything, all-in-one (full pipeline)
- **s** : start, run only step 1 (create OrthoMAP Seurat Object)
- **c** : continue, run step 2 (standard Seurat standard workflow)
- **p** : peek, run visualization of pre-processed OrthoMAP object
- **v** : visualize, run visualization of processed OrthoMAP result