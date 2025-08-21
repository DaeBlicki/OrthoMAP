# OrthoMAP 

This project aims to implement a pipeline for create an embedding for cross-species comparisons. It depends on orthologous genes identified by any orthology interference software and data from scRNA-seq.

## Abstract

## Experimental Results
![alt text](https://github.com/DaeBlicki/OrthoMAP/blob/main/results.png)

## 1 Guideline

### 1.1 Download Project Repository

Run the following commands in your console. The command will download the project code.

```
# clone the git repository
git clone https://github.com/DaeBlicki/OrthoMAP.git
cd OrthoMAP

# create `data/` structure
mkdir data
mkdir data/DB data/single_cell_atlas data/annotation_table
```

The project structure after running the pipeline looks like this.

```
OrthoMAP
├── data
│    ├── annotation_table/    * store annotation table
│    ├── DB/                  * store FASTA files
│    ├── sc_objects/          * store Seurat Objects
│    ├── config.yaml          * User configuration file
│    ├── parameter.drw        * OMA parameter file
│    ├── scsRNA_metadata.csv  * User dataset file
│    └── tissue_palette.csv   * (optional) color palette
│ 
├── doc
│    ├── Guideline.md         * Tutorial for OrthoMAP
│    └── Project_Journal.md   * Description for reproducibility 
│
├── results
│    ├── OMA/                           * results for OMA
│    ├── OrthoFinder/diamond            * results for Diamond
│    ├── OrthoFinder/mmseqs2            * results for MMseqs2
│    ├── OrthoMAP_Data_Visualization/   * visual results
│    ├── OrthoMAP_Seurat_Objects/       * seurat results
│    └── OrthoMAP_Statistical_Results/  * empirical results
│ 
├─── scripts
│    ├── oma_part1.sh              * Part 1 of OMA
│    ├── oma_part2.sh              * Part 2 of OMA
│    ├── oma_part3.sh              * Part 3 of OMA
│    ├── orthofinder_diamond.sh    * run OrthoFinder using Diamond
│    ├── orthofinder_mmseqs2.sh    * run OrthoFinder using MMseqs2
│    └── run_orthomap.sh           * run OrthoMAP
│
└── logs.tar.gz     * Project output and error files
```

### 1.2. Requirement

The pipeline requires **three** inputs for each species: Transcriptomes stored as `.fa` in `data/DB`, scRNA-seq data stored as `.Robj` or `.rds` in `data/sc_objects`, and the annotation table from FASTA file (geneID) to Robj file (gene.name) stored as `.csv` in `data/annotation_table`. The required software and R packages are shown in Supplementary Material - Software. (Optional) Color palette for coloring celltypes stored in `IDs` can be integrated in the pipeline using `data/tissue_palette.csv`

**⚠️ IMPORTANT: Each species inputs must have the same filename**

## 2. OrthoMAP Workflow

### STEP 1: Identify comparable genes between species

This part was exectued with [OMA standalone (Orthologous MAtrix)](https://doi.org/10.1101/gr.243212.118)[1], [OrthoFinder](https://doi.org/10.1186/s13059-019-1832-y)[2,3], using [MMseqs2](https://doi.org/10.1038/nbt.3988)[4] and [Diamond](https://doi.org/10.1038/nbt.3988)[5]. In the project repository, run following command.

⚠️ Make sure to adapt **'parameters.drw'** in the `data/` subdirectory for the OMA run.

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

In case you want to update the data analysis using more species, add `.fa` file in `data/DB/` subfolder (symbolic links works aswell). Make sure to update 'parameter.drw' in results/ and $species in scripts/orthofinder_update.sh. Afterwards, run following command line in the project repository.

```
# Update OMA 
jid1=$(sbatch --parsable scripts/oma_part1.sh) &&
jid2=$(sbatch --parsable --dependency=afterok:$jid1 scripts/oma_part2.sh) &&
sbatch --dependency=afterok:$jid2 scripts/oma_part3.sh

# Update OrthoFinder
sbatch < scripts/orthofinder_update.sh
```

### STEP 2: Produce orthologous gene embeddings
This was developed in [*R* (v.4.5.0)](https://www.R-project.org/)[6]. The workflow create a [Seurat](https://doi.org/10.1016/j.cell.2021.04.048) [7] Object with an orthologous gene cluster (orthogroups) against cells expression. See Supplementary Material for the used R packages. In the project repository, run following command: ```Rscript src/main.R (opt)[f]```

**[f] Pipeline flag**
- **x** : execute everything, all-in-one (full pipeline)
- **s** : start, run only step 1 (create OrthoMAP Seurat Object)
- **c** : continue, run step 2 (standard Seurat standard workflow)
- **p** : peek, run visualization of pre-processed OrthoMAP object
- **v** : visualize, run visualization of processed OrthoMAP result


**⚠️ When use another data set, make sure to make changes before running the pipeline!**

- In `data/scsRNA_metadata.csv`: Metadata used in the pipeline
- In `data/DB`: Transcriptome objects as `.fa`  
- In `data/single_cell_atlas`: Single-cell data as `.Robj`
- In `data/annotation_table` : Maps for gene.id (`.fa`) to gene.name (`.Robj`)
- In `data/config.yaml`: Contains configuration for the resulting path in the pipeline and which input the user wants to use
- In `data/parameters.drw`: Parameter file for the OMA standard workflow, change $SpeciesTree when available
- In `data/parameters.drw`: (Optional) color palette for coloring the UMAP and Donut Charts and so on.

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
This file must need **geneID** (gene name in `.fa`) and **gene.name** (gene name in `.Robj`). In case that the proteome transcript in FASTA files contains multiple transcript for one gene (usually marked at the end with `.t1`,`.t2`, ...) and the annotation table does not distinct the isoforms, the pipeline invalidate the projection from gene to orthogroup when multiple isoforms map to different orthogroups. In the following, this approach is discussed:

#### ❓ How to handle gene with multiple isoforms that maps to different orthogroups?

💡 **Solution:** Map the gene to no orthogroup

✅ **Best Case Scenario:** No isoforms or all isoforms maps to the same orthogroup
- No loss of information
- No gain of missleading information

❌ **Worst Case Scenario:** High number of isoforms maps to different orthogroup
- Loss of information
- Gene that could represent an orthogroup was replaced by next best gene

🛠️ **Code Implementation**
```
Function: collapse_transcripts
Input:
    - Hashmap (isoform → orthogroup)
Output:
    - Hashmap (gene → orthogroup)
Steps:
    1. Get isoform names from Hashmap
    2. Get gene names from isoform names
    3. Map isoform to orthogroups
    4. Map gene to isoform and get orthogroups
    5. Keep gene mapping thats maps to one orthogroup
```

💬 **Discussion**

The main problem occurs when genes of a species can express multiple proteins due to various results from the splicing process. The difference between the proteins can quantitatively predicted after running the orthologous gene clustering. However, most scRNA-seq data does not differentiates them and merge the expression together as gene. When the translated proteins from the gene also map to different orthologous gene cluster, then it's not possible to predict which of the mapped orthogroups represents the gene. Furthermore, the most expressed gene represents an orthogroup. Hence, including a gene in an orthologous cluster that translates multiple proteins may expressed more than genes that only express one protein. When the multiple of the resulting proteins are also outside of the orthogroup the error to make exactly this gene the representative of the orthologous cluster may increase the uncertainty more than remove it.  

### config.yaml
This file contains information for the result path. By default, it generates the suggested structure described in the project `README`. Only change the parameters in **general** and **version** when needed. The whole part 1 of the pipeline print out the necessary information to complete **version**. This part will **NOT** be generated automatically and the user **must take full responsibility** for that (Automatize this would be an excellent practice). The parameter in **general** describes which orthogroup finder can be used as input and if verbose is activated (highly recommended).

### parameteres.drw
Relevant for the OMA analysis. Add species tree when available

## Supplementary Material

### 💾 Dataset information

| File ID | Scientific Name                | Phylum         | Class           | Notes              |
| ------- | ------------------------------ | -------------- | --------------- | ------------------ |
 |Hv      | *Hydra vulgaris*               | Cnidaria       | Hydrozoa        | (HVAEP)            |
| Nv      | *Nematostella vectensis*       | Cnidaria       | Anthozoa        | (NV2, local)       |
| Ac      | *Aurelia coerulea*             | Cnidaria       | Scyphozoa       | (local 51K)        |

## 🖥️ Hardware 

Computational results of this work have been achieved using the *Life Science Compute Cluster* (LiSC) of the University of Vienna.

- AMD EPYC™ 7452 — 32-core, *2.35 GHz*
- AMD EPYC™ 9554 — 64-core, *3.10 GHz*
- AMD EPYC™ 7662 — 64-core, *2.00 GHz*
- AMD EPYC™ 7763 — 64-core, *2.45 GHz*

## 🖥️  Software

### 🧬 Orthology interference tool and algorithms
| Software / Package          | Version    | Description                         |
|-----------------------------|------------|-------------------------------------|
| [OMA standalone](https://github.com/DessimozLab/OmaStandalone)          | 2.6.0      | Detect orthologous gene cluster     |
| [OrthoFinder](https://github.com/davidemms/OrthoFinder)             | 3.0.1b1    | Detect orthologous gene cluster     |
| [Diamond](https://github.com/bbuchfink/diamond)                 | 2.1.12     | Part of OrthoFinder application     |
| [MMseqs2](https://github.com/soedinglab/MMseqs2)                 | 17-b804f   | Part of OrthoFinder application     |


### 📦 R-Environment and packages

| Software / Package | Version      | Description                          |
|--------------------|--------------|--------------------------------------|
| [R](https://www.r-project.org/)              | 4.5.0        | Core language for data analysis      |
| [Seurat](https://satijalab.org/seurat/)         | 5.2.1        | Single-cell RNA-seq analysis toolkit |
| [tidyverse](https://www.tidyverse.org/)      | 2.0.0        | Data manipulation and visualization  |
| [yaml](https://cran.r-project.org/web/packages/yaml/index.html)           | 2.3.10       | Create project configuration         |
| [patchwork](https://cran.r-project.org/web/packages/patchwork/index.html)      | 1.3.0        | Data visualization extension         |
| [Matrix](https://cran.r-project.org/web/packages/Matrix/index.html)         | 1.7-1        | Sparse Matrix manipulation           |
| [Rcpp](https://cran.r-project.org/web/packages/Rcpp/index.html)           | 1.0.14       | High-performance C++ interface for R |

## Authors
David Blickenstorfer, *D-MATH*, davidbl@student.ethz.ch

## Licence
This project was developed at the **Molecular Evolution and Development Division**, Department of Neuroscience and developmental biology, *University of Vienna* under the supervision of Dr. Alison Cole & Dr. Juan Daniel Montenegro Cabrera as part of a Bachelor's thesis submitted to *ETH Zurich*, under the supervision of Prof. Dr. Valentina Boeva, **Boeva Lab**.

## 📚 References

[1]: [Altenhoff et al. OMA standalone: orthology inference among public and custom genomes and transcriptomes. *Genome Research*, 2019, 29:1152-1163](https://genome.cshlp.org/content/29/7/1152)

[2]: [Emms, D.M. and Kelly, S. (2015) OrthoFinder: solving fundamental biases in whole genome comparisons dramatically improves orthogroup inference accuracy. *Genome Biology* 16:157](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0721-2)

[3]: [Emms, D.M. and Kelly, S. (2019) OrthoFinder: phylogenetic orthology inference for comparative genomics. *Genome Biology* 20:238](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1832-y)

[4]: [Steinegger M and Soeding J. MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. *Nature Biotechnology*, doi: 10.1038/nbt.3988 (2017)](https://www.nature.com/articles/nbt.3988)

[5]: [Benjamin Buchfink, Chao Xie, and Daniel H. Huson. “Fast and sensitive protein alignment using DIAMOND”. In: *Nature Methods* 12.1 (2015), doi: 10.1038/nmeth.3176.](https://www.nature.com/articles/s41467-018-04964-5)

[6]: [R Core Team (2024). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria](https://www.R-project.org/)

[7] [Hao, Y., Hao, S., Andersen-Nissen, E., Mauck III, W.M., Zheng, S., Butler, A., Lee, M.J., Wilk, A.J., Darby, C., Zager, M., et al. (2021). Integrated analysis of multimodal single-cell data. *Cell*, 184(13), 3573–3587.e29. DOI: 10.1016/j.cell.2021.04.048](https://doi.org/10.1016/j.cell.2021.04.048)