# Project Journal

This **Markdown** have the purpose to replicate the results of the Thesis and provide a step-by-step guide.

## 1. Download and pre-process single-cell RNA sequencing atlas

### Acces single-cell atlas and FASTA files

Please contact Dr. Cole A.G (alison.cole@univie.ac.at) or Dr. Montenegro J. (juan.montenegro@univie.ac.at) for data access for *Nematostella vectensis* and *Aurelia coerulea*. Ask for the Annotation table `.csv`, Fasta file `.fa`, and single-cell atlas `.Robj` or `.rds`. Place the annotiation table in `data/annotation_table`, the `.fa` in `data/DB/` and the `.Robj` in `data/single_cell_atlas`. After the data access was sucessful, run following in the project root:

```
source scripts/OrthoMAP/01_download.sh
source scripts/OrthoMAP/02_preprocess_data.sh
```

1) It downloads and pre-processs the scRNA-seq data from *Hydra vulgaris*. The script will download [`aepAtlasDubInclude.rds`](https://pubmed.ncbi.nlm.nih.gov/36639202/) and the list of the non-doublets. In the next step, possible doublets are getting removed (see [05_hydraAtlasReMap.md](https://github.com/cejuliano/brown_hydra_genomes/blob/main/05_hydraAtlasReMap.md)) as described in their paper. Additionally, the annotation table is created to match `.fa` and `.Robj`. Last but not least, it creates a new meta.data `IDs` by removing the prefix: `I_`, `Ec_`, and `En_` in `CuratedIdent` and update the Seurat Object to the newest version. 

2) It removes the isoforms in Aurelias protein model by choosing the longest sequence per geneID. 

## 2. Run orthology interference analysis
The orthology interference was analyzed with OMA standalone and OrthoFinder (Diamond and MMseqs2 algorithm). Run ```source scripts/OrthoMAP/03_run_orthology.sh``` in the project repository. This can take a while (2-10h, for each analyis using 64 AMD cores). **After** the results are produced, run ```source scripts/OrthoMAP/04_preprocess_orthology.sh``` to evaluate the orthologs table used in the OrthoMAP workflow.

## 3. Run OrthoMAP
In `data/configuration/` are `config.yaml` used to generate the best results in the project.

- `05a` : OMA analysis on OMA groups
- `05b` : OrthoFinder MMseqs2 on HOG (Hierarchical Orthologous Groups)
- `05c` : OrthoFinder MMseqs2 on OOG (1:1 Orthologs)
- `05d` : OrthoFinder Diamond on HOG
- `05e` : OrthoFinder Diamond on OOG

For example, if you want to load OMA configuration, run following command in project root ```source scripts/OrthoMAP/05a_load_oma_config.sh```. After the configuration is loaded, run ```sbatch < scripts/OrthoMAP/05_run_orthomap.sh``` in SLURM or ````source scripts/OrthoMAP/05_run_orthomap.sh``` on the local device.

The project used 9 runs to generate the results and plots in the thesis. Please, use the same configuration as shown in the `.log` files. The code for the plots and results is hardcoded, make sure to respect the directory order and numbering.

## 4. Run SAMap
The comparison with SAMap requires the result from SAMap. For that, go to the `SAMap/` subdirectory and copy (or use symbolic links) the `.fa` files in `data/DB/`, the `.Robj` files in `data/single_cell_atlas/`, and the `.csv` files in `data/annotation_table/`. Next, use the first bash script to install SAMap and prepare it dependencies. Below is a short code snipped.

```
cd SAMap

# Make symbolic links with Seurat files
ln -s ../data/single_cell_atlas/Ac.Robj data/single_cell_atlas/Ac.Robj
ln -s ../data/single_cell_atlas/Nv.Robj data/single_cell_atlas/Nv.Robj
ln -s ../data/single_cell_atlas/Hv.Robj data/single_cell_atlas/Hv.Robj

# Make symbolic links with FASTA files
ln -s ../data/DB/Ac.fa data/DB/Ac.fa
ln -s ../data/DB/Nv.fa data/DB/Nv.fa
ln -s ../data/DB/Hv.fa data/DB/Hv.fa

# Make symbolic links with annotation tables files
ln -s ../data/annotation_table/Ac.csv data/annotation_table/Ac.csv
ln -s ../data/annotation_table/Nv.csv data/annotation_table/Nv.csv
ln -s ../data/annotation_table/Hv.csv data/annotation_table/Hv.csv

# install SAMap from GitHub
source scripts/01_install.sh
```

After SAMap is installed, run the command to convert .Robj to .h5ad files. Python uses scanpy and thus, needs AnnData format.
```
sbatch < 02_convert_to_h5ad.sh
```

SAMap depends on BLAST results. Run following command.
```
sbatch < 03_run_blast.sh
```

After the BLAST run is finished and the files are converted, run following command for the SAMap analysis. This will generate a SAMap Object.
```
sbatch < 04_run_samap.sh
```

The benchmark uses the best SAMap result in dependency of Leiden clustering. However, Leiden uses a random seed and thus, it measures the Leiden clusters 20 times and stores the results in `SAMap/results/leiden_clusters/`
```
sbatch < 05_run_measurements.sh
```

To estimate the best SAMap results, it uses the sum of cARI, iARI, cNMI, and iNMI of the 20 measurements. In the .out file, the results is shown. In my case, it was the Leiden Clustering using the resolution: 2.0 
```
sbatch < 06_estimate_best.sh
```

Run the next line to create a Seurat Object of the SAMap result. It contains the meta data and UMAP embeddings.
```
sbatch < 07_convert_to_h5seurat.sh
```

## 4. Comparison with SAMap
Copy paste the measurements of the best SAMap result in `results/benchmark/leiden_clusters` and the `samap_md.csv` and `samap.Robj` in `results/SAMap/`. Next, go to the root of the project repository and run the scripts in `scripts/benchmark/` in the described order. Make sure, that each job is finished before going to the next number. The measurements can be done in parallel without any issues.