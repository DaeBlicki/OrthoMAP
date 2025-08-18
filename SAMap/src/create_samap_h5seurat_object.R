#' ----------------------------------------------------------------------------
#' @title:      Create SAMap Object as h5Seurat Object
#'
#' @description
#' Create R h5seurat object using path of Seurat Object. Intended running as
#' Rscript: `Rscript create_samap_h5seurat_object.R`
#'
#' @author  David Blickenstorfer, Technau Group (2025)
#' @created 07/08/2025
#'
#' @copyright Technau Group, University of Vienna, 2025
#' @copyright Boeva Lab, ETH Zurich, 2025
# -----------------------------------------------------------------------------
# convertion
library(SeuratDisk)
library(Seurat)

SeuratDisk::Convert("results/samap.h5ad", dest = "h5seurat", overwrite = TRUE)
samap_obj <- SeuratDisk::LoadH5Seurat("results/samap.h5seurat",
                                      images = FALSE,
                                      meta.data = FALSE,
                                      misc = FALSE)
samap_md <- read.csv("results/raw_md.csv")
samap_md$species <- factor(as.character(samap_md$species),
                           levels = unique(as.character(samap_md$species)))
samap_md$IDs <- factor(as.character(samap_md$IDs),
                       levels = unique(as.character(samap_md$IDs)))
samap_md$ID.separate <- factor(as.character(samap_md$ID.separate),
                               levels = unique(samap_md$ID.separate))
samap_counts <- samap_obj[["RNA"]]@counts

# Create new samap object as seurat object
new_samap_obj <- Seurat::CreateSeuratObject(counts = samap_counts,
                                            meta.data = samap_md)
umap <- Seurat::Embeddings(samap_obj, reduction = "umap")
umap_reduction <- Seurat::CreateDimReducObject(embeddings = umap, key = "umap")
new_samap_obj[["umap"]] <- umap_reduction
Loadings(object = new_samap_obj[["umap"]]) <- umap

save(new_samap_obj, file = "results/samap.Robj")