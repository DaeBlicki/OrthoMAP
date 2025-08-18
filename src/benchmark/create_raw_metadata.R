#' ----------------------------------------------------------------------------
#' @title:      Create Raw Metadata as .csv Object
#'
#' @description
#' Create csv metadata for the SAMap Object. Intended running as Rscript
#'
#' Usage: `Rscript create_samap_h5seurat_object.R`
#'
#' @author  David Blickenstorfer, Technau Group (2025)
#' @created 11/08/2025
#'
#' @copyright Technau Group, University of Vienna, 2025
#' @copyright Boeva Lab, ETH Zurich, 2025
# -----------------------------------------------------------------------------
library(Seurat)
library(tidyverse)

# Load Ac, Hv, and Nv Seurat Objects
seurat_ac <- get(load("data/single_cell_atlas/Ac.Robj"))
seurat_hv <- get(load("data/single_cell_atlas/Hv.Robj"))
seurat_nv <- get(load("data/single_cell_atlas/Nv.Robj"))

# Get needed metadata
md_ac <- data.frame(
  species = "Ac",
  orig.ident = seurat_ac$orig.ident,
  IDs = paste0("Ac.", seurat_ac$IDs),
  ID.separate = seurat_ac$ID.separate
)
md_hv <- data.frame(
  species = "Hv",
  orig.ident = seurat_hv$orig.ident,
  IDs = paste0("Hv.", seurat_hv$IDs),
  ID.separate = seurat_hv$curatedIdent
)
md_nv <- data.frame(
  species = "Nv",
  orig.ident = seurat_nv$orig.ident,
  IDs = paste0("Nv.", seurat_nv$IDs),
  ID.separate = seurat_nv$ID.separate
)

# meta data
md <- rbind(md_ac, md_hv, md_nv)
write.csv(md, "results/benchmark/raw_md.csv", row.names = TRUE)
