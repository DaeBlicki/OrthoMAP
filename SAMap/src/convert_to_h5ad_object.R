#' ----------------------------------------------------------------------------
#' @title:      Convert to python h5ad object
#'
#' @description
#' Create python h5ad object using path of Seurat Object. Intended running as
#' Rscript: `Rscript convert_to_h5ad_object.R [-input] [-output]`
#'    [-input]  : Path to the Seurat Object as `.Robj`
#'    [-output] : Path to the resulting `.h5ad`
#'
#' Rscript convert_to_h5ad_object.R Hv.Robj Hv.h5Seurat
#'
#' @author  David Blickenstorfer, Technau Group (2025)
#' @created 04/08/2025
#'
#' @copyright Technau Group, University of Vienna, 2025
#' @copyright Boeva Lab, ETH Zurich, 2025
# -----------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript convert_to_h5ad_object.R [-input] [-output]")
}
input_path <- args[1]
output_path <- args[2]

# convertion
library(SeuratDisk)
library(Seurat)

# Load and change gene names
seurat_obj <- get(load(input_path))

# CleanUp
DefaultAssay(seurat_obj) <- "RNA"
seurat_obj@assays <- list(RNA = seurat_obj@assays$RNA)
seurat_obj <- Seurat::NormalizeData(seurat_obj,
                                    normalization.method = "LogNormalize",
                                    verbose = TRUE)
seurat_obj@reductions <- list()
seurat_obj@commands <- list()
seurat_obj@tools <- list()
seurat_obj@graphs <- list()
seurat_obj@assays$RNA@scale.data <- matrix()
Seurat::VariableFeatures(seurat_obj) <- rownames(seurat_obj)

# save as .h5ad
SeuratDisk::SaveH5Seurat(seurat_obj, filename = output_path)
SeuratDisk::Convert(output_path, dest = "h5ad", overwrite = TRUE,
                    X = "data", counts = "counts")