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
  message(
    "❌ Invalid usage of create_oma_table.R [directory]\n",
    "[directory] : input path to `results/OMA/`\n"
  )
  stop("Usage: create_oma_table.R [directory]")
}
input_directory <- args[1]
output_directory <- args[2]

# convertion
library(Seurat)

# Create Clean Seurat Object
cat("[Loading] started at", date(), "\n")
seurat_obj <- get(load(input_directory))
cat("[Counts] started at", date(), "\n")
counts <- Seurat::GetAssayData(seurat_obj, layer = "counts")
cat("[Convert] started at", date(), "\n")
counts_df <- as.data.frame(as.matrix(counts))
cat("[Save] started at", date(), "\n")
write.csv(counts_df, file = output_directory)