#'-----------------------------------------------------------------------------
#' @title   Repeat OrthoMAP Clustering
#'
#' @concept
#' It use the current configuration and resulting OrthoMAP object and replicate
#' the results for 20 different random states in Leiden algorithm. The Output
#' stores the orthomap clusters for each cell as .csv file.
#'
#' Usage: `Rscript measure_orthomap_cluster [-i] [-o]`
#'    [-i]: Input path  as .yaml, example: "oma.yaml"
#'    [-o]: Filename as .csv, example: "oma_clusters.csv"
#'
#' @author  David Blickenstorfer, Technau Group (2025)
#' @created 11/08/2025
#'
#' Version: v1.0.0
#' Changelog:
#'  - v1.0.0: Initial release
#'
#' @copyright Technau Group, University of Vienna, 2025
#' @copyright Boeva Lab, ETH Zurich, 2025
#' ----------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage as Rscript: Rscript measure_orthomap_cluster.R [-i] [-o]\n",
       "   [-i]  : Input path as .yaml, example: 'oma.yaml'\n",
       "   [-o]  : Input path as .csv, example: 'oma_clusters.csv'\n\n")
}
# Load library
library(yaml)       # yaml (v.2.3.10)
library(Seurat)     # seurat5 (v.5.2.1), SeuratObject (5.1.0)
library(tidyverse)  # tidyverse (v.2.0.0), use dplyr, purrr

# read input variables
config <- yaml::read_yaml(args[1])
csv_name <- args[2]
output_directory <- file.path("results", "benchmark", "leiden_clusters")
dir.create(output_directory,
           recursive = TRUE,
           showWarnings = FALSE)
filename <- file.path(output_directory, csv_name)

# source top-down approach
source("src/OrthoMAP_print.R")
source("src/OrthoMAP_toolbox.R")
source("src/create_orthomap_object.R")
source("src/visualize_orthomap_object.R")
source("src/visualize_orthomap_result.R")
source("src/analyze_orthomap_result.R")
source("src/run_standard_routine.R")
source("src/run_topdown_routine.R")

# print session info
sessionInfo()
cat(yaml::as.yaml(config))

# Get orthomap result
tmp_name <- config$output$orthomap_tmp
orthomap_obj_path <- config$output$orthomap_seurat_object_path
file_path <- file.path(orthomap_obj_path, tmp_name)
orthomap_obj <- get(load(file_path))
num_iterations <- 20

cat("------------------------------------------------------------------\n")
cat("Measure OrthoMAP Clusters started ", date(), "\n")
cat("------------------------------------------------------------------\n\n")
# Step 2: Loop over random seed
cluster_list <- vector("list", num_iterations)
for (x in 1:num_iterations) {
  cat("Current Random State : ", x, "started at", date(), "\n")
  config$orthomap_top_down$fine$cluster.seed <- x
  result_obj <- run_topdown_routine(orthomap_obj = orthomap_obj,
                                    config = config,
                                    verbose = FALSE)
  cluster_list[[x]] <- result_obj$orthomap_clusters
  rm(result_obj)
}
colnames <- c("Cell_ID", paste0("Iterations_", 1:20))
df <- as.data.frame(do.call(cbind, cluster_list))
df <- cbind(rownames(df), df)
colnames(df) <- colnames

# Save csv
write.csv(df, filename, row.names = FALSE)

cat("------------------------------------------------------------------\n")
cat("Measurements succesfully finished ", date(), "\n")
cat("------------------------------------------------------------------\n")