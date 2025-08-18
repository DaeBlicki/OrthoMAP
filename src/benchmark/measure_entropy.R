#'-----------------------------------------------------------------------------
#' @title   Measure normalized Shannon Entropy
#'
#' @concept
#' This script measures the normalized Shannon Entropy from the given OrthoMAP
#' Object and the estimated orthomap clusters.
#'
#' Usage: `Rscript measure_ari_nmi [-i] [-o]`
#'    [-i]: Filename as .csv, example: "oma_clusters.csv"
#'    [-o]: Output path as .csv, example: "oma_entropy.csv"
#'
#' @details
#' The normalized Shannon entropy is measured for each orthomap cluster
#' estimated in each iteration. The normalized Shannon entropy is returned
#' as one vector that can be presented as boxplot.
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
       "   [-i]  : Input path as .csv, example: 'oma_clusters.csv'\n",
       "   [-o]  : Filename as .csv, example: 'oma_entropy.csv'\n\n")
}
library(tidyverse)
library(entropy)
# -----------------------------------------------------------------------------
# Step 1: Load Data and Preparation
# -----------------------------------------------------------------------------
file <- args[1]
csv_name <- args[2]
cat("[Step 1] Load and merge data, started at ", date(), "\n")
output_directory <- file.path("results", "benchmark", "metrics")
dir.create(output_directory,
           recursive = TRUE,
           showWarnings = FALSE)
filename <- file.path(output_directory, csv_name)
md <- read.csv(file.path("results", "benchmark", "raw_md.csv"), row.names = 1)
md$Cell_ID <- rownames(md)
data <- read.csv(file)
data <- merge(data, md, by = "Cell_ID", all.x = TRUE)

# -----------------------------------------------------------------------------
# Step 1: Load Data and Preparation
# -----------------------------------------------------------------------------
cat("[Step 2] Estimate normalized shannon entropy")
md_colnames <- c("Cell_ID", "IDs", "ID.separate", "species", "orig.ident")
# -----------------------------------------------------------------------------
#' @title get_norm_entropy
#' @brief generate harmonized/normalized Shannon Entropy
#' @param cluster_vec Vector, contains cluster assignement
#' @param tissue_vec Vector, contains tissue assignement
#' @return Harmonized Shannon Entropy [0, 1]
# -------------------------------------------------------------------------
get_norm_entropy <- function(cluster_vec, tissue_vec) {
  # calculate percentage of tissue for each cluster
  tissue_raw_table <- as.data.frame.matrix(table(cluster_vec, tissue_vec))
  size_vector <- apply(tissue_raw_table, 1, sum)
  tissue_pct_table <- tissue_raw_table / size_vector
  # evaluate harmonized shannon entropy for each cluster
  entropy_vector <- apply(tissue_pct_table, 1, function(pct_in_cluster) {
    entropy::entropy.empirical(pct_in_cluster, unit = "log")
  })
  harmonized_entropy_vector <- entropy_vector / log(size_vector)
  harmonized_entropy_vector[size_vector == 1] <- 0.0
  harmonized_entropy_vector
}
entropy <- apply(data[, !(names(data) %in% md_colnames)], 2, function(clust) {
  get_norm_entropy(clust, data$IDs)
})
# -----------------------------------------------------------------------------
# Step 3: Save as .csv file
# -----------------------------------------------------------------------------
entropy_df <- data.frame(
  entropy = as.numeric(unlist(entropy))
)
write.csv(entropy_df, filename)