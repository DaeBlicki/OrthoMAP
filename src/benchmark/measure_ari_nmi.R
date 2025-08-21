#'-----------------------------------------------------------------------------
#' @title   Measure ARI and NMI
#'
#' @concept
#' This script measures Adjusted Rand Index and Normalized Mutual Information
#' for given OrthoMAP object.
#'
#' Usage: `Rscript measure_ari_nmi [-i] [-o]`
#'    [-i]: Filename as .csv, example: "oma_clusters.csv"
#'    [-o]: Output path as .csv, example: "oma_ari_nmi.csv"
#'
#' @details
#' The ARI and NMI is measures the overlapping between the estimated orthomap
#' clusters and the celltype ("IDs") and the batch ("orig.ident").
#'
#' see more https://doi.org/10.1093/nar/gkae1316
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
       "   [-o]  : Filename as .csv, example: 'oma_ari_nmi.csv'\n\n")
}
library(tidyverse)
library(aricode)
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
# Step 2: Estimation of Metrics using 20 measurements
# - ARI (Adjusted Random Index)
# - NMI (Normalized Mutual Index)
# -----------------------------------------------------------------------------
cat("[Step 2] Estimate ARI and NMI, started at ", date(), "\n")
md_colnames <- c("Cell_ID", "IDs", "ID.separate", "species", "orig.ident")
# calculate the table
metrics <- apply(data[!names(data) %in% md_colnames], 2, function(x) {
  ari_bio <- aricode::ARI(x, data[, "IDs"])
  ari_batch <- 1 - aricode::ARI(x, data[, "orig.ident"])
  nmi_bio <- aricode::NMI(x, data[, "IDs"])
  nmi_batch <- 1 - aricode::NMI(x, data[, "orig.ident"])
  names <- c("cARI", "iARI", "cNMI", "iNMI")
  metrics <- c(ari_bio, ari_batch, nmi_bio, nmi_batch)
  names(metrics) <- names
  metrics
})
# -----------------------------------------------------------------------------
# Step 3: Save as .csv file
# -----------------------------------------------------------------------------
write.csv(metrics, filename)