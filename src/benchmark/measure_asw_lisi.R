#'-----------------------------------------------------------------------------
#' @title   Measure ASW and LISI
#'
#' @concept
#' This script measures Adjusted Silhouette Score and Local Inverse Simpson
#' Index for given OrthoMAP object.
#'
#' Usage: `Rscript measure_asw_lisi [-o]`
#'    [-i]: Input path as .Robj, example: "oma.Robj"
#'    [-o]: Output path as .csv, example: "oma_asw_lisi.csv"
#'
#' @details
#' The computational intensity and cost forces to reduce the number of cells
#' in the Embeddings and Distance matrix. This scripts samples randomely 80%
#' of the total data set for ASW and LISI. This is repeated 20 times and the
#' average for the silhouette width is measured and the median for the inverse
#' Simpson index is evaluated and stored in .csv file.
#'
#' see more https://doi.org/10.1093/nar/gkae1316
#'
#' @author  David Blickenstorfer, Technau Group (2025)
#' @created 15/08/2025
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
  stop("Usage as Rscript: Rscript benchmark_orthomap.R [-o]\n",
       "   [-i]  : Input path as .yaml, example: 'oma.Robj'\n",
       "   [-o]  : Filename as .csv, example: 'oma_asw_lisi.csv'\n\n")
}
library(tidyverse)
library(Seurat)
library(cluster)
library(proxy)
if (!requireNamespace("lisi")) {
  library(devtools, lib.loc = "~/Rlibs")
  library(lisi, lib.loc = "~/Rlibs")
  requireNamespace("lisi")
}
# -----------------------------------------------------------------------------
# Step 1: Load Data and Preparation
# -----------------------------------------------------------------------------
cat("[Step 1] Loading, started at ", date(), "\n")
orthomap_obj_path <- args[1]
csv_name <- args[2]
output_directory <- file.path("results", "benchmark", "metrics")
dir.create(output_directory,
           recursive = TRUE,
           showWarnings = FALSE)
filename <- file.path(output_directory, csv_name)
orthomap_obj <- get(load(orthomap_obj_path))

# -----------------------------------------------------------------------------
# Step 2: Estimates ASW and LISI
# -----------------------------------------------------------------------------
cat("[Step 2] Estimate ASW and LISI, staretd at ", date(), "\n")
# -----------------------------------------------------------------------------
#' @title estimate_asw_and_lisi
#' @brief estimate the ASW and LISI using percent_extract of total cells
#' @param orthomap_obj OrthoMAP Object as Seurat Object
#' @param percent_extract Numeric, extract percentage of cells
#' @param num_iterations Integer, number of measurements
#' @param random_seed Integer, reproducibility of experiment/ measurements
#' @return Matrix (metrics x iterations)
# -------------------------------------------------------------------------
estimate_asw_and_lisi <- function(orthomap_obj,
                                  percent_extract = 0.25,
                                  num_iterations = 20,
                                  random_seed = 42) {
  # get meta data and remove memory of orthomap object
  if (percent_extract > 1.0 || percent_extract <= 0.0) {
    stop("Invalid parameter in 'estimate_asw_and_clisi'\n")
  }
  set.seed(random_seed)
  num_cells <- ncol(orthomap_obj)
  num_pseudo_cells <- ceiling(percent_extract * num_cells)
  cat("Number of used cells: ", num_pseudo_cells, "/", num_cells, "\n")
  metrics <- lapply(1:num_iterations, function(iteration) {
    cat("Current Iteration: ", iteration, "\n")
    # sample cells randomely
    pseudo_cells <- sample(1:num_cells, num_pseudo_cells, replace = FALSE)
    embeddings <- Seurat::Embeddings(orthomap_obj,
                                     reduction = "harmony")[pseudo_cells, ]
    # Process metadata on the pseudo cells
    metadata <- orthomap_obj@meta.data[pseudo_cells, ]
    metadata$IDs <- factor(metadata$IDs,
                           levels = unique(metadata$IDs))
    metadata$IDs <- as.numeric(metadata$IDs)
    metadata$orig.ident <- factor(metadata$orig.ident,
                                  levels = unique(metadata$orig.ident))
    metadata$orig.ident <- as.numeric(metadata$orig.ident)
    # Estimate ASW scores
    dist <- proxy::dist(embeddings, method = "euclidean")
    sil_bio <- cluster::silhouette(metadata$IDs, dist)
    sil_batch <- cluster::silhouette(metadata$orig.ident, dist)
    asw_bio <- (mean(sil_bio[, 3]) + 1) / 2
    asw_batch <- 1 - ((mean(sil_batch[, 3]) + 1) / 2)
    rm(dist)
    # Estimate LISI scores
    labels <- c("IDs", "orig.ident")
    lisi <- lisi::compute_lisi(embeddings,
                               metadata,
                               label_colnames = labels,
                               perplexity = 40,
                               nn_eps = 0.0)
    # normalize the LISI score for each cell
    n_celltype <- length(unique(metadata$IDs))
    n_batch <- length(unique(metadata$orig.ident))
    norm_lisi_bio <- (n_celltype - lisi$IDs) / (n_celltype - 1)
    norm_lisi_batch <- (lisi$orig.ident - 1) / (n_batch - 1)
    lisi_bio <- median(norm_lisi_bio)
    lisi_batch <- median(norm_lisi_batch)
    # return the results as follows: cASW, iASW, cLISI, iLISI
    rm(embeddings, metadata, norm_lisi_bio, norm_lisi_batch)
    gc()
    names <- c("cASW", "iASW", "cLISI", "iLISI")
    metrics <- c(asw_bio, asw_batch, lisi_bio, lisi_batch)
    names(metrics) <- names
    metrics
  })
  as.data.frame(do.call(rbind, metrics))
}
metrics <- estimate_asw_and_lisi(orthomap_obj)
# -----------------------------------------------------------------------------
# Step 3: Save as .csv file
# -----------------------------------------------------------------------------
write.csv(metrics, filename)