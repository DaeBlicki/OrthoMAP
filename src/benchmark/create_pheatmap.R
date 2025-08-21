#'-----------------------------------------------------------------------------
#' @title   Run Pheatmap for OrthoMAP and SAMap
#'
#' @concept
#' This script visualize the Pheatmap (Cluster vs. IDs) from OrthoMAP and
#' SAMap.
#'
#' @author  David Blickenstorfer, Technau Group (2025)
#' @created 18/08/2025
#'
#' Version: v1.0.0
#' Changelog:
#'  - v1.0.0: Initial release
#'
#' @copyright Technau Group, University of Vienna, 2025
#' @copyright Boeva Lab, ETH Zurich, 2025
#' ----------------------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(patchwork)
library(pheatmap)
# -----------------------------------------------------------------------------
# Step 1: Load Data and Preparation
# -----------------------------------------------------------------------------
seurat_obj_path <- file.path("results", "OrthoMAP_Seurat_Objects")
filenames <- c(
  file.path(seurat_obj_path, "oma", "result.Robj"),
  file.path(seurat_obj_path, "mmseqs2", "HOG", "result.Robj"),
  file.path(seurat_obj_path, "diamond", "HOG", "result.Robj"),
  file.path(seurat_obj_path, "mmseqs2", "OOG", "result.Robj"),
  file.path(seurat_obj_path, "diamond", "OOG", "result.Robj"),
  file.path("results", "SAMap", "samap.Robj")
)
output_directory <- file.path("results", "benchmark", "pheatmaps")
dir.create(output_directory,
           recursive = TRUE,
           showWarnings = FALSE)
plotnames <- c("oma.pdf", "mmseqs2_hog.pdf", "diamond_hog.pdf",
               "mmseqs2_oog.pdf", "diamond_oog.pdf", "samap.pdf")
plotnames <- file.path(output_directory, plotnames)
coarsenames <- c("oma_coarse.pdf", "mmseqs2_hog_coarse.pdf",
                 "diamond_hog_coarse.pdf", "mmseqs2_oog_coarse.pdf",
                 "diamond_oog_coarse.pdf", "")
coarsenames <- file.path(output_directory, coarsenames)
samap_clusters_path <- file.path("results", "benchmark",
                                 "leiden_clusters", "samap_2.0.csv")
samap_clusters <- read.csv(samap_clusters_path)
set.seed(42)
# -----------------------------------------------------------------------------
#' @title create_shannon_plot
#' @brief Create shannon entropy plots for given file
#' @param filename Path to the Seurat Object
#' @param plotname Title of the subplot - orthomap clusters
#' @param coarsename Title of the subplot - coarse clusters
#' @return ggplot2 UMAP plot
# -------------------------------------------------------------------------
create_pheatmap <- function(file, plotname, coarsename) {
  # Step 1: Loading the file
  cat("Current Process: ", plotname, "\n")
  cat("Started at ", date(), "\n\n")
  obj <- get(load(file))
  if (plotname == plotnames[[6]]) {
    cluster_iteration <- sample(2:21, 1)
    representative_cluster <- samap_clusters[, cluster_iteration]
    obj$orthomap_clusters <- factor(representative_cluster,
                                    levels = unique(representative_cluster))
  }
  ct_raw_table <- as.data.frame.matrix(table(obj$orthomap_clusters,
                                             obj$IDs))

  # Tissue vs Cluster: scaled-tissue
  gt_col_scale <- pheatmap::pheatmap(ct_raw_table,
                                     cluster_rows = TRUE, cluster_cols = TRUE,
                                     treeheight_row = 0,
                                     treeheight_col = 0,
                                     scale = "column",
                                     display_numbers = TRUE,
                                     fontsize_number = 4,
                                     fontsize_row = 4,
                                     main = "column-scaled")$gtable

  gt_row_scale <- pheatmap::pheatmap(ct_raw_table,
                                     cluster_rows = TRUE, cluster_cols = TRUE,
                                     treeheight_row = 0,
                                     treeheight_col = 0,
                                     display_numbers = TRUE,
                                     scale = "row",
                                     fontsize_number = 4,
                                     fontsize_row = 4,
                                     main = "row-scaled")$gtable
  gt <- gridExtra::grid.arrange(gt_row_scale, gt_col_scale, nrow = 2)
  ggplot2::ggsave(plotname, plot = gt, width = 12, height = 14)

  if (plotname != plotnames[[6]]) {
    coarse_raw_table <- as.data.frame.matrix(table(obj$coarse_clusters,
                                                   obj$IDs))
    gt_col_scale <- pheatmap::pheatmap(coarse_raw_table,
                                       cluster_rows = TRUE,
                                       cluster_cols = TRUE,
                                       treeheight_row = 0,
                                       treeheight_col = 0,
                                       scale = "column",
                                       display_numbers = TRUE,
                                       fontsize_number = 4,
                                       fontsize_row = 4,
                                       main = "column-scaled")$gtable

    gt_row_scale <- pheatmap::pheatmap(coarse_raw_table,
                                       cluster_rows = TRUE, cluster_cols = TRUE,
                                       treeheight_row = 0,
                                       treeheight_col = 0,
                                       display_numbers = TRUE,
                                       scale = "row",
                                       fontsize_number = 4,
                                       fontsize_row = 4,
                                       main = "row-scaled")$gtable
    gt <- gridExtra::grid.arrange(gt_row_scale, gt_col_scale, nrow = 2)
    ggplot2::ggsave(coarsename, plot = gt, width = 12, height = 14)
  }
  rm(obj)
}
invisible(mapply(create_pheatmap, filenames, plotnames, coarsenames))