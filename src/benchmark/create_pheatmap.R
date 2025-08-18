#'-----------------------------------------------------------------------------
#' @title   Run Pheatmap for OrthoMAP and SAMap
#'
#' @concept
#' This script visualize the Pheatmap (Cluster vs. IDs) from OrthoMAP and
#' SAMap.
#'
#' @author  David Blickenstorfer, Technau Group (2025)
#' @created 16/08/2025
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
plotnames <- c("OMA", "MMseqs2 [HOG]", "Diamond [HOG]",
               "MMseqs2 [OOG]", "Diamond [HOG]", "SAMap")
tissue_config_path <- file.path("data", "tissue_palette.csv")
palette <- read.csv(tissue_config_path, stringsAsFactors = TRUE)
samap_clusters_path <- file.path("results", "benchmark",
                                 "leiden_clusters", "samap_2.0.csv")
samap_clusters <- read.csv(samap_clusters_path)
tissue_colors <- setNames(as.character(palette$Color), palette$Tissue)
coverage <- 0.8
set.seed(42)
# -----------------------------------------------------------------------------
#' @title create_shannon_plot
#' @brief Create shannon entropy plots for given file
#' @param filename Path to the Seurat Object
#' @param plotname Title of the subplot
#' @return ggplot2 UMAP plot
# -------------------------------------------------------------------------
create_shannon_plot <- function(file, plotname) {
  # Step 1: Loading the file
  cat("Current Process: ", plotname, "\n")
  cat("Started at ", date(), "\n\n")
  obj <- get(load(file))
  if (plotname == "SAMap") {
    cluster_iteration <- sample(2:21, 1)
    representative_cluster <- samap_clusters[, cluster_iteration]
    obj$orthomap_clusters <- factor(representative_cluster,
                                    levels = unique(representative_cluster))
  }
  ct_raw_table <- as.data.frame.matrix(table(obj$orthomap_clusters,
                                             obj$IDs))
  ct_raw_table <- ct_raw_table[order(as.integer(rownames(ct_raw_table))), ]

    # Tissue vs Cluster: scaled-tissue
  plotname <- paste0("Pheatmap_Tissue_Scaling", result_format)
  filename <- file.path(output_directory, plotname)
  gt <- pheatmap::pheatmap(ct_raw_table,
                           cluster_rows = TRUE, cluster_cols = TRUE,
                           scale = "column",
                           display_numbers = TRUE,
                           fontsize_number = 4,
                           fontsize_row = 4,
                           main = "Tissue vs. Cluster: column-scaled",
                           filename = filename)$gtable
  ggplot2::ggsave(filename, plot = gt)
  # Tissue vs Cluster: scaled-cluster
  plotname <- paste0("Pheatmap_Cluster_Scaling", result_format)
  filename <- file.path(output_directory, plotname)
  gt <- pheatmap::pheatmap(ct_raw_table,
                           cluster_rows = TRUE, cluster_cols = TRUE,
                           display_numbers = TRUE,
                           scale = "row",
                           fontsize_number = 4,
                           fontsize_row = 4,
                           main = "Tissue vs. Cluster: row-scaled",
                           filename = filename)$gtable
  ggplot2::ggsave(filename, plot = gt)
  # Tissue vs Cluster: percentage tissue in cluster
  plotname <- paste0("Pheatmap_pct_clus_in_tissue", result_format)
  filename <- file.path(output_directory, plotname)
  gt <- pheatmap::pheatmap(ct_pct_table,
                           cluster_rows = TRUE, cluster_cols = TRUE,
                           display_numbers = TRUE,
                           fontsize_number = 4,
                           main = "Cluster [%] in Tissue",
                           filename = filename)$gtable
  ggplot2::ggsave(filename, plot = gt)
  # Tissue vs Cluster: percentage tissue in cluster
  plotname <- paste0("Pheatmap_pct_tissue_in_clust", result_format)
  filename <- file.path(output_directory, plotname)
  gt <- pheatmap::pheatmap(clust_ct_pct_table,
                           cluster_rows = TRUE, cluster_cols = TRUE,
                           display_numbers = TRUE,
                           fontsize_number = 4,
                           main = "Tissue [%] in Cluster",
                           filename = filename)$gtable
  ggplot2::ggsave(filename, plot = gt)
  # plot species per cluster pheatmap
  plotname <- paste0("Pheatmap_species_in_clust", result_format)
  filename <- file.path(output_directory, plotname)
  gt <- pheatmap::pheatmap(clust_sp_pct_table,
                           cluster_rows = TRUE, cluster_cols = TRUE,
                           display_numbers = TRUE,
                           main = "Species Distribution in Cluster",
                           filename = filename)$gtable
  ggplot2::ggsave(filename, plot = gt)
  
  pie_plots <- apply(ct_raw_table, 1, estimate_cluster_plots)
  # return the shannon_plots
  patchwork::wrap_plots(pie_plots, ncol = 3)
}
result_plots <- mapply(create_shannon_plot, filenames, plotnames)
# saving the plot
output_directory <- file.path("results", "benchmark", "donut_plots")
dir.create(output_directory,
           recursive = TRUE,
           showWarnings = FALSE)
ggplot2::ggsave(file.path(output_directory, "oma.pdf"), result_plots[[1]],
                width = 12, height = 14)
ggplot2::ggsave(file.path(output_directory, "mmseqs2_hog.pdf"),
                result_plots[[2]], width = 12, height = 14)
ggplot2::ggsave(file.path(output_directory, "diamond_hog.pdf"),
                result_plots[[3]], width = 12, height = 14)
ggplot2::ggsave(file.path(output_directory, "mmseqs2_oog.pdf"),
                result_plots[[4]], width = 12, height = 14)
ggplot2::ggsave(file.path(output_directory, "diamond_oog.pdf"),
                result_plots[[5]], width = 12, height = 14)
ggplot2::ggsave(file.path(output_directory, "samap.pdf"),
                result_plots[[6]], width = 12, height = 14)