#'-----------------------------------------------------------------------------
#' @title   Run Donut Chart Visualization for OrthoMAP and SAMap
#'
#' @concept
#' This script visualize the results from OrthoMAP and SAMap for Species, IDs
#' and orig.ident as UMAP. It shows the 80% coverage of the cluster. For SAMap,
#' it chooses randomely from the 20 iterations (use seed for reproducibility)
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
  # -------------------------------------------------------------------------
  #' @title create_cluster_plots
  #' @brief evaluate the tissue percentage of the clusters
  #' @param ct_raw_table_line line in ct_raw_table
  #' @return Character Vector, line in the empirical results table
  # -------------------------------------------------------------------------
  estimate_cluster_plots <- function(ct_raw_table_line) {
    # -------------------------------------------------------------------------
    ## 1 Create Pie Chart and evaluate number of celltypes
    # -------------------------------------------------------------------------
    n_table <- as.integer(ct_raw_table_line)
    cluster_id <- rownames(ct_raw_table)[which(apply(ct_raw_table,
                                                     1, identical,
                                                     ct_raw_table_line))]
    total_cells <- sum(n_table)
    pct_table <- n_table / total_cells
    ct_table <- colnames(ct_raw_table)
    # get order and calculate the minimal number of cells to achieve coverage
    order <- order(as.integer(n_table), decreasing = TRUE)
    ordered_n_table <- n_table[order]
    ordered_pct_table <- pct_table[order]
    ordered_ct_table <- ct_table[order]
    # create vector of most important celltype that covers "coverage"
    cumsum_coverage <- cumsum(ordered_pct_table)
    num_celltype <- which(cumsum_coverage >= coverage)[1]
    n_vector <- ordered_n_table[1:num_celltype]
    pct_vector <- ordered_pct_table[1:num_celltype]
    tissue_vector <- ordered_ct_table[1:num_celltype]
    # add the rest as cumultative labelled as "Others"
    n_vector <- c(n_vector, total_cells - sum(n_vector))
    pct_vector <- c(pct_vector, 1.0 - sum(pct_vector))
    pct_vector <- 100 * round(pct_vector, digits = 2)
    tissue_vector <- c(tissue_vector, "Others")
    # generate doughnut chart with ggplot2
    df <- data.frame(
      Tissue = tissue_vector,
      n = n_vector,
      pct = pct_vector,
      label = paste0("\n ", n_vector,
                     "\n ", pct_vector, "%")
    )
    df$ymax <- cumsum(df$pct)
    df$ymin <- c(0, head(df$ymax, n = -1))
    df$labelPosition <- (df$ymax + df$ymin) / 2
    df$Tissue <- factor(df$Tissue, levels = df$Tissue)

    pie_plot <- ggplot2::ggplot(df, aes(ymax=ymax, ymin=ymin,
                                      xmax=4, xmin=2.5,
                                      fill=Tissue)) +
      ggplot2::theme_void() +
      ggplot2::scale_fill_manual(values = tissue_colors) +
      guides(fill = guide_legend(override.aes = list(shape = 16, size = 4))) +
      ggplot2::theme(plot.background = element_rect(fill = "white",
                                                    color = NA),
                     legend.justification = "left",
                     legend.title = ggplot2::element_blank(),
                     legend.key = ) +
      ggplot2::annotate(geom = "text", x = -1, y = 0, label = cluster_id,
                        size = 10, fontface = "bold", color = "black") +
      ggplot2::geom_rect() +
      ggplot2::geom_text(x = 3.25, aes(y=labelPosition, label=label), size=2) +
      ggplot2::coord_polar(theta="y", start = 0) +
      ggplot2::xlim(c(-1, 4))
    # return pie_plot
    ggplot2::ggsave("test/donut/single.png", pie_plot,
                  width = 12, height = 14)
    pie_plot
  }
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
ggplot2::ggsave(file.path(output_directory, "oma.png"), result_plots[[1]],
                width = 12, height = 14)
ggplot2::ggsave(file.path(output_directory, "mmseqs2_hog.png"),
                result_plots[[2]], width = 12, height = 14)
ggplot2::ggsave(file.path(output_directory, "diamond_hog.png"),
                result_plots[[3]], width = 12, height = 14)
ggplot2::ggsave(file.path(output_directory, "mmseqs2_oog.png"),
                result_plots[[4]], width = 12, height = 14)
ggplot2::ggsave(file.path(output_directory, "diamond_oog.png"),
                result_plots[[5]], width = 12, height = 14)
ggplot2::ggsave(file.path(output_directory, "samap.png"),
                result_plots[[6]], width = 12, height = 14)