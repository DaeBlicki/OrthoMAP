#' ----------------------------------------------------------------------------
#' @title visualize_orthomap_object
#'
#' @description
#' Generate visualization results for the final OrthoMAP Seurat Object using
#' Seurat5 visualization tool.
#'
#' @details
#'  - ElbowPlot: find out which dimensionality PCA needs
#'  - DimPlot PCA and Harmony
#'  - DimPlot UMAP.pca and UMAP.harmony
#'
#' @param orthomap_obj The Seurat Object generated from step 1 of the pipeline.
#' @param config RObject of .yaml, parameter for coloring
#' @param output_directory path to the resulting picture output
#' @param verbose Logical, if TRUE then print information to the user
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_histogram
#' @importFrom ggplot2 geom_vline
#' @importFrom scales
#'
#' @examples
#' \dontrun{
#'   visualize_orthomap_object(
#'     orthomap_obj = orthomap_obj,
#'     config = config,
#'     output_directory = method_orthomap_qc_path,
#'     verbose = verbose,
#'   )
#' }
#'
#' @author  David Blickenstorfer, Technau Group (2025)
#' @created 15/07/2025
#'
#' @copyright Technau Group, University of Vienna, 2025
#' @copyright Boeva Lab, ETH Zurich, 2025
#'
#' @export
#' ----------------------------------------------------------------------------
visualize_orthomap_result <- function(orthomap_obj,
                                      config,
                                      output_directory,
                                      verbose = TRUE) {
  # Read configuration
  tissue_colors <- Seurat::DiscretePalette(length(levels(orthomap_obj$IDs)))
  tissue_config_path <- file.path("data", "tissue_palette.csv")
  if (file.exists(tissue_config_path)) {
    palette <- read.csv(tissue_config_path, stringsAsFactors = TRUE)
    tissue_colors <- setNames(as.character(palette$Color), palette$Tissue)
    orthomap_obj$IDs <- factor(orthomap_obj$IDs, levels = palette$Tissue)
  }
  result_format <- config$visualization$result_format

  # Create DimPlot for Harmony
  if (verbose) cat("Create Dimplot for harmony and pca \n")
  dim_harmony_plot <- Seurat::DimPlot(
    orthomap_obj,
    group.by = "orthomap_clusters",
    reduction = "harmony"
  ) + Seurat::NoLegend()
  plotname <- paste0("DimPlot_harmony", result_format)
  plotname <- file.path(output_directory, plotname)
  ggplot2::ggsave(plotname, dim_harmony_plot, width = 6, height = 4)

  # Create DimPlot for PCA
  dim_pca_plot <- Seurat::DimPlot(
    orthomap_obj,
    group.by = "orthomap_clusters",
    reduction = "pca"
  ) + Seurat::NoLegend()
  plotname <- paste0("DimPlot_pca", result_format)
  plotname <- file.path(output_directory, plotname)
  ggplot2::ggsave(plotname, dim_pca_plot, width = 6, height = 4)

  # Create ElbowPlot for Harmony
  if (verbose) cat("Create ElbowPlot for harmony and pca \n")
  elbow_harmony_plot <- Seurat::ElbowPlot(
    orthomap_obj,
    reduction = "harmony",
    ndims = config$seurat_standard$harmony.dimensions
  )
  sd <- config$seurat_standard$reduction.stdev
  elbow_harmony_plot <- elbow_harmony_plot +
    ggplot2::geom_hline(yintercept = sd,
                        linewidth = 0.5,
                        color = "black") +
    ggplot2::theme_bw()
  plotname <- paste0("ElbowPlot_harmony", result_format)
  plotname <- file.path(output_directory, plotname)
  ggplot2::ggsave(plotname, elbow_harmony_plot, width = 6, height = 4)

  # Create ElbowPlot for PCA
  elbow_pca_plot <- Seurat::ElbowPlot(
    orthomap_obj,
    reduction = "pca",
    ndims = config$seurat_standard$pca.dimensions
  )
  sd <- config$seurat_standard$reduction.stdev
  elbow_pca_plot <- elbow_pca_plot +
    ggplot2::geom_hline(yintercept = sd,
                        linewidth = 0.5,
                        color = "black") +
    ggplot2::theme_bw()
  plotname <- paste0("ElbowPlot_pca", result_format)
  plotname <- file.path(output_directory, plotname)
  ggplot2::ggsave(plotname, elbow_pca_plot, width = 6, height = 4)

  # ---------------------------------------------------------------------------
  # RunUMAP visualizations
  #   1) Cluster from findNeighbors and findClusters
  #   2) Speceis
  #   3) Orig.ident
  #   4) IDs (CellType)
  # ---------------------------------------------------------------------------
  if (verbose) cat("Create Dimplot for Cluster\n")
  umap_harmony_plot <- Seurat::DimPlot(
    orthomap_obj,
    group.by = "orthomap_clusters",
    reduction = "umap.harmony"
  ) + Seurat::NoAxes() + ggplot2::theme_bw()

  plotname <- paste0("umap_h_cluster", result_format)
  plotname <- file.path(output_directory, plotname)
  ggplot2::ggsave(plotname, umap_harmony_plot, width = 12, height = 6)

  if (verbose) cat("Create Dimplot for Species\n")
  umap_harmony_plot <- Seurat::DimPlot(
    orthomap_obj,
    group.by = "species",
    reduction = "umap.harmony"
  ) + Seurat::NoAxes() + ggplot2::theme_bw()
  plotname <- paste0("umap_h_species", result_format)
  plotname <- file.path(output_directory, plotname)
  ggplot2::ggsave(plotname, umap_harmony_plot, width = 6, height = 6)

  if (verbose) cat("Create Dimplot for orig.ident\n")
  umap_harmony_plot <- Seurat::DimPlot(
    orthomap_obj,
    group.by = "orig.ident",
    reduction = "umap.harmony"
  ) + Seurat::NoAxes() + ggplot2::theme_bw()
  plotname <- paste0("umap_h_origIdent", result_format)
  plotname <- file.path(output_directory, plotname)
  ggplot2::ggsave(plotname, umap_harmony_plot, width = 20, height = 12)

  if (verbose) cat("Create Dimplot for IDs\n")
  umap_harmony_plot <- Seurat::DimPlot(
    orthomap_obj,
    group.by = "IDs",
    reduction = "umap.harmony",
    cols = tissue_colors,
  ) + Seurat::NoAxes() + ggplot2::theme_bw()
  plotname <- paste0("umap_h_IDs", result_format)
  plotname <- file.path(output_directory, plotname)
  ggplot2::ggsave(plotname, umap_harmony_plot, width = 20, height = 12)

  # RunUMAP for pca
  umap_pca_plot <- Seurat::DimPlot(
    orthomap_obj,
    group.by = "orthomap_clusters",
    reduction = "umap.pca"
  ) + Seurat::NoAxes() + ggplot2::theme_bw()
  plotname <- paste0("UMAP_pca", result_format)
  plotname <- file.path(output_directory, plotname)
  ggplot2::ggsave(plotname, umap_pca_plot, width = 12, height = 6)
}
