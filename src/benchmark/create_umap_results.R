#'-----------------------------------------------------------------------------
#' @title   Run UMAP visualization for OrthoMAP and SAMap
#'
#' @concept
#' This script visualize the results from OrthoMAP and SAMap for Species, IDs
#' and orig.ident as UMAP.
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
               "MMseqs2 [OOG]", "Diamond [OOG]", "SAMap")
tissue_config_path <- file.path("data", "tissue_palette.csv")
palette <- read.csv(tissue_config_path, stringsAsFactors = TRUE)
tissue_colors <- setNames(as.character(palette$Color), palette$Tissue)

# -----------------------------------------------------------------------------
# Step 2: Create UMAP for all results
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
#' @title create_umap
#' @brief Create UMAP for given file for celltype, species, and batch
#' @param filename Path to the Seurat Object
#' @param plotname Title of the subplot
#' @return ggplot2 UMAP plot
# -------------------------------------------------------------------------
create_umap <- function(filename, plotname) {
  cat("Current Process: ", plotname, "\n")
  cat("Started at ", date(), "\n\n")
  # check for umap reduction (SAMap has only umap)
  reduction <- "umap.harmony"
  if (plotname == "SAMap") reduction <- "umap"
  # Load the Object
  obj <- get(load(filename))
  obj$IDs <- factor(obj$IDs, levels = palette$Tissue)
  species_plot <- Seurat::DimPlot(obj, group.by = "species",
                                  reduction = reduction) +
    Seurat::NoAxes() + ggplot2::theme_void() +
    ggplot2::theme(legend.position = "bottom",
                   plot.title = ggplot2::element_blank(),
                   plot.margin = ggplot2::margin(1, 2, 1, 2))
  celltype_plot <- Seurat::DimPlot(obj, group.by = "IDs",
                                   reduction = reduction,
                                   cols = tissue_colors) +
    Seurat::NoAxes() + ggplot2::theme_void() +
    ggplot2::theme(legend.position = "bottom",
                   plot.title = ggplot2::element_blank(),
                   plot.margin = ggplot2::margin(1, 2, 1, 2))
  batch_plot <- Seurat::DimPlot(obj, group.by = "orig.ident",
                                reduction = reduction) +
    Seurat::NoAxes() + ggplot2::theme_void() +
    ggplot2::theme(legend.position = "bottom",
                   plot.title = ggplot2::element_blank(),
                   plot.margin = ggplot2::margin(1, 2, 1, 2))
  method_plot <- celltype_plot + species_plot + batch_plot +
    patchwork::plot_layout(ncol = 3, width = c(1, 1, 1)) +
    patchwork::plot_annotation(title = plotname,  tag_levels = NULL)
  rm(batch_plot, celltype_plot, species_plot, obj)
  gc()
  method_plot
}
umap_list <- mapply(create_umap, filenames, plotnames)

# Reorder the List for nice overview
order <- c(2, 3, 1, 4, 5, 6)
umap_list_ordered <- umap_list[order]
plotnames_ordered <- plotnames[order]

# -----------------------------------------------------------------------------
# Plot the celltype UMAP
celltype_plot_list <- lapply(umap_list_ordered, function(umaps) umaps[[1]])
celltype_plot_list_mod <- mapply(function(plot, name) {
  plot + ggplot2::ggtitle(name) +
    ggplot2::theme(plot.title = element_text(hjust = 0.5)) +
    ggplot2::guides(color = guide_legend(nrow = 4, bycol = TRUE))
}, celltype_plot_list, plotnames_ordered, SIMPLIFY = FALSE)
celltype_umap <- wrap_plots(celltype_plot_list_mod, nrow = 2, ncol = 3) +
  plot_layout(guides = "collect", tag_level = NULL) &
  theme(legend.position = "bottom", panel.border = ggplot2::element_blank())

# -----------------------------------------------------------------------------
# Plot the species UMAP
species_plot_list <- lapply(umap_list_ordered, function(umaps) umaps[[2]])
species <- unique(unlist(lapply(species_plot_list,
                                function(p) levels(p$data$species))))
species_plot_list_mod <- mapply(function(plot, name) {
  plot$data$species <- factor(plot$data$species, levels = species)
  plot + ggplot2::ggtitle(name) +
    ggplot2::theme(plot.title = element_text(hjust = 0.5)) +
    ggplot2::guides(color = guide_legend(nrow = 1, bycol = TRUE))
}, species_plot_list, plotnames_ordered, SIMPLIFY = FALSE)
species_umap <- wrap_plots(species_plot_list_mod, nrow = 2, ncol = 3) +
  plot_layout(guides = "collect", tag_level = NULL) &
  theme(legend.position = "bottom", panel.border = ggplot2::element_blank())

# -----------------------------------------------------------------------------
# Plot the batch UMAP
batch_plot_list <- lapply(umap_list_ordered, function(umaps) umaps[[3]])
batches <- unique(unlist(lapply(batch_plot_list,
                                function(p) levels(p$data$orig.ident))))
batch_map <- setNames(seq_along(batches), batches)
batch_df <- data.frame(orig.ident = names(batch_map),
                       UMAP.id = unname(batch_map))
colnames(batch_df) <- c("orig.ident", "UMAP.id")
filename <- file.path("results", "benchmark", "batch_map.csv")
write.csv(batch_df, filename, row.names = FALSE)
batch_plot_list_mod <- mapply(function(plot, name) {
  plot$data$orig.ident <- factor(batch_map[as.character(plot$data$orig.ident)],
                                 levels = seq_along(batch_map))
  plot + ggplot2::ggtitle(name) +
    ggplot2::theme(plot.title = element_text(hjust = 0.5)) +
    ggplot2::guides(color = guide_legend(nrow = 3, bycol = TRUE))
}, batch_plot_list, plotnames_ordered, SIMPLIFY = FALSE)
batch_umap <- wrap_plots(batch_plot_list_mod, nrow = 2, ncol = 3) +
  plot_layout(guides = "collect", tag_level = NULL) &
  theme(legend.position = "bottom", panel.border = ggplot2::element_blank())

# -----------------------------------------------------------------------------
# Merge and save result
umap_result <- patchwork::wrap_elements(species_umap) /
  patchwork::wrap_elements(celltype_umap) /
  patchwork::wrap_elements(batch_umap) +
  patchwork::plot_layout(heights = c(1, 1.1, 1)) +
  patchwork::plot_annotation(tag_levels = "A") &
  ggplot2::theme(plot.tag = element_text(size = 20, face = "bold"))

# Save
filename <- file.path("results", "benchmark", "UMAP.png")
ggplot2::ggsave(filename, umap_result, width = 12, height = 14)