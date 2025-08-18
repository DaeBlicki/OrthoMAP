#'-----------------------------------------------------------------------------
#' @title   OrthoMAP_validator.R
#' @concept Validate given Input and User configuration
#'
#' @description
#' This library contains error handling and functions that may stop the
#' program.
#'
#' Intended as a sourceable script: load with `source("OrthoMAP_validator.R")`
#'
#' @author  David Blickenstorfer, Technau Group (2025)
#' @created 11/07/2025
#'
#' Version: v1.0.0
#' Last Updated: 16/07/2025
#' Changelog:
#'  - v1.0.0: Initial release
#'
#' @copyright Technau Group, University of Vienna, 2025
#' @copyright Boeva Lab, ETH Zurich, 2025
#' ----------------------------------------------------------------------------

#' ----------------------------------------------------------------------------
#' @name      validate_argument_list
#'
#' @description
#' Checks the argument list and usage of the program, using command:
#' Rscript script/main.R (opt)[f]
#' [f] Pipeline flag, {'x','s','c','v'}
#'    'x': execute everything, all-in-one (full pipeline)
#'    's': start, run only step 1 (create OrthoMAP Seurat Object)
#'    'c': continue, run step 2 (standard Seurat standard workflow)
#'    'p': peek, run visualization of pre-processed OrthoMAP Object
#'    'v': visualize, run visualization of processed OrthoMAP Object
#'
#' @param args R argument list, commandArgs(trailingOnly = TRUE)
#'
#' @return pipeline flag [f]
#'
#' @export
#' ----------------------------------------------------------------------------
validate_argument_list <- function(args) {
  # list of valid pipeline flags as described
  valid_flags <- c("x", "s", "p", "c", "v")
  f <- if (length(args) == 1) args[1] else "x"
  # checks if f is valid and argument list is size 1
  if (!(f %in% valid_flags) || length(args) > 1) {
    message("❌ Invalid usage of main.R (opt)[f]\n")
    message(
      "[f] Pipeline flag\n",
      "   'x': execute everything, all-in-one (full pipeline)\n",
      "   's': start, run only step 1 (create OrthoMAP Seurat Object)\n",
      "   'c': continue, run step 2 (standard Seurat standard workflow)\n",
      "   'v': visualize, run visualization with the given arguments\n"
    )
    stop("Usage: Rscript main.R (opt)[f]\n\n")
  }
  # return pipeline flag [f]
  f
}

#' ----------------------------------------------------------------------------
#' @title validate_configuration
#'
#' @description
#' Checks the Configuration file in config.yaml for the Seurat Workflow. It
#' terminate the program in beforehand.
#'
#' @param config yaml file containing config.yaml
#'
#' @importFrom Seurat UpdateSeuratObject
#'
#' @export
#' ----------------------------------------------------------------------------
validate_configuration <- function(config) {
  # ---------------------------------------------------------------------------
  # Validate the configuration for filtering the cells and orthogroups
  # ---------------------------------------------------------------------------
  preprocessing_config <- config$seurat_preprocessing
  # quantitative properties
  min_cells_in_og <- as.numeric(preprocessing_config$min_cells_in_orthogroup)
  max_cells_in_og <- as.numeric(preprocessing_config$max_cells_in_orthogroup)
  min_ogs_in_cell <- as.numeric(preprocessing_config$min_orthogroups_in_cell)
  max_ogs_in_cell <- as.numeric(preprocessing_config$max_orthogroups_in_cell)
  # qualitative properties
  min_count_in_og <- as.numeric(preprocessing_config$min_count_in_orthogroup)
  max_count_in_og <- as.numeric(preprocessing_config$max_count_in_orthogroup)
  min_count_in_cell <- as.numeric(preprocessing_config$min_count_in_cell)
  max_count_in_cell <- as.numeric(preprocessing_config$max_count_in_cell)
  # convergence criterium
  max_iteration <- as.integer(preprocessing_config$max_iteration)
  # No negative coefficient in filtering the cells and orthogroups
  if (min_cells_in_og < 0 || max_cells_in_og < 0 ||
        min_ogs_in_cell < 0 || max_ogs_in_cell < 0) {
    message("❌ Invalid value in pre-processing OrthoMAP!\n")
    message("Current: negative values in pre-processing\n")
    stop("Usage: all coefficients must be zero or positive!\n\n")
  }
  if (min_count_in_og < 0 || max_count_in_og < 0 ||
        min_count_in_cell < 0 || max_count_in_cell < 0) {
    message("❌ Invalid value in pre-processing OrthoMAP!\n")
    message("Current: negative values in pre-processing\n")
    stop("Usage: all coefficients must be zero or positive!\n\n")
  }
  # Min. value is bigger than max. value but max. value is not ignored
  if (min_cells_in_og > max_cells_in_og &&
        max_cells_in_og != 0) {
    message("❌ Invalid value in pre-processing OrthoMAP!\n")
    message("Current: min.coefficient > max.coefficient!\n")
    stop("Usage: min.coefficient < max.coefficient!\n\n")
  }
  if (min_ogs_in_cell > max_ogs_in_cell &&
        max_ogs_in_cell != 0) {
    message("❌ Invalid value in pre-processing OrthoMAP!\n")
    message("Current: min.coefficient > max.coefficient!\n")
    stop("Usage: min.coefficient < max.coefficient!\n\n")
  }
  if (min_count_in_og > max_count_in_og &&
        max_count_in_og != 0) {
    message("❌ Invalid value in pre-processing OrthoMAP!\n")
    message("Current: min.coefficient > max.coefficient!\n")
    stop("Usage: min.coefficient < max.coefficient!\n\n")
  }
  if (min_count_in_cell > max_count_in_cell &&
        max_count_in_cell != 0) {
    message("❌ Invalid value in pre-processing OrthoMAP!\n")
    message("Current: min.coefficient > max.coefficient!\n")
    stop("Usage: min.coefficient < max.coefficient!\n\n")
  }
  # check iteration
  if (max_iteration < 1) {
    message("❌ Invalid value in pre-processing OrthoMAP!\n")
    stop("Usage: max_iteration >= 1")
  }
  # ---------------------------------------------------------------------------
  # Validate the configuration for Seurat Standard Workflow
  # ---------------------------------------------------------------------------
  seurat_config <- config$seurat_standard
  valid <- c("LogNormalize", "CLR", "RC")
  if (!(seurat_config$normalization.method %in% valid)) {
    message("❌ Invalid parameter in Seurat Standard Workflow!\n")
    message(
      "normalization.method\n",
      "   'LogNormalize'\n",
      "   'CLR'\n",
      "   'RC'\n",
      "see https://satijalab.org/seurat/reference/normalizedata\n"
    )
    stop("Invalid parameter in Seurat Standard Workflow!")
  }
  valid <- c("vst", "mean.var.plot", "dispersion")
  if (!(seurat_config$selection.method %in% valid)) {
    message("❌ Invalid parameter in Seurat Standard Workflow!\n")
    message(
      "selection.method\n",
      "   'vst'\n",
      "   'mean.var.plot'\n",
      "   'dispersion'\n",
      "see https://satijalab.org/seurat/reference/findvariablefeatures\n"
    )
    stop("Invalid parameter in Seurat Standard Workflow!")
  }
  valid <- c("pca", "harmony")
  if (!(seurat_config$reduction.space %in% valid)) {
    message("❌ Invalid parameter in Seurat Standard Workflow!\n")
    message(
      "reduction.space\n",
      "   'pca': Use standard principal components\n",
      "   'harmony': Use harmonized principal components\n",
    )
    stop("Invalid parameter in Seurat Standard Workflow!")
  }
  valid <- c("cosine", "euclidian", "manhatten", "hamming")
  if (!(seurat_config$neighbors.annoy.metric %in% valid)) {
    message("❌ Invalid parameter in Seurat Standard Workflow!\n")
    message(
      "neighbors.annoy.metric\n",
      "   'cosine'\n",
      "   'euclidian'\n",
      "   'manhatten'\n",
      "   'hamming'\n",
      "see https://satijalab.org/seurat/reference/findneighbors\n"
    )
    stop("Invalid parameter in Seurat Standard Workflow!")
  }
  valid <- c(1, 2, 3, 4)
  if (!(seurat_config$cluster.algorithm %in% valid)) {
    message("❌ Invalid parameter in Seurat Standard Workflow!\n")
    message(
      "cluster.algorithm\n",
      "   '1': Original Louvain Algorithm\n",
      "   '2': Louvain Algorithm with Multi-Level Refinement\n",
      "   '3': SLM Algorithm\n",
      "   '4': Leiden Algorithm\n",
      "see https://satijalab.org/seurat/reference/findclusters\n"
    )
    stop("Invalid parameter in Seurat Standard Workflow!")
  }
  if (as.numeric(seurat_config$normalization.scale.factor) <= 0) {
    message("❌ Invalid parameter in Seurat Standard Workflow!\n")
    message("normalization.scale.factor > 0\n")
    stop("Invalid parameter in Seurat Standard Workflow!")
  }
  if (as.numeric(seurat_config$selection.nfeatures) <= 0) {
    message("❌ Invalid parameter in Seurat Standard Workflow!\n")
    message("selection.nfeatures > 0\n")
    stop("Invalid parameter in Seurat Standard Workflow!")
  }
  if (as.numeric(seurat_config$pca.dimensions) <= 0) {
    message("❌ Invalid parameter in Seurat Standard Workflow!\n")
    message("pca.dimensions > 0\n")
    stop("Invalid parameter in Seurat Standard Workflow!")
  }
  if (as.numeric(seurat_config$harmony.dimensions) <= 0) {
    message("❌ Invalid parameter in Seurat Standard Workflow!\n")
    message("harmony.dimensions > 0\n")
    stop("Invalid parameter in Seurat Standard Workflow!")
  }
  if (as.numeric(seurat_config$reduction.stdev) <= 0) {
    message("❌ Invalid parameter in Seurat Standard Workflow!\n")
    message("reduction.stdev > 0\n")
    stop("Invalid parameter in Seurat Standard Workflow!")
  }
  if (as.numeric(seurat_config$neighbors.k) <= 0) {
    message("❌ Invalid parameter in Seurat Standard Workflow!\n")
    message("neighbors.k > 0\n")
    stop("Invalid parameter in Seurat Standard Workflow!")
  }
  if (as.numeric(seurat_config$umap.spread) <= 0) {
    message("❌ Invalid parameter in Seurat Standard Workflow!\n")
    message("umap.spread > 0\n")
    stop("Invalid parameter in Seurat Standard Workflow!")
  }
  if (as.numeric(seurat_config$umap.min.dist) <= 0) {
    message("❌ Invalid parameter in Seurat Standard Workflow!\n")
    message("umap.min.dist > 0\n")
    stop("Invalid parameter in Seurat Standard Workflow!")
  }
  if (as.numeric(seurat_config$umap.local.connectivity) <= 0) {
    message("❌ Invalid parameter in Seurat Standard Workflow!\n")
    message("umap.local.connectivity > 0\n")
    stop("Invalid parameter in Seurat Standard Workflow!")
  }
  if (as.numeric(seurat_config$tsne.perplexity) <= 0) {
    message("❌ Invalid parameter in Seurat Standard Workflow!\n")
    message("tsne.perplexity > 0\n")
    stop("Invalid parameter in Seurat Standard Workflow!")
  }
  # ---------------------------------------------------------------------------
  # Validate the configuration for ggplot2
  # ---------------------------------------------------------------------------
  visual_config <- config$visualization
  # list of all color configuration
  visual_vector <- c(
    visual_config$histogram_color,
    visual_config$scatterplot_color,
    visual_config$vline_color,
    visual_config$hline_color
  )
  is_ggplot_color <- lapply(visual_vector, function(color){
    tryCatch({
      # Try to convert the color to RGB
      grDevices::col2rgb(color)
      TRUE
    }, error = function(e) {
      message(color)
      FALSE
    })
  })
  is_ggplot_color <- unlist(is_ggplot_color)
  if (!any(is_ggplot_color)) {
    message("❌ Invalid parameter in Visualization configuration!\n")
    message("Color can not convert to RGB\n")
    stop("Invalid parameter in Visualization configuration!")
  }
}
