#' ----------------------------------------------------------------------------
#' @title visualize_orthomap_object
#'
#' @description
#' Estimate the number of nCount and nFeature for orthogroup and cells. It
#' generates statistical summary and distribution plots for the OrthoMAP
#' Seurat Object.
#'
#' @details
#' Distribution of Cells
#'    Histogram: Plot distribution of the cells in term of nFeature and nCount
#'    Log-Log plot: Double logarithmic plot for nFeature and nCount
#'    Scatter-Plot: nFeature against nCount plot for each cell
#' Distribution of Orthogroup
#'    Histogram: Plot distribution of the OG in term of nFeature and nCount
#'    Log-Log plot: Double logarithmic plot for nFeature and nCount
#'    Scatter-Plot: nFeature against nCount plot for each Orthogroup
#' Summary Table: A statistical table including metric variables
#'
#' @param orthomap_obj The Seurat Object generated from step 1 of the pipeline.
#' @param config yaml configuration file
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
#'     seurat_configuration = seurat_config
#'     output_directory = method_orthomap_qc_path,
#'     output_format = ".png",
#'     verbose = verbose,
#'   )
#' }
#'
#' @author  David Blickenstorfer, Technau Group (2025)
#' @created 07/07/2025
#'
#' @copyright Technau Group, University of Vienna, 2025
#' @copyright Boeva Lab, ETH Zurich, 2025
#'
#' @export
#' ----------------------------------------------------------------------------
visualize_orthomap_object <- function(orthomap_obj,
                                      config,
                                      output_directory,
                                      verbose = TRUE) {
  # Inform the user
  if (verbose) {
    cat("Start visualize distribution of OrthoMAP variables\n")
    cat("Begin at ", date(), "\n\n")
  }
  # ===========================================================================
  # Step 0: Get Configuration
  # ===========================================================================
  output_format <- config$visualization$result_format
  histogram_color <- config$visualization$histogram_color
  scatterplot_color <- config$visualization$scatterplot_color
  vline_color <- config$visualization$vline_color
  # quantitative
  og_low <- config$seurat_preprocessing$min_cells_in_orthogroup
  og_high <- config$seurat_preprocessing$max_cells_in_orthogroup
  cell_low <- config$seurat_preprocessing$min_orthogroups_in_cell
  cell_high <- config$seurat_preprocessing$max_orthogroups_in_cell
  #qualitative
  og_count_low <- config$seurat_preprocessing$min_count_in_orthogroup
  og_count_high <- config$seurat_preprocessing$max_count_in_orthogroup
  cell_count_low <- config$seurat_preprocessing$min_count_in_cell
  cell_count_high <- config$seurat_preprocessing$max_count_in_cell

  # ===========================================================================
  # Step 1: Get Data from OrthoMAP Seurat Object
  # ===========================================================================
  og_counts <- Seurat::GetAssayData(orthomap_obj, layer = "counts")
  og_counts <- as(og_counts, "CsparseMatrix")
  nog_in_cell <- diff(og_counts@p)
  ncount_in_cell <- Matrix::colSums(og_counts)
  cell_names <- colnames(og_counts)
  df_cell <- data.frame(Cells = cell_names,
                        nOG = nog_in_cell,
                        nCount = ncount_in_cell)

  # get information about orthogroups
  og_counts <- as(og_counts, "RsparseMatrix")
  ncell_in_og <- diff(og_counts@p)
  ncount_in_og <- Matrix::rowSums(og_counts)
  og_names <- rownames(og_counts)
  df_og <- data.frame(Orthogroup = og_names,
                      nCell = ncell_in_og,
                      nCount = ncount_in_og)

  # ===========================================================================
  # Step 2: Print Summary and Statistical Table
  # ===========================================================================
  dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
  data_list <- list(nog_in_cell, ncount_in_cell, ncell_in_og, ncount_in_og)
  data_names <- c("nOrthogroup_in_cell", "nCount_in_cell",
                  "nCell_in_orthogroup", "nCount_in_orthogroup")
  generate_summary_table(data_list = data_list,
                         data_names = data_names,
                         output_directory = output_directory,
                         verbose = verbose)

  # ===========================================================================
  # Step 3: Plot Results
  # ===========================================================================
  # create directories for pdf
  dir.create(file.path(output_directory, "Histogram"),
             recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(output_directory, "Logscale"),
             recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(output_directory, "Scatter"),
             recursive = TRUE, showWarnings = FALSE)
  if (verbose) {
    cat("Generate visualisation results: Histogram and Logscale\n")
  }
  # colors


  # nFeature against cell distribution
  vlines_og <- add_vlines(og_low, og_high, vline_color)
  vline_count_og <- add_vlines(og_count_low, og_count_high, vline_color)
  vlines_cell <- add_vlines(cell_low, cell_high, vline_color)
  vline_count_cell <- add_vlines(cell_count_low, cell_count_high, vline_color)
  histogram_og_in_cell <- ggplot2::ggplot(df_cell, aes(x = nOG)) +
    ggplot2::geom_histogram(binwidth = 1,
                            fill = histogram_color,
                            color = histogram_color) +
    vlines_cell +
    ggplot2::labs(x = "Orthogroup per cell", y = "Frequency") +
    ggplot2::theme_bw()
  filename <- file.path(output_directory, "Histogram", "nOG_in_cell")
  filename <- paste0(filename, output_format)
  ggplot2::ggsave(filename, histogram_og_in_cell, width = 6, height = 4)

  histogram_og_in_cell <- ggplot2::ggplot(df_cell, aes(x = nOG)) +
    ggplot2::geom_histogram(binwidth = 0.01,
                            fill = histogram_color,
                            color = histogram_color) +
    vlines_cell +
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    scale_y_continuous(trans = pseudo_log_trans(base = 10),
                       breaks = trans_breaks("log10", function(x) 10^x),
                       labels = trans_format("log10", math_format(10^.x))) +
    annotation_logticks() + expand_limits(x = 0.1, y = 0.1) +
    ggplot2::labs(x = "Orthogroup per cell", y = "Frequency") +
    ggplot2::theme_bw()
  filename <- file.path(output_directory, "Logscale", "nOG_in_cell")
  filename <- paste0(filename, output_format)
  ggplot2::ggsave(filename, histogram_og_in_cell, width = 6, height = 4)

  # nCount against cell distribution
  histogram_count_in_cell <- ggplot(df_cell, aes(x = nCount)) +
    ggplot2::geom_histogram(binwidth = 10,
                            fill = histogram_color,
                            color = histogram_color) +
    vline_count_cell +
    ggplot2::labs(x = "Count per cell", y = "Frequency") +
    ggplot2::theme_bw()
  filename <- file.path(output_directory, "Histogram", "nCount_in_cell")
  filename <- paste0(filename, output_format)
  ggplot2::ggsave(filename, histogram_count_in_cell, width = 6, height = 4)

  histogram_count_in_cell <- ggplot(df_cell, aes(x = nCount)) +
    ggplot2::geom_histogram(binwidth = 0.01,
                            fill = histogram_color,
                            color = histogram_color) +
    vline_count_cell +
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    scale_y_continuous(trans = pseudo_log_trans(base = 10),
                       breaks = trans_breaks("log10", function(x) 10^x),
                       labels = trans_format("log10", math_format(10^.x))) +
    annotation_logticks() + expand_limits(x = 0.1, y = 0.1) +
    ggplot2::labs(x = "Count per cell", y = "Frequency") +
    ggplot2::theme_bw()
  filename <- file.path(output_directory, "Logscale", "nCount_in_cell")
  filename <- paste0(filename, output_format)
  ggplot2::ggsave(filename, histogram_count_in_cell, width = 6, height = 4)

  # nFeature against orthogroup distribution
  histogram_cell_in_og <- ggplot2::ggplot(df_og, aes(x = nCell)) +
    ggplot2::geom_histogram(binwidth = 10,
                            fill = histogram_color,
                            color = histogram_color) +
    vlines_og +
    labs(x = "Cell per orthogroup", y = "Frequency") +
    theme_bw()
  filename <- file.path(output_directory, "Histogram", "nCell_in_OG")
  filename <- paste0(filename, output_format)
  ggsave(filename, histogram_cell_in_og, width = 6, height = 4)
  histogram_cell_in_og <- ggplot2::ggplot(df_og, aes(x = nCell)) +
    ggplot2::geom_histogram(binwidth = 0.01,
                            fill = histogram_color,
                            color = histogram_color) +
    vlines_og +
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    scale_y_continuous(trans = pseudo_log_trans(base = 10),
                       breaks = trans_breaks("log10", function(x) 10^x),
                       labels = trans_format("log10", math_format(10^.x))) +
    annotation_logticks() + expand_limits(x = 0.1, y = 0.1) +
    labs(x = "Cell per orthogroup", y = "Frequency") +
    theme_bw()
  filename <- file.path(output_directory, "Logscale", "nCell_in_OG")
  filename <- paste0(filename, output_format)
  ggsave(filename, histogram_cell_in_og, width = 6, height = 4)

  # nFeature against orthogroup distribution
  histogram_ncount_in_og <- ggplot2::ggplot(df_og, aes(x = nCount)) +
    ggplot2::geom_histogram(binwidth = 100,
                            fill = histogram_color,
                            color = histogram_color) +
    vline_count_og +
    labs(x = "Count per orthogroup", y = "Frequency") +
    theme_bw()
  filename <- file.path(output_directory, "Histogram", "nCount_in_OG")
  filename <- paste0(filename, output_format)
  ggsave(filename, histogram_ncount_in_og, width = 6, height = 4)

  histogram_ncount_in_og <- ggplot2::ggplot(df_og, aes(x = nCount)) +
    ggplot2::geom_histogram(binwidth = 0.01,
                            fill = histogram_color,
                            color = histogram_color) +
    vline_count_og +
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    scale_y_continuous(trans = pseudo_log_trans(base = 10),
                       breaks = trans_breaks("log10", function(x) 10^x),
                       labels = trans_format("log10", math_format(10^.x))) +
    annotation_logticks() + expand_limits(x = 0.1, y = 0.1) +
    labs(x = "Count per orthogroup", y = "Frequency") +
    theme_bw()
  filename <- file.path(output_directory, "Logscale", "nCount_in_OG")
  filename <- paste0(filename, output_format)
  ggsave(filename, histogram_ncount_in_og, width = 6, height = 4)

  if (verbose) cat("Generate Feature-Feature scatter plot\n")
  # Create Scatterplot: Cell vs Count
  og_scatter_plot <- ggplot(df_og, aes(y = nCount, x = nCell)) +
    geom_point(color = scatterplot_color, size = 0.5, alpha = 0.7) +
    labs(x = "Orthogroup", y = "Count") +
    scale_x_log10() + scale_y_log10() +
    annotation_logticks() + expand_limits(x = 0.1, y = 0.1) +
    theme_bw()
  filename <- file.path(output_directory, "Scatter", "OG_vs_Count")
  filename <- paste0(filename, output_format)
  ggsave(filename, og_scatter_plot, width = 6, height = 4)

  # Create Scattterplot: OG vs Count
  cell_scatter_plot <- ggplot(df_cell, aes(y = nCount, x = nOG)) +
    geom_point(color = scatterplot_color, size = 0.5, alpha = 0.7) +
    labs(x = "Cell", y = "Count") +
    scale_x_log10() + scale_y_log10() +
    annotation_logticks() + expand_limits(x = 0.1, y = 0.1) +
    theme_bw()
  filename <- file.path(output_directory, "Scatter", "Cell_vs_Count")
  filename <- paste0(filename, output_format)
  ggsave(filename, cell_scatter_plot, width = 6, height = 4)
  if (verbose) {
    cat("OrthoMAP Seurat Object successful visualized!\n")
    cat("Finished at ", date(), "\n\n")
  }
}