#'-----------------------------------------------------------------------------
#' @title   Produce cell-type distance matrices for each species
#' @concept Main file to generate cell-type distance matrices for given species
#'
#' @description
#' This script is the second part of the standard pipeline. It is structured
#' into two parts:
#'
#'  1.  Use the Orthogroups detected in part 1 of the pipeline and the
#'      experession matrices (genes vs. cells) stored in the seurat file to
#'      generate a single expression matrix (Orthogroups vs. cells) for all
#'      species stored in seurat format.
#'
#'  2.  Analyse the generated orthogroup matrix to identify filter parameters
#'      to ensure reliability (enhance memory and runtime management). Apply
#'      Seurat Standard single-cell RNA sequencing workflow.
#'
#' @author  David Blickenstorfer, Technau Group (2025)
#' @created 24/06/2025
#'
#' Version: v1.0.0
#' Changelog:
#'  - v1.0.0: Initial release
#'
#' @copyright Technau Group, University of Vienna, 2025
#' @copyright Boeva Lab, ETH Zurich, 2025
#' ----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Step 0: Preprocessing and Load Parameters
# -----------------------------------------------------------------------------
# Load OrthoMAP "package"
source("src/OrthoMAP_toolbox.R")    # OrthoMAP utility library
args <- commandArgs(trailingOnly = TRUE)
f <- validate_argument_list(args)

# load input variables and check validation file
config <- yaml::read_yaml(file.path("data", "config.yaml"))
sc_objects_path <- file.path("data", "single_cell_atlas")
datatable <- readr::read_csv(file.path("data", "scsRNA_metadata.csv"))
annotation_table_path <- file.path("data", "annotation_table")
validate_configuration(config)

# User development information
verbose <- as.logical(config$general$verbose)
version <- config$general$version

# Input and Output
orthologous_groups_path <- config$input$ortholog_analysis_input
tmp_result_name <- config$output$orthomap_tmp
result_name <- config$output$orthomap_result
orthomap_obj_path <- config$output$orthomap_seurat_object_path
orthomap_visualization_path <- config$output$orthomap_visualization_path
orthomap_statistical_path <- config$output$orthomap_statistical_path

# Print info and welcome message
sessionInfo()
print_config(config)
welcome_message()
print_orthology_interference(version)

# -------------------------------------------------------------------------
# Step 1: Create Orthogroup vs. Cells Expression Matrix
# -------------------------------------------------------------------------
if (f == "x" || f == "s") {
  if (verbose) step_1_enter()
  # Generate input path depending on the flag and name
  dir.create(orthomap_obj_path,
             recursive = TRUE,
             showWarnings = FALSE)
  orthomap_obj <- create_orthomap_object(
    orthologous_groups_path = orthologous_groups_path,
    datatable = datatable,
    single_cell_directory = sc_objects_path,
    annotation_table_directory = annotation_table_path,
    output_directory = orthomap_obj_path,
    output_name = tmp_result_name,
    verbose = verbose
  )

  # create pre-processed orthomap object directory
  raw_visualization_path <- file.path(orthomap_visualization_path, "raw")
  dir.create(raw_visualization_path,
             recursive = TRUE,
             showWarnings = FALSE)

  # visualize the orthomap object
  if (verbose) cat("[Visualization] started at ", date(), "\n")
  visualize_orthomap_object(
    orthomap_obj = orthomap_obj,
    config = config,
    output_directory = raw_visualization_path,
    verbose = verbose
  )

  # Exit Step 1
  if (verbose) step_1_exit()
}

# -------------------------------------------------------------------------
# Load the OrthoMAP Seurat Object when use "p" or "c"
# -------------------------------------------------------------------------
if (f == "p" || f == "c") {
  if (verbose) cat("[Loading] started at ", date(), "\n\n")
  filename <- file.path(orthomap_obj_path, tmp_result_name)
  orthomap_obj <- get(load(filename))
  # ---------------------------------------------------------------------------
  # Visualization the OrthoMAP Seurat Object when using [f] = "p"
  # ---------------------------------------------------------------------------
  if (f == "p") {
    if (verbose) cat("[Visualize Distribution] started at ", date(), "\n\n")
    # need to load the seurat object
    filename <- file.path(orthomap_obj_path, tmp_result_name)
    orthomap_obj <- get(load(filename))
    raw_visualization_path <- file.path(orthomap_visualization_path, "raw")
    # visualize the orthomap object
    visualize_orthomap_object(
      orthomap_obj = orthomap_obj,
      config = config,
      output_directory = raw_visualization_path,
      verbose = verbose
    )
    if (verbose) cat("Unprocessed OrthoMAP Seurat Object visualized!\n\n")
  }

}
# -------------------------------------------------------------------------
# Step 2: Apply Seurat Pipeline Analysis
# -------------------------------------------------------------------------
if (f == "x" || f == "c") {
  routine <- config$general$workflow
  # Apply Standard routine
  if (routine == "standard") {
    orthomap_obj <- run_standard_routine(orthomap_obj = orthomap_obj,
                                         config = config,
                                         verbose = verbose)
  } else if (routine == "top_down") {
    orthomap_obj <- run_topdown_routine(orthomap_obj = orthomap_obj,
                                        config = config,
                                        verbose = verbose)
  }

  # Save the resulting OrthoMAP Seurat Object
  if (is.character(result_name)) {
    filename <- file.path(orthomap_obj_path, result_name)
    if (verbose) cat("[Save Result] save result as ", filename, "\n\n")
    save(orthomap_obj, file = filename)
  }
}

# -------------------------------------------------------------------------
# Step 3: Visualize and Analyze
# -------------------------------------------------------------------------
if (f == "v") {
  # Load the Orthomap Result Object when necessary
  if (verbose) cat("[Loading] started at ", date(), "\n\n")
  filename <- file.path(orthomap_obj_path, result_name)
  orthomap_obj <- get(load(filename))
}
if (f == "x" || f == "c" || f == "v") {
  if (verbose) step_3_enter()
  result_visualization_path <- file.path(
    orthomap_visualization_path, "result"
  )
  dir.create(result_visualization_path,
             recursive = TRUE,
             showWarnings = FALSE)
  # Step 2: Visualize UMAP and Seurat pipeline relevance
  if (verbose) cat("[Visualize UMAP] started at ", date(), "\n\n")
  visualize_orthomap_result(
    orthomap_obj = orthomap_obj,
    config = config,
    output_directory = result_visualization_path,
    verbose = verbose
  )
  # Step 3: Statistical analysise and empirical estimation
  if (verbose) cat("[Evaluate Result] started at ", date(), "\n")
  dir.create(orthomap_statistical_path,
             recursive = TRUE,
             showWarnings = FALSE)
  orthomap_obj <- analyze_orthomap_result(
    orthomap_obj = orthomap_obj,
    output_directory = orthomap_statistical_path,
    config = config,
    verbose = verbose
  )
  # Apply Standard routine
  if (config$general$workflow == "standard") {
    resolutions <- seq(0.1, 1, by = 0.1)
    for (res in resolutions) {
      if (verbose) cat("Current Resolution: ", res, "\n")
      output_directory <- file.path(orthomap_statistical_path,
                                    as.character(res))
      dir.create(output_directory,
                 recursive = TRUE,
                 showWarnings = FALSE)
      colname <- paste0("RNA_snn_res.", res)
      orthomap_obj$orthomap_clusters <- orthomap_obj[[colname]]
      orthomap_obj <- analyze_orthomap_result(
        orthomap_obj = orthomap_obj,
        output_directory = output_directory,
        config = config,
        verbose = verbose
      )
    }
  }
}

# Exit Message
exit_message()