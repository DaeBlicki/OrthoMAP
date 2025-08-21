#'-----------------------------------------------------------------------------
#' @title   OrthoMAP_print.R
#' @concept R utility file to print information in case of verbose
#'
#' @importFrom yaml as.yaml
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
#' ----------------------------------------------------------------------------
#' @name      print_session
#' @concept   prints the loaded library and configuration
#'
#' @param config config.yaml as RObject
#'
#' @export
#' ----------------------------------------------------------------------------
print_config <- function(config) {
  cat("\n")
  cat("OrthoMAP Configuration:\n")
  cat(yaml::as.yaml(config))
  cat("\n")
}

#' ----------------------------------------------------------------------------
#' @name      welcome_message
#' @concept   print out the welcome message of the pipeline
#'
#' @export
#' ----------------------------------------------------------------------------
welcome_message <- function() {
  cat("------------------------------------------------------------------\n")
  cat("OrthoMAP - Interpecific Cell-Type Tree Construction Method\n")
  cat("Produce cell-type distance matrices for each species\n")
  cat("David Blickenstorfer \n\n")

  cat("(c) Technau Group, University of Vienna, 2025\n")
  cat("(c) Boeva Lab, ETH Zurich, 2025 \n")
  cat("------------------------------------------------------------------\n\n")
}

#' ----------------------------------------------------------------------------
#' @name      exit_message
#' @concept   print out the exit message of the pipeline
#'
#' @export
#' ----------------------------------------------------------------------------
exit_message <- function() {
  cat("------------------------------------------------------------------\n")
  cat("Interpecific Cell-Type Tree Construction Method Pipeline succeeded!\n")
  cat("Finished at ", date(), "\n")
  cat("------------------------------------------------------------------\n")
}

#' ----------------------------------------------------------------------------
#' @name      print_orthology_interference
#' @concept   prints the current input of the analysis
#'
#' @param version, string, version parameter in the config file
#'
#' @export
#' ----------------------------------------------------------------------------
print_orthology_interference <- function(version) {
  cat("------------------------------------------------------------------\n")
  cat("Used Orthology Interference Analysis Tool\n")
  cat("Version: ", version, "\n")
  cat("------------------------------------------------------------------\n\n")
}

#' ----------------------------------------------------------------------------
#' @name      step_1_enter
#' @concept   enter the step: "Create OrthoGroup against Celltype matrix".
#' @param orthogroup_analysis: "oma"          : use OMA input
#'                             "diamond"      : use OrthoFinder Diamond input
#'                             "mmseqs2"      : use OrthoFinder MMseqs2 input
#' @param config config.yaml as RObject
#'
#' @export
#' ----------------------------------------------------------------------------
step_1_enter <- function() {
  # print out the Orthogroup analyser
  cat("------------------------------------------------------------------\n")
  cat("Step 1: Create OrthoGroup against Celltype matrix\n")
  cat("Started at ", date(), "\n")
  cat("------------------------------------------------------------------\n\n")
}

#' ----------------------------------------------------------------------------
#' @name      step_1_exit
#' @concept   exit the step: "Create OrthoGroup against Celltype matrix".
#'
#' @export
#' ----------------------------------------------------------------------------
step_1_exit <- function() {
  cat("OrthoGroup against Celltype matrix for all species created\n")
  cat("Finished at ", date(), "\n")
  cat("------------------------------------------------------------------\n\n")
}

#' ----------------------------------------------------------------------------
#' @name      step_3_enter
#' @concept   enter the step: "Create OrthoGroup against Celltype matrix".
#' @param version: prints out the version
#' @param config config.yaml as RObject
#'
#' @export
#' ----------------------------------------------------------------------------
step_3_enter <- function() {
  cat("------------------------------------------------------------------\n")
  cat("Step 3: Visualization and Analyse \n")
  cat("Started at ", date(), "\n")
  cat("------------------------------------------------------------------\n\n")
}

#' ----------------------------------------------------------------------------
#' @name      visualize_exit
#' @concept   exit the visualization.
#'
#' @export
#' ----------------------------------------------------------------------------
step_3_exit <- function() {
  cat("OrthoMAP Seurat Object visualized and analyzed!\n")
  cat("Finished at ", date(), "\n")
  cat("------------------------------------------------------------------\n\n")
}

#' ----------------------------------------------------------------------------
#' @name      step_2_enter
#' @concept   enter the step: "Analyse OG Expression and Apply Transformation".
#' @param pipeline_approach: "standard"   : Seurat standard approach
#'                           "bottom_up"  : Fine filtering into merging
#'                           "top_down"   : Coarse filtering into clustering
#'
#' @export
#' ----------------------------------------------------------------------------
step_2_enter <- function(pipeline_approach) {
  # print out the Orthogroup analyser
  cat("------------------------------------------------------------------\n")
  cat("Step 2: Apply Seurat Object Standard Workflow\n")
  if (pipeline_approach == "standard") {
    cat("Method: Seurat Standard Pipeline\n")
  } else if (pipeline_approach == "bottom_up") {
    cat("Method: OrthoMAP Bottom-Up\n")
  } else if (pipeline_approach == "top_down") {
    cat("Method: OrthoMAP Top-Down\n")
  }
  cat("Started at ", date(), "\n")
  cat("------------------------------------------------------------------\n\n")
}

#' ----------------------------------------------------------------------------
#' @name      step_2_exit
#' @concept   exit the step: "Analyse OG Expression and Apply Transformation".
#'
#' @export
#' ----------------------------------------------------------------------------
step_2_exit <- function() {
  cat("OrthoMAP Seurat Object updated!\n")
  cat("Finished at ", date(), "\n")
  cat("------------------------------------------------------------------\n\n")
}

#' ----------------------------------------------------------------------------
#' @name      print_preprocessing_config
#' @concept   prints the OrthoMAP preprocessing parameter
#'
#' @param config config.yaml as RObject
#' @export
#' ----------------------------------------------------------------------------
print_preprocessing_config <- function(config) {
  # Access to information
  preprocessing_config <- config$seurat_preprocessing
  # quantitative properties
  min_cell_in_og <- as.numeric(preprocessing_config$min_cells_in_orthogroup)
  max_cell_in_og <- as.numeric(preprocessing_config$max_cells_in_orthogroup)
  min_og_in_cell <- as.numeric(preprocessing_config$min_orthogroups_in_cell)
  max_og_in_cell <- as.numeric(preprocessing_config$max_orthogroups_in_cell)
  # qualitative properties
  min_count_in_og <- as.numeric(preprocessing_config$min_count_in_orthogroup)
  max_count_in_og <- as.numeric(preprocessing_config$max_count_in_orthogroup)
  min_count_in_cell <- as.numeric(preprocessing_config$min_count_in_cell)
  max_count_in_cell <- as.numeric(preprocessing_config$max_count_in_cell)
  # convergence criterium
  max_iter <- as.integer(preprocessing_config$max_iteration)
  # Print information
  cat("Start preprocessing OrthoMAP Seurat Object!\n")
  cat("Begin at ", date(), "\n\n")
  cat(
    "Maximal number of iterations: ", max_iter, "\n\n",
    "Remove orthogroups with less than ", min_cell_in_og, "cells!\n",
    "Remove orthogroups with more than ", max_cell_in_og, "cells!\n",
    "Remove cells with less than ", min_og_in_cell, "orthogroups!\n",
    "Remove cells with more than ", max_og_in_cell, "orthogroups!\n\n",
    "Remove orthogroups with less than ", min_count_in_og, "reads!\n",
    "Remove orthogroups with more than ", max_count_in_og, "reads!\n",
    "Remove cells with less than ", min_count_in_cell, "reads!\n",
    "Remove cells with more than ", max_count_in_cell, "reads!\n\n",
    "\n"
  )
}
