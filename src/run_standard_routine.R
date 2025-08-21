#' ----------------------------------------------------------------------------
#' @title run_standard_routine.R
#'
#' @description
#' Runs Seurat standard pipeline modified for cross-species comparison
#'
#' @details
#' Short description of the workflow. Uses the parameter defined in config
#' 2.1 : Pre-processing workflow
#' 2.2 : Normalizing the data
#' 2.2 : Identify and merge of highly variable features for each species
#' 2.3 : Scaling the data
#' 2.4 : Perform linear dimensional reduction
#' 2.5 : Find k-nearest neighbors
#' 2.6 : Find cluster
#' 2.7 : Run UMAP
#'
#' @param orthomap_obj OrthoMAP Seurat Object
#' @param config yaml configuration file
#' @param verbose Logical; if TRUE, prints progress messages.
#'
#' @return Analyzed OrthoMAP Seurat Object
#'
#' @importFrom Seurat NormalizeData
#' @importFrom Seurat GetAssayData
#'
#' @seealso \code{\link[Seurat]{CreateSeuratObject}}
#' @references \url{https://cran.r-project.org/web/packages/Seurat/index.html}
#'
#' @examples
#' \dontrun{
#' orthomap_obj <- run_standard_routine(
#'     orthomap_obj = orthomap_obj,
#'     config = config,
#'     verbose = TRUE
#' )
#' }
#'
#' @author  David Blickenstorfer, Technau Group (2025)
#' @created 31/07/2025
#'
#' Version: v1.0.0
#' Changelog:
#'  - v1.0.0: Initial release
#'
#' @copyright Technau Group, University of Vienna, 2025
#' @copyright Boeva Lab, ETH Zurich, 2025
#'
#' @export
#' ----------------------------------------------------------------------------
run_standard_routine <- function(orthomap_obj, config, verbose) {
  if (verbose) step_2_enter("standard")
  config_routine <- config$seurat_standard
  # Step 2.0: preprocess the orthomap object
  orthomap_obj <- preprocessing_orthomap_object(
    orthomap_obj = orthomap_obj,
    config = config,
    verbose = verbose
  )

  # Step 2.1 Log-Normalize
  if (verbose) cat("[Log-Normalize] started at ", date(), "\n")
  orthomap_obj <- Seurat::NormalizeData(
    object = orthomap_obj,
    normalization.method = config_routine$normalization.method,
    scale.factor = config_routine$normalization.scale.factor,
    verbose = verbose
  )
  if (verbose) cat("Finished normalization!\n\n")

  # Step 2.2 Feature-Selection
  # Choose most variable features per species then concartenate
  if (verbose) cat("[Feature selection] started at ", date(), "\n")
  species <- unique(as.character(orthomap_obj$species))
  most_variable_features <- lapply(species, function(sp) {
    if (verbose) {
      cat("Find most variable features for ", sp, "\n")
      cat("Started at ", date(), "\n")
    }
    species_obj <- subset(orthomap_obj, subset = species == sp)
    species_obj <- Seurat::FindVariableFeatures(
      object = species_obj,
      selection.method = config_routine$selection.method,
      nfeatures = config_routine$selection.nfeatures,
      verbose = verbose
    )
    Seurat::VariableFeatures(species_obj)
  })
  # Get get concartenated most variable feature list as unique vector
  most_variable_features <- unlist(most_variable_features)
  merged_variable_features <- unique(most_variable_features)
  Seurat::VariableFeatures(orthomap_obj) <- merged_variable_features
  if (verbose) {
    cat("Number of VariableFeatures : ",
        length(Seurat::VariableFeatures(orthomap_obj)), "\n\n")
  }

  # Step 2.3 Scale the data
  if (verbose) cat("[Data scaling] started at ", date(), "\n")
  orthomap_obj <- Seurat::ScaleData(
    object = orthomap_obj,
    split.by = config_routine$scale.features,
    verbose = verbose
  )
  if (verbose) cat("Data scaling finished!\n\n")

  # Step 2.4 Run pca and harmony
  sd <- config_routine$reduction.stdev
  # Run PCA
  if (verbose) cat("[RunPCA] started at ", date(), "\n")
  orthomap_obj <- Seurat::RunPCA(
    object = orthomap_obj,
    pcs.compute = config_routine$pca.dimensions,
    verbose = verbose
  )
  pca_dim <- as.integer(
    which(orthomap_obj@reductions$pca@stdev > sd)
  )
  if (verbose) {
    cat("Number of PCA dimensions : ", length(pca_dim), "\n\n")
  }

  # Run Harmony
  if (verbose) cat("[RunHarmony] started at ", date(), "\n")
  orthomap_obj <- harmony::RunHarmony(
    object = orthomap_obj,
    group.by.vars = config_routine$reduction.group.by,
    verbose = verbose
  )
  harmony_dim <- as.integer(
    which(orthomap_obj@reductions$harmony@stdev > sd)
  )
  if (verbose) {
    cat("Number of harmony dimensions : ", length(harmony_dim), "\n\n")
  }

  # Step 2.5 Find k-nearest neighbors
  if (verbose) cat("[FindNeighbors] started at ", date(), "\n")
  reduction <- config_routine$reduction.space
  dims <- pca_dim
  if (reduction == "harmony") dims <- harmony_dim
  orthomap_obj <- Seurat::FindNeighbors(
    object = orthomap_obj,
    dims = dims,
    reduction = config_routine$reduction.space,
    annoy.metric = config_routine$neighbors.annoy.metric,
    k.param = config_routine$neighbors.k,
    verbose = verbose
  )
  if (verbose) cat("Finished kNN algorithm!\n\n")

  # Step 2.6 Find cluster, iterate until convergence
  if (verbose) cat("[Build Cluster] started at ", date(), "\n")
  resolutions <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
  for (cur_resolution in resolutions) {
    if (verbose) cat("Current Resolution : ", cur_resolution, "\n")
    orthomap_obj <- Seurat::FindClusters(
      object = orthomap_obj,
      resolution = cur_resolution,
      random.seed = config$seurat_standard$cluster.seed,
      algorithm = config$seurat_standard$cluster.algorithm,
      verbose = verbose
    )
  }
  orthomap_obj$orthomap_clusters <- factor(
    orthomap_obj$seurat_clusters,
    levels = unique(orthomap_obj$seurat_clusters)
  )
  if (verbose) cat("Finished cluster building!\n\n")

  # Step 2.7 Run UMAP for harmony and pca
  # Run UMAP
  if (verbose) cat("[RunUMAP] started at ", date(), "\n")
  orthomap_obj <- Seurat::RunUMAP(
    orthomap_obj,
    dims = harmony_dim,
    reduction = "harmony",
    reduction.name = "umap.harmony", reduction.key = "umap.harmony",
    n.neighbors = config_routine$neighbors.k,
    spread = config_routine$umap.spread,
    min.dist = config_routine$umap.min.dist,
    local.connectivity = config_routine$umap.local.connectivity,
    verbose = verbose
  )
  orthomap_obj <- Seurat::RunUMAP(
    orthomap_obj,
    dims = pca_dim,
    reduction = "pca",
    reduction.name = "umap.pca", reduction.key = "umap.pca",
    n.neighbors = config_routine$neighbors.k,
    spread = config_routine$umap.spread,
    min.dist = config_routine$umap.min.dist,
    local.connectivity = config_routine$umap.local.connectivity,
    verbose = verbose
  )
  if (verbose) cat("Finished UMAP\n\n")

  # exit step 2
  if (verbose) step_2_exit()

  # return orthomap
  orthomap_obj
}
