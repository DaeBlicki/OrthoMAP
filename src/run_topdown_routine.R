#' ----------------------------------------------------------------------------
#' @title run_topdown_routine.R
#'
#' @description
#' Runs OrthoMAP Top-Down workflow. It uses coarse (low resolution) clustering
#' and then fine clustering on the subset.
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
#' @copyright Technau Group, University of Vienna, 2025
#' @copyright Boeva Lab, ETH Zurich, 2025
#'
#' @export
#' ----------------------------------------------------------------------------
run_topdown_routine <- function(orthomap_obj, config, verbose) {
  if (verbose) step_2_enter("top_down")
  # Pre-process the orthomap object
  orthomap_obj <- preprocessing_orthomap_object(
    orthomap_obj = orthomap_obj,
    config = config,
    verbose = verbose
  )
  coarse_clustering_routine <- 1
  fine_clustering_routine <- 2
  # Step 2.1 Coarse clustering
  config_routine <- config$orthomap_top_down$coarse
  if (verbose) {
    cat("----------------------------------------------------------------\n")
    cat("Step 2.1 : Coarse Clustering\n")
    cat("started at ", date(), "\n")
    cat("----------------------------------------------------------------\n\n")
  }
  orthomap_obj <- run_topdown_routine_helper(orthomap_obj = orthomap_obj,
                                             config_routine = config_routine,
                                             step = coarse_clustering_routine,
                                             verbose = verbose)
  orthomap_obj$coarse_clusters <- orthomap_obj$seurat_clusters
  clusters <- levels(orthomap_obj)
  if (verbose) cat("Coarse clusters detected : ", length(clusters), "\n")
  if (verbose) cat("Coarse clustering finished!\n\n")
  #----------------------------------------------------------------------------
  # Step 2.2: Fine Clustering
  #----------------------------------------------------------------------------
  config_routine <- config$orthomap_top_down$fine
  if (verbose) {
    cat("----------------------------------------------------------------\n")
    cat("Step 2.2 : Fine Clustering\n")
    cat("started at ", date(), "\n")
    cat("----------------------------------------------------------------\n\n")
  }
  fine_clusters_list <- vector("list", length(clusters))
  for (i in seq_along(clusters)) {
    # Run clustering on subset
    clust <- clusters[i]
    if (verbose) cat("Current Cluster: ", as.character(clust), "\n")
    if (verbose) cat("started at ", date(), "\n\n")
    orthomap_cluster <- subset(x = orthomap_obj, idents = clust)
    orthomap_cluster <- run_topdown_routine_helper(
      orthomap_obj = orthomap_cluster,
      config_routine = config_routine,
      step = fine_clustering_routine,
      verbose = verbose
    )
    fine_clusters_list[[i]] <- orthomap_cluster$seurat_clusters
    if (verbose) cat("Finished Cluster ", as.character(clust), "\n\n")
    # clear memory
    rm(orthomap_cluster)
    gc()
  }
  fine_clusters <- unlist(fine_clusters_list)
  orthomap_obj$fine_clusters <- fine_clusters
  # map to orthomap clusters
  if (verbose) cat("[OrthoMAP clusters] started at ", date(), "\n")
  num_fine_clust_in_coarse <- unlist(lapply(clusters, function(coarse_cluster) {
    coarse_idx <- which(orthomap_obj$coarse_clusters == coarse_cluster)
    length(unique(orthomap_obj$fine_clusters[coarse_idx]))
  }))
  fine_clust_indexing <- cumsum(num_fine_clust_in_coarse)
  # -------------------------------------------------------------------------
  #' @title map_orthomap_clusters
  #' @brief map coarse and fine clustering to global OrthoMAP cluster
  #' @param coarse_cluster Integer, cluster idx in coarse clustering
  #' @param fine_cluster Integer, cluster idx in fine clustering
  #' @return Integer, OrthoMAP cluster idx
  # -------------------------------------------------------------------------
  map_orthomap_clusters <- function(coarse_cluster, fine_cluster) {
    cluster_idx <- fine_cluster
    if (coarse_cluster > 1) {
      cluster_idx <- fine_cluster + fine_clust_indexing[coarse_cluster - 1]
    }
    cluster_idx
  }
  orthomap_clusters <- mapply(map_orthomap_clusters,
                              as.integer(orthomap_obj$coarse_clusters),
                              as.integer(orthomap_obj$fine_clusters))
  orthomap_obj$orthomap_clusters <- factor(
    orthomap_clusters,
    levels = unique(orthomap_clusters)
  )
  # Run UMAP
  config_routine <- config$orthomap_top_down$coarse
  if (verbose) cat("[RunUMAP] started at ", date(), "\n")
  sd <- config_routine$reduction.stdev
  harmony_dim <- as.integer(
    which(orthomap_obj@reductions$harmony@stdev > sd)
  )
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
  pca_dim <- as.integer(
    which(orthomap_obj@reductions$pca@stdev > sd)
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
  # exit step 2
  if (verbose) step_2_exit()

  # return orthomap_obj
  orthomap_obj
}
