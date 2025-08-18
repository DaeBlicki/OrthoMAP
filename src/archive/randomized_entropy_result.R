#'-----------------------------------------------------------------------------
#' @title   randomized_entropy_result.R
#'
#' @concept
#' Creates an randomized clustering and calculated it harmonized Shannon
#' Entropy. It serves as control group. The number of cells and number of
#' cell tissue can be given as parameter.
#'
#' @param num_cells Integer, total number of cells in sample
#' @param num_tissues Integer, number of cell tissue in sample
#' @param num_clusters Integer, number of identified cluster
#' @param rng_seed Integer, randoms seed number for reproducibility
#'
#' @return Numeric Vector, contains harmonized Shannon Entropy for meta data
#' that randomely clusters the samples.
#'
#' @author  David Blickenstorfer, Technau Group (2025)
#' @created 05/08/2025
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
randomized_entropy_result <- function(num_cells = 150000,
                                      num_tissues  = 50,
                                      num_clusters = 50,
                                      rng_seed = 42) {
  ## create randomized vectors
  # tissue samples from exponential distribution
  # cluster samples from uniform distribution
  set.seed(rng_seed)
  tissue_vals <- rexp(num_cells, rate = 1)
  tissue_norm <- tissue_vals / max(tissue_vals)
  tissue_label <- ceiling(tissue_norm * num_tissues)
  tissue_label[tissue_label > num_tissues] <- num_tissues
  tissue_label <- paste0("CT", tissue_label)
  cluster_label <- sample(1:num_clusters, size = num_cells, replace = TRUE)
  cluster_label <- paste0("Cluster-", cluster_label)
  raw_table <- as.data.frame.matrix(table(cluster_label, tissue_label))
  sum_table <- apply(raw_table, 1, sum)
  pct_table <- raw_table / sum_table
  # Estimate the Shannon Entropy for each cluster
  entropy_vector <- apply(pct_table, 1, function(pct_table_line) {
    entropy::entropy.empirical(pct_table_line, unit = "log")
  })
  harmonized_entropy_vector <- entropy_vector / log(sum_table)
  harmonized_entropy_vector[sum_table == 1] <- 0.0
  harmonized_entropy_vector
}

#'-----------------------------------------------------------------------------
#' @title   simulate_entropy_result.R
#'
#' @concept
#' Use the given OrthoMAP Object to calculate the exact tissue and cluster
#' percentage. It uniformly and randomely distributes the tissue and cluster
#' assignement and then estimates the shannon entropy. It's used to simulate
#' the shannon entropy of the OrthoMAP Object when the clustering is random.
#'
#' @param orthomap_obj OrthoMAP Object
#' @param rng_seed Integer, randoms seed number for reproducibility
#'
#' @return Numeric Vector, contains harmonized Shannon Entropy for meta data
#' that randomely clusters the samples.
#'
#' @author  David Blickenstorfer, Technau Group (2025)
#' @created 05/08/2025
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
simulate_entropy_result <- function(orthomap_obj, rng_seed = 42) {
  # estimate tissue type percentage
  set.seed(rng_seed)
  ct_names <- unique(orthomap_obj$IDs)
  num_ct <- length(ct_names)
  ct_size_table <- unlist(lapply(ct_names, function(ct) {
    length(which(orthomap_obj$IDs == ct))
  }))
  num_cells <- sum(ct_size_table)
  ct_pct_table <- ct_size_table / num_cells
  ct_label <- sample(1:num_ct, size = num_cells,
                     replace = TRUE, prob = ct_pct_table)
  ct_label <- paste0("CT-", ct_label)
  # estimate clustering percentage
  cluster_names <- unique(orthomap_obj$orthomap_clusters)
  num_cluster <- length(cluster_names)
  clust_size_table <- unlist(lapply(cluster_names, function(clust) {
    length(which(orthomap_obj$orthomap_clusters == clust))
  }))
  clust_pct_table <- clust_size_table / num_cells
  clust_label <- sample(1:num_cluster, size = num_cells,
                        replace = TRUE, prob = clust_pct_table)
  # evaluate Shannon Entropy
  rand_ct_raw_table <- as.data.frame.matrix(table(clust_label, ct_label))
  rand_size_table <- apply(rand_ct_raw_table, 1, sum)
  rand_ct_pct_table <- rand_ct_raw_table / rand_size_table

  entropy_vector <- apply(rand_ct_pct_table, 1, function(ct_pct_table_line) {
    entropy::entropy.empirical(ct_pct_table_line, unit = "log")
  })
  harmonized_entropy_vector <- entropy_vector / log(rand_size_table)
  # handling cluster with size == 1 (log(1) == 0) which cause error!
  harmonized_entropy_vector[rand_size_table == 1] <- 0.0
  harmonized_entropy_vector
}