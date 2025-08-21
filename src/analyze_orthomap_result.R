#' ----------------------------------------------------------------------------
#' @title analyze_orthomap_result.R
#'
#' @description
#' Analyze the generated clusters in the OrthoMAP result object. It calculates
#' The distribution of species and celltype for each clusters. Furthermore, the
#' clusters are investigated using local Shannon Entropy, Adjusted Rand Index,
#' Normalized Mutual Information and Homogeneity Score. The cluster purity
#' are visualized using donut charts and pheatmap
#'
#' @details
#' This function integrates intra-species single-cell expression data with
#' orthology mappings from orthology interference tools to build a Seurat
#' Object compatible for cross-species comparisons or collapsing redundant
#' gene families into representative orthologous features.
#'
#' @param orthomap_obj The Seurat Object generated from step 1 of the pipeline.
#' @param output_directory path to statistical output
#' @param config RObject of .yaml, parameter for coloring
#' @param verbose Logical, if TRUE then print information to the user
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_histogram
#' @importFrom ggplot2 geom_vline
#' @importFrom scales
#'
#' @examples
#' \dontrun{
#'   analyze_orthomap_result(
#'     orthomap_obj = orthomap_obj,
#'     output_directory = orthomap_statistical_result
#'     config = config,
#'     verbose = verbose,
#'   )
#' }
#'
#' @author  David Blickenstorfer, Technau Group (2025)
#' @created 25/07/2025
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
analyze_orthomap_result <- function(orthomap_obj,
                                    output_directory,
                                    config,
                                    verbose) {
  # ---------------------------------------------------------------------------
  # Step 1: Create pheatmap for percentage of species and celltype per cluster
  # ---------------------------------------------------------------------------
  if (verbose) cat("[Pheatmap] started at", date(), "\n")
  result_format <- config$visualization$result_format
  ct_raw_table <- as.data.frame.matrix(table(orthomap_obj$orthomap_clusters,
                                             orthomap_obj$IDs))
  sp_raw_table <- as.data.frame.matrix(table(orthomap_obj$orthomap_clusters,
                                             orthomap_obj$species))
  size_summary <- apply(sp_raw_table, 1, sum)
  ct_summary <- apply(ct_raw_table, 2, sum)
  clust_ct_pct_table <- ct_raw_table / size_summary
  clust_sp_pct_table <- sp_raw_table / size_summary
  ct_pct_table <- t(ct_raw_table) / ct_summary
  dom_ct_table <- apply(ct_raw_table, 1, function(clust) {
    dom_idx <- which(clust == max(clust))
    colnames(ct_raw_table[dom_idx])
  })

  # check if the sum of the percentage is 1.9
  epsilon <- 1e-6
  invisible(apply(clust_sp_pct_table, 1, function(x) {
    if (sum(x) - 1.0 > epsilon) {
      message("⚠️ Species percentage failed in cluster!\n")
      cat("⚠️ Species percentage failed in cluster!\n")
    }
  }))
  invisible(apply(clust_sp_pct_table, 1, function(x) {
    if (sum(x) - 1.0 > epsilon) {
      message("⚠️ Celltype percentage failed in cluster!\n")
      cat("⚠️ Celltype percentage failed in cluster!\n")
    }
  }))
  # Tissue vs Cluster: scaled-tissue
  plotname <- paste0("Pheatmap_Tissue_Scaling", result_format)
  filename <- file.path(output_directory, plotname)
  gt <- pheatmap::pheatmap(ct_raw_table,
                           cluster_rows = TRUE, cluster_cols = TRUE,
                           scale = "column",
                           display_numbers = TRUE,
                           fontsize_number = 4,
                           fontsize_row = 4,
                           main = "Tissue vs. Cluster: column-scaled",
                           filename = filename)$gtable
  ggplot2::ggsave(filename, plot = gt)
  # Tissue vs Cluster: scaled-cluster
  plotname <- paste0("Pheatmap_Cluster_Scaling", result_format)
  filename <- file.path(output_directory, plotname)
  gt <- pheatmap::pheatmap(ct_raw_table,
                           cluster_rows = TRUE, cluster_cols = TRUE,
                           display_numbers = TRUE,
                           scale = "row",
                           fontsize_number = 4,
                           fontsize_row = 4,
                           main = "Tissue vs. Cluster: row-scaled",
                           filename = filename)$gtable
  ggplot2::ggsave(filename, plot = gt)
  # Tissue vs Cluster: percentage tissue in cluster
  plotname <- paste0("Pheatmap_pct_clus_in_tissue", result_format)
  filename <- file.path(output_directory, plotname)
  gt <- pheatmap::pheatmap(ct_pct_table,
                           cluster_rows = TRUE, cluster_cols = TRUE,
                           display_numbers = TRUE,
                           fontsize_number = 4,
                           main = "Cluster [%] in Tissue",
                           filename = filename)$gtable
  ggplot2::ggsave(filename, plot = gt)
  # Tissue vs Cluster: percentage tissue in cluster
  plotname <- paste0("Pheatmap_pct_tissue_in_clust", result_format)
  filename <- file.path(output_directory, plotname)
  gt <- pheatmap::pheatmap(clust_ct_pct_table,
                           cluster_rows = TRUE, cluster_cols = TRUE,
                           display_numbers = TRUE,
                           fontsize_number = 4,
                           main = "Tissue [%] in Cluster",
                           filename = filename)$gtable
  ggplot2::ggsave(filename, plot = gt)
  # plot species per cluster pheatmap
  plotname <- paste0("Pheatmap_species_in_clust", result_format)
  filename <- file.path(output_directory, plotname)
  gt <- pheatmap::pheatmap(clust_sp_pct_table,
                           cluster_rows = TRUE, cluster_cols = TRUE,
                           display_numbers = TRUE,
                           main = "Species Distribution in Cluster",
                           filename = filename)$gtable
  ggplot2::ggsave(filename, plot = gt)
  # ---------------------------------------------------------------------------
  # Step 2: Empirical evaluaion of the data
  #   1) Celltypes
  #       number of celltype contributes to statistical.coverage in cluster
  #   2) (Shannon) Entropy score
  #       measures diversity and uncertainty of distribution
  #   3) Adjusted Rand Index (ARI)
  #       single-cell standard statistical clustering index
  #   4) Silhouette score (for harmony, pca, and UMAPs)
  #       metric to evaluate quality of clustering in unsupervised ML algorithm
  #
  # Colnames of the empirical results table
  # Cluster, Celltypes, Entropy, ARI, Silhouette
  # ---------------------------------------------------------------------------
  stat_config <- config$statistic
  tissue_config_path <- file.path("data", "tissue_palette.csv")
  tissue_color_flag <- FALSE
  if (file.exists(tissue_config_path)) {
    palette <- read.csv(tissue_config_path, stringsAsFactors = TRUE)
    tissue_colors <- setNames(as.character(palette$Color), palette$Tissue)
    tissue_color_flag <- TRUE
  }
  cluster_charts_path <- file.path(output_directory, "Charts")
  dir.create(cluster_charts_path,
             recursive = TRUE,
             showWarnings = FALSE)
  # -------------------------------------------------------------------------
  #' @title create_cluster_plots
  #' @brief evaluate the tissue percentage of the clusters
  #' @param ct_raw_table_line line in ct_raw_table
  #' @return Character Vector, line in the empirical results table
  # -------------------------------------------------------------------------
  create_cluster_plots <- function(ct_raw_table_line) {
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
    coverage <- as.numeric(stat_config$coverage)
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
    # Make the plot
    palette_line <- ggplot2::scale_fill_brewer(palette = "Set3")
    if (tissue_color_flag) {
      palette_list <- lapply(tissue_vector, function(x) tissue_colors[x])
      palette <- unlist(palette_list)
      palette <- c(palette, "Others" = 	"#DCDCDC")
      palette_line <- ggplot2::scale_fill_manual(values = palette)
    }
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
      palette_line +
      ggplot2::theme(plot.background = element_rect(fill = "white",
                                                    color = NA)) +
      ggplot2::annotate(geom = "text", x = -1, y = 0, label = cluster_id,
                        size = 10, fontface = "bold", color = "black") +
      ggplot2::geom_rect() +
      ggplot2::geom_text(x = 3.25, aes(y=labelPosition, label=label), size=2) +
      ggplot2::coord_polar(theta="y", start = 0) +
      ggplot2::xlim(c(-1, 4))
    # save the dougnut plots
    plotname <- paste0("Dougnut_cluster_", cluster_id, result_format)
    plotname <- file.path(cluster_charts_path, plotname)
    ggplot2::ggsave(plotname, pie_plot, width = 6, height = 4)
    # return num_celltypes for coverage
    num_celltype
  }
  if (verbose) cat("[Doughnut Charts] started at ", date(), "\n")
  num_cell_vector <- apply(ct_raw_table, 1, create_cluster_plots)
  # -------------------------------------------------------------------------
  ## 2 Evaluate normalized Shannon Entropy
  # -------------------------------------------------------------------------
  if (verbose) cat("[Shannon Entropy] started at ", date(), "\n")
  entropy_vector <- apply(clust_ct_pct_table, 1, function(ct_pct_table_line) {
    entropy::entropy.empirical(ct_pct_table_line, unit = "log")
  })
  harmonized_entropy_vector <- entropy_vector / log(size_summary)
  # handling cluster with size == 1 (log(1) == 0) which cause error!
  harmonized_entropy_vector[size_summary == 1] <- 0.0

  stat_table <- cbind(dom_ct_table, num_cell_vector, harmonized_entropy_vector)
  colnames(stat_table) <- c("Dominant_Tissue", "NumTissues", "Shannon_Entropy")
  filename <- file.path(output_directory, "local_statistic.csv")
  write.csv(stat_table, filename, row.names = TRUE)
  # -------------------------------------------------------------------------
  ## 3 Differential gene expression
  # -------------------------------------------------------------------------
  if (FALSE) {
    if (verbose) cat("[Diff. Gene Expression] started at ", date(), "\n")
    features <- Seurat::VariableFeatures(orthomap_obj)
    if (as.character(config$statistic$dge.features) == "all") {
      features <- rownames(orthomap_obj)
    }
    all.markers <- Seurat::FindAllMarkers(
      orthomap_obj,
      group.by = "orthomap_clusters",
      test.use = config$statistic$dge.test.use,
      logfc.threshhold = config$statistic$dge.logfc.threshhold,
      random.seed = config$statistic$dge.random.seed,
      slot = config$statistic$dge.slot,
      features = features,
      return.thresh = config$statistic$dge.return.thresh,
      min.pct = config$statistic$dge.min.pct,
      only.pos = TRUE,
      verbose = verbose
    )
    clusters <- levels(orthomap_obj$orthomap_clusters)
    top_genes <- as.integer(config$statistic$dge.top.genes)
    if (verbose) cat("Select Top genes: ", top_genes, "\n")
    cluster_markers_list <- lapply(clusters, function(clust) {
      cluster_idx <- which(as.integer(all.markers$cluster) == clust)
      if (length(cluster_idx) > top_genes) {
        cluster_idx <- head(cluster_idx, n = top_genes)
      }
      all.markers[cluster_idx, "gene"]
    })
    dge_genes <- unique(unlist(cluster_markers_list))
    dplyr::glimpse(orthomap_obj$orthomap_clusters)
    titlename <- paste("Top ", top_genes, " DEGs")
    subtitlename <- paste("p-val < ", config$statistic$dge.return.thresh)
    dge_plot <- Seurat::DotPlot(
      orthomap_obj,
      group.by = "orthomap_clusters",
      features = dge_genes,
      scale.by = "radius",
      col.min = 0,
      col.max = 3
    ) + ggplot2::theme_bw() +
      Seurat::FontSize(6, 6) +
      ggplot2::theme(axis.text.x = element_text(angle = 90,
                                                vjust = 0.5,
                                                hjust = 1)) +
      ggplot2::labs(title = titlename, subtitle = subtitlename)

    plotname <- paste0("Differential_Gene_Expression", result_format)
    filename <- file.path(output_directory, plotname)
    ggplot2::ggsave(filename, plot = dge_plot, width = 12, height = 8)
  }
  orthomap_obj
}