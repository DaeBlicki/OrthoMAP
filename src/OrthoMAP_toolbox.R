#'-----------------------------------------------------------------------------
#' @title   OrthoMAP_toolbox.R
#' @concept Core utility functions for the OrthoMAP pipeline
#'
#' @description
#' This script defines a collection of utility functions used across the
#' OrthoMAP pipeline. While the pipeline’s behavior can vary depending on
#' the input data, these functions represent the invariant logic — stable,
#' reusable components that support consistent downstream analysis across
#' all input configurations.
#'
#' Intended as a sourceable script: load with `source("OrthoMAP_toolbox.R")`
#'
#' @author  David Blickenstorfer, Technau Group (2025)
#' @created 24/06/2025
#'
#' Version: v1.0.0
#' Last Updated: 16/07/2025
#' Changelog:
#'  - v1.0.0: Initial release
#'
#' @copyright Technau Group, University of Vienna, 2025
#' @copyright Boeva Lab, ETH Zurich, 2025
#' ----------------------------------------------------------------------------

# =============================================================================
# STEP 1: Create OrthoGroup against Celltype matrix
#     - get_single_cell_data
#     - get_annotation_table
#     - generate_species_metadata
#     - collapse_transcript
#     - merge_into_seurat_object
#     - visualize_qc_metric
# =============================================================================

#' ----------------------------------------------------------------------------
#' @title get_single_cell_data
#'
#' @description
#' Read single cell data with sc_version from species with the file_id.
#'
#' Created 01/07/2025 by David Blickenstorfer
#'
#' @param file_id String or char, name of current species in the list
#' @param sc_version String, single-cell data file version
#' @param input_directory path to the seurat objects
#'
#' @return return updated single-cell data
#'
#' @importFrom Seurat UpdateSeuratObject
#'
#' @export
#' ----------------------------------------------------------------------------
get_single_cell_data <- function(file_id, sc_version, input_directory) {
  # load the seurat object
  seurat_obj_path <- file.path(input_directory, paste0(file_id, ".Robj"))
  if (sc_version == "old") {
    seurat_obj <- readRDS(seurat_obj_path)
  } else {
    seurat_obj <- get(load(seurat_obj_path))
  }
  seurat_obj <- Seurat::UpdateSeuratObject(seurat_obj)
  seurat_obj
}

#' ----------------------------------------------------------------------------
#' @title get_annotation_table
#'
#' @description
#' Read the annotation table of file_id and return an hash table from gene.name
#' to gene.id. In other words, it maps from gene names in Seurat to the gene
#' names in FASTA.
#'
#' Created 03/07/2025 by David Blickenstorfer
#'
#' @param file_id String or char, name of current species in the list
#' @param input_directory path to the annotation table
#' @return return hashtable [seurat_names] -> [fasta_names]
#'
#' @export
#' ----------------------------------------------------------------------------
get_annotation_table <- function(file_id, input_directory) {
  # load the seurat object
  input_path <- file.path(input_directory, paste0(file_id, ".csv"))
  annotation_data <- read.csv(input_path, stringsAsFactors = FALSE)
  fasta_names <- annotation_data$geneID
  seurat_names <- annotation_data$gene.name
  gene_map <- setNames(fasta_names, seurat_names)
  gene_map
}

#' ----------------------------------------------------------------------------
#' @title generate_species_metadata
#'
#' @description
#' Generate species meta.data of the given single-cell data respecting the
#' common structure:
#'  - $species: species name stored as file_id
#'  - $orig.ident: data source or batch
#'  - $ID.separate: description of the cell type
#'
#' Created 02/07/2025 by David Blickenstorfer
#'
#' @param sc_object Updated Seurat Obj for the species
#' @param file_id String or char, name of current species in the list
#' @param sc_version String, single-cell data file version
#' @param orig_ident String, source of the data (laboratory)
#' @param cell_type String, cell type of cell
#' @param cell_type_separate String, more precise cell type of cell
#' @return concartenated metadata, list of seurat5.metadata
#'
#' @export
#' ----------------------------------------------------------------------------
generate_species_metadata <- function(sc_object,
                                      file_id,
                                      orig_ident,
                                      cell_type,
                                      cell_type_separate) {
  # create meta.data from scratch
  md <- data.frame(matrix(ncol = 4,
                          nrow = length(rownames(sc_object@meta.data))))
  rownames(md) <- rownames(sc_object@meta.data)
  colnames(md) <- c("species", "orig.ident", "ID.separate", "IDs")
  md$species <- file_id
  md$orig.ident <- as.character(sc_object@meta.data[[orig_ident]])
  md$ID.separate <- as.character(sc_object@meta.data[[cell_type_separate]])
  md$IDs <- paste0(file_id, ".", as.character(sc_object@meta.data[[cell_type]]))
  md
}

#' ----------------------------------------------------------------------------
#' @title collapse_transcript
#'
#' @description
#' Some FASTA files contains gene mapping to various proteins (isoforms). The
#' orthogroup analysis does not discriminate them and may map them to different
#' orthogroups. The geneID gets dismissed in case that multiple isoforms to
#' different orthogroups. The dismissed geneID is getting labelled as INVALID,
#' it returns the percentage of invalid genes in the geneID dataset.
#'
#' Created 04/07/2025 by David Blickenstorfer
#'
#' @param file_id String or char, name of current species in the list
#' @param genes_to_og_table Hashtable, maps from geneID.t% to orthogroup
#' @param verbose Logical, if TRUE then print information to the user
#'
#' @return Updated version Hashtable, maps from geneID to orthogroup
#'
#' @export
#' ----------------------------------------------------------------------------
collapse_transcript <- function(file_id,
                                genes_to_og_table,
                                verbose = TRUE) {
  # inform the user
  message("Invalidate genes in species ", file_id, "!\n")
  if (verbose) {
    cat("[Invalidate genes] Begin at ", date(), "\n")
  }
  # get isoform (geneID.t%) and gene (geneID)
  isoform_names <- names(genes_to_og_table)
  gene_names <- sub("\\.t[0-9]+$", "", isoform_names)
  genes_to_og_map <- tapply(
    genes_to_og_table,
    gene_names,
    function(orthogroup) {
      unique(orthogroup[!is.na(orthogroup)])
    }
  )
  genes_to_og_map <- unlist(genes_to_og_map)
  percent <- 1 - (length(genes_to_og_map) / length(gene_names))
  percent <- 100 * percent
  message("⚠️ Percentage of removed genes: ", percent, "%\n\n")
  if (verbose) {
    cat("⚠️ Percentage of removed genes: ", percent, "%\n")
  }
  genes_to_og_map
}

#' ----------------------------------------------------------------------------
#' @title merge_orthogroup_expression
#'
#' @description
#' Concartenate the orthogroup expression matrix of each species into one
#' Seurat Object. The resulting Seurat Object has following structure:bins
#'    - count: Orthogroup vs Cell expression matrix
#'      - rownames: names of the orthogroups
#'      - colnames: names of the cell
#'    - meta.data: Containing meta data of the single cell object
#'      - species: name of the species (file_id)
#'      - orig.ident: source or batch of the data
#'      - ID.separate: cell type of each cell
#'
#' Created 02/07/2025 by David Blickenstorfer
#'
#' @param species_result_list List of species containing:
#'        - sp_oma_matrix: dgTMatrix, OMA vs. cells expression for species
#'        - sp_metadata  : Meta.data, meta data of the species
#'        - sp_cells     : name of the cells
#' @param orthogroup_names List of strings with orthogroup names
#' @param verbose Logical, when TRUE then print information to the user
#'
#' @return Concartenated Seurat object
#'
#' @importFrom Seurat CreateSeuratObject
#'
#' @export
#' ----------------------------------------------------------------------------
merge_into_seurat_object <- function(species_result_list,
                                     orthogroup_names,
                                     verbose = TRUE) {
  # merge matrix
  if (verbose) cat("[Concatenate Matrix] Begin at ", date(), "\n")
  sp_counts_list <- purrr::map(species_result_list, "sp_og_matrix")
  counts <- do.call(cbind, sp_counts_list)

  # in case of single-species analysis
  if (is.null(counts)) {
    counts <- sp_counts_list
  }

  # merge metadata
  if (verbose) cat("[Concatenate Metadata] Begin at ", date(), "\n")
  sp_metadata_list <- purrr::map(species_result_list, "sp_metadata")
  metadata <- dplyr::bind_rows(sp_metadata_list)

  # merge cells
  if (verbose) cat("[Concatenate Cells] Begin at ", date(), "\n")
  sp_cell_names_list <- purrr::map(species_result_list, "sp_cell_names")
  cell_names <- unlist(sp_cell_names_list)

  # analyze and modify counts and metadata
  if (verbose) cat("[Analyze Count and Feature] Begin at ", date(), "\n")
  counts <- as(counts, "CsparseMatrix")
  metadata$nCount_OG <- Matrix::colSums(counts)
  metadata$nFeature_OG <- diff(counts@p)
  rownames(counts) <- orthogroup_names
  colnames(counts) <- cell_names

  # generate seurat objects
  if (verbose) cat("[CreateSeuratObject] Begin at ", date(), "\n")
  alldata <- CreateSeuratObject(counts = counts, meta.data = metadata)

  # filling
  alldata
}

# =============================================================================
# STEP 2: Apply Seurat Standard-Worflow on OrthoMAP Seurat Object
#   - generate_summary_table
#   - add_vlines
#   - preprocessing_orthomap_object
# =============================================================================

#' ----------------------------------------------------------------------------
#' @title generate_summary_table
#'
#' @description
#' Generate summary table of the given observables for all data. It calculates
#'  1) Median, Min, Max, 1st and 3rd Qu.
#'  2) Mean, Standard Deviation and standard Error of mean
#'  3) Confidence Intervals and Error for 95%
#'
#' @param data_list Numerical Vector of size n, Observable vector
#' @param data_names, Character Vector of size n, Name of the data
#' @param output_directory Path to the resulting .csv output
#' @param verbose Logical, when TRUE then print information to the user
#'
#' @return summary table, data.frame of the statistical results
#'
#' @export
#' ----------------------------------------------------------------------------
generate_summary_table <- function(data_list,
                                   data_names,
                                   output_directory,
                                   verbose) {
  if (verbose) cat("Generate summary table\n")
  summary <- lapply(data_list, function(x) {
    n <- length(x)
    # Generate Min, 1st. Qu, Median, 3rd Qu., and Max
    set_summary <- summary(x)
    # Generate Mean and Standard Deviation
    mean <- mean(x)
    sd <- sd(x)
    standard_error <- sd / sqrt(n)
    # Lowever and Upper 95%-Confidence
    conf_interval <- t.test(x)$conf.int
    # create Data.frame
    set_summary$N <- n
    set_summary$SD <- sd
    set_summary$SE <- standard_error
    set_summary$CI_lower <- conf_interval[1]
    set_summary$CI_upper <- conf_interval[2]
    as.data.frame(set_summary)
  })
  summary <- do.call(rbind, summary)
  summary <- summary[, c(7, 1, 2, 3, 4, 5, 6, 8, 9, 10, 11)]
  row.names(summary) <- data_names
  # store summary table in output directory
  filename <- file.path(output_directory, "summary_table.csv")
  write.csv(summary, filename, row.names = TRUE)
  if (verbose) cat("Stored summary as ", filename, "\n\n")
}

#' ----------------------------------------------------------------------------
#' @title add_vlines
#'
#' @description
#' Generate vlines depending on the choosen input of the summary table
#'
#' @param low Numeric, lower vline
#' @param high Numeric, higher vline
#' @param color Character Vector, name of color in ggplot2
#'
#' @return ggplot2::vlines
#'
#' @export
#' ----------------------------------------------------------------------------
add_vlines <- function(low, high, color) {
  vlines <- list()
  if (low == 0 && high > 0) {
    vlines <- ggplot2::geom_vline(xintercept = high,
                                  color = color,
                                  linetype = "dashed", size = 0.8)
  } else if (low > 0 && low < high) {
    vlines <- list(
      ggplot2::geom_vline(xintercept = low,
                          color = color,
                          linetype = "dashed", size = 0.8),
      ggplot2::geom_vline(xintercept = high,
                          color = color,
                          linetype = "dashed", size = 0.8)
    )
  } else if (low > 0 && high == 0) {
    vlines <- ggplot2::geom_vline(xintercept = low,
                                  color = color,
                                  linetype = "dashed", size = 0.8)
  }
  vlines
}

#' ----------------------------------------------------------------------------
#' @title preprocessing_orthomap_object
#'
#' @description
#' Remove the upper and lower threshhold in cells and orthogroups. This process
#' is iterative until convergence or maximal iteration is reached.
#'
#' @param orthomap_obj Orthogroup Seurat Object
#' @param config yaml configuration file
#' @param verbose Logical, when TRUE then print information to the user
#' @return Seurat Object
#'         - pruned expression matrix, dgCMatrix
#'         - pruned meta.data, data.frame
#'
#' @export
#' ----------------------------------------------------------------------------
preprocessing_orthomap_object <- function(orthomap_obj,
                                          config,
                                          verbose = TRUE) {
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
  if (verbose) print_preprocessing_config(config)

  # First Iteration (upper line has to cut only once)
  # Get information
  counts <- Seurat::GetAssayData(orthomap_obj, layer = "counts")
  before_num_cells <- ncol(counts)
  before_num_og <- nrow(counts)

  # filtering the number of orthogroups
  counts <- as(counts, "RsparseMatrix")
  # remove quantitative properties
  og_to_keep <- which(rowSums(counts > 0) > min_cell_in_og)
  counts <- counts[og_to_keep, ]
  orthomap_obj <- orthomap_obj[og_to_keep, ]
  # remove qualitative properties
  og_to_keep <- which(rowSums(counts) > min_count_in_og)
  counts <- counts[og_to_keep, ]
  orthomap_obj <- orthomap_obj[og_to_keep, ]
  # remove max values
  if (max_cell_in_og > 0) {
    og_to_keep <- which(rowSums(counts > 0) < max_cell_in_og)
    counts <- counts[og_to_keep, ]
    orthomap_obj <- orthomap_obj[og_to_keep, ]
  }
  if (max_count_in_og > 0) {
    og_to_keep <- which(rowSums(counts) < max_count_in_og)
    counts <- counts[og_to_keep, ]
    orthomap_obj <- orthomap_obj[og_to_keep, ]
  }
  counts <- as(counts, "CsparseMatrix")
  orthomap_obj$nFeature_OG <- Matrix::colSums(counts > 0)
  orthomap_obj$nCount_OG <- Matrix::colSums(counts)

  # filtering the number of cells
  orthomap_obj <- subset(x = orthomap_obj,
                         subset = nFeature_OG > min_og_in_cell &
                                  nCount_OG > min_count_in_cell) #nolint
  if (max_og_in_cell > 0) {
    orthomap_obj <- subset(x = orthomap_obj,
                           subset = nFeature_OG < max_og_in_cell)
  }
  if (max_count_in_cell > 0) {
    orthomap_obj <- subset(x = orthomap_obj,
                           subset = nCount_OG < max_count_in_cell)
  }

  # prepare iteration
  counts <- Seurat::GetAssayData(orthomap_obj, layer = "counts")
  diff_cells <- before_num_cells - ncol(counts)
  diff_orthogroups <- before_num_og - nrow(counts)
  it <- 2

  # print statistics
  if (verbose) {
    cat("Current iteration: 1/", max_iter, "\n")
    cat("Number relevant cells: ", ncol(counts), "/", before_num_cells, "\n")
    cat("Number relevant OG: ", nrow(counts), "/", before_num_og, "\n\n")
  }

  # iterate until convergence or hit max_iter
  while (it <= max_iter && (diff_cells) + (diff_orthogroups) > 0) {
    # Get information
    counts <- Seurat::GetAssayData(orthomap_obj, layer = "counts")
    before_num_cells <- ncol(counts)
    before_num_og <- nrow(counts)

    # filtering the number of orthogroups
    counts <- as(counts, "RsparseMatrix")
    og_to_keep <- which(rowSums(counts > 0) > min_cell_in_og)
    counts <- counts[og_to_keep, ]
    orthomap_obj <- orthomap_obj[og_to_keep, ]
    # remove qualitative properties
    og_to_keep <- which(rowSums(counts) > min_count_in_og)
    counts <- counts[og_to_keep, ]
    orthomap_obj <- orthomap_obj[og_to_keep, ]
    # update orthomap object
    counts <- as(counts, "CsparseMatrix")
    orthomap_obj$nFeature_OG <- Matrix::colSums(counts > 0)
    orthomap_obj$nCount_OG <- Matrix::colSums(counts)
    # filtering the number of cells
    orthomap_obj <- subset(x = orthomap_obj,
                           subset = nFeature_OG > min_og_in_cell &
                                    nCount_OG > min_count_in_cell) #nolint
    # statistics
    counts <- Seurat::GetAssayData(orthomap_obj, layer = "counts")
    diff_cells <- before_num_cells - ncol(counts)
    diff_orthogroups <- before_num_og - nrow(counts)

    # inform the user
    if (verbose) {
      cat("Current iteration: ", it, "/", max_iter, "\n")
      cat("Number relevant cells: ", ncol(counts), "/", before_num_cells, "\n")
      cat("Number relevant OG: ", nrow(counts), "/", before_num_og, "\n\n")
      if (diff_cells + diff_orthogroups == 0) {
        cat("Convergenced in ", it, " Iterations!\n\n")
      }
    }
    it <- it + 1
  }
  if (verbose && it > max_iter) {
    cat("⚠️ OrthoMAP did not converged in ", max_iter, "iterations!\n\n")
    message("⚠️ OrthoMAP did not converged in ", max_iter, "iterations!\n\n")
  }
  # return orthomap_obj
  orthomap_obj@meta.data$nFeature_RNA <- orthomap_obj@meta.data$nFeature_OG
  orthomap_obj@meta.data$nCount_RNA <- orthomap_obj@meta.data$nCount_OG
  orthomap_obj
}

#' ----------------------------------------------------------------------------
#' @title iterative_seurat_clustering
#'
#' @description
#' Iteratively optimize the resolution in Seurat::FindCluster until the number
#' of found clusters is equals the number of the meta.data of interest in the
#' given OrthoMAP Seurat Object. The presision is limited to 0.01 units.
#'
#' @details
#' The resolution is optimized iteratively by comparing the number of clusters
#' with the size of interest. When the number of clusters is more, the
#' resolution is decreased, when the number is less, the resolution increases.
#' To enhance the process, resolutions of 1, 0.1, and 0.01 are
#' included.
#'
#' @param orthomap_obj OrthoMAP Seurat Object
#' @param config yaml configuration file
#' @param verbose Logical, when TRUE then print information to the user
#'
#' @return Updated OrthoMAP Seurat Object
#'
#' @export
#' ----------------------------------------------------------------------------
iterative_clustering <- function(orthomap_obj, config, verbose) {
  # get metadata of interest
  interest <- orthomap_obj@meta.data[[config$seurat_standard$cluster.metadata]]
  interested_size <- length(unique(interest))
  # create precision vector and current resolution idx
  precision <- c(1, 0.1, 0.01)
  precision_idx <- 1
  iteration <- 1
  cur_resolution <- 1.0
  visited_resolutions <- c(0)
  # create while loop
  flag <- TRUE
  while (flag) {
    # perform the clustering
    if (verbose) {
      cat("Current iteration: ", iteration, "\n")
      cat("Current resolution: ", cur_resolution, "\n")
    }
    orthomap_obj <- Seurat::FindClusters(
      object = orthomap_obj,
      resolution = cur_resolution,
      random.seed = config$seurat_standard$cluster.seed,
      algorithm = config$seurat_standard$cluster.algorithm,
      verbose = verbose
    )
    # -------------------------------------------------------------------------
    # Case 1: Optimal size => Stop the iteration
    # -------------------------------------------------------------------------
    current_cluster_size <- length(unique(Seurat::Idents(orthomap_obj)))
    if (verbose) cat("Current cluster size: ", current_cluster_size, "\n\n")
    if (current_cluster_size == interested_size) flag <- FALSE
    # -------------------------------------------------------------------------
    # Case 2: Continue iteration
    #   Check if the next resolution was already calculated
    #   (A) Yes
    #     (A.1) Use res = 0.5*[res + (res ± precision)]
    #     (A.2) Increase precision and stop iteration when it's not possible
    #     (A.3) Apply clustering algorithm
    #   (B) No
    #     (B.1) check if number of cluster must increase or decrease
    #       increase => res + precision
    #       decrease => res - precision
    #     (B.2) Apply clustering algorithm
    # -------------------------------------------------------------------------
    visited_resolutions <- c(visited_resolutions, cur_resolution)
    # Apply B.1 and B.2
    if (current_cluster_size < interested_size) {
      next_resolution <- cur_resolution + precision[precision_idx]
    } else {
      next_resolution <- cur_resolution - precision[precision_idx]
    }
    # Check if you need Enter (A)
    if (next_resolution %in% visited_resolutions) {
      next_resolution <- (cur_resolution + next_resolution) / 2.0
      precision_idx <- precision_idx + 1
      # break if the precision can not be more precisely
      if (precision_idx > length(precision)) {
        # warn the user and break the repeat
        message("⚠️ Clustering did not converged using ",
                precision[precision_idx - 1], " precision! \n")
        if (verbose) {
          cat("⚠️ Clustering did not converged using ",
              precision[precision_idx - 1], " precision! \n")
        }
        flag <- FALSE
      }
    }
    # update resolution and number iteration
    if (flag) {
      cur_resolution <- next_resolution
      iteration <- iteration + 1
    }
  }
  cat("Estimated resolution for clustering: ", cur_resolution, "\n\n")
  orthomap_obj
}

#' ----------------------------------------------------------------------------
#' @title run_topdown_routine_helper
#'
#' @description
#' Helper function in the top-down and benchmark approach. It performs Seurat
#' standard clustering depending on the given configuration and step (coarse
#' or fine).
#'
#' @details
#' In the `coarse` clustering, the most variable features for each species are
#' concatenated. In the `fine` clustering, the most variable features for all
#' cells are used.
#'
#' @param orthomap_obj OrthoMAP Seurat Object
#' @param config_routine yaml configuration for the clustering
#' @param step Integer
#'    1: Coarse-Clustering on the whole orthomap object
#'    2: Fine-Clustering on the subset of the coarse clustering
#' @param verbose Logical, when TRUE then print information to the user
#'
#' @return Updated OrthoMAP Seurat Object
#'
#' @export
#' ----------------------------------------------------------------------------
run_topdown_routine_helper <- function(orthomap_obj,
                                       config_routine,
                                       step,
                                       verbose) {
  # Normalize data
  if (verbose) cat("[Log-Normalize] started at ", date(), "\n")
  orthomap_obj <- Seurat::NormalizeData(
    object = orthomap_obj,
    normalization.method = config_routine$normalization.method,
    scale.factor = config_routine$normalization.scale.factor,
    verbose = verbose
  )
  # Choose most variable features per species then concartenate
  if (verbose) cat("[Feature selection] started at ", date(), "\n")
  if (step == 1) {
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
  } else {
    orthomap_obj <- Seurat::FindVariableFeatures(
      object = orthomap_obj,
      selection.method = config_routine$selection.method,
      nfeatures = config_routine$selection.nfeatures,
      verbose = verbose
    )
  }
  if (verbose) {
    cat("Number of VariableFeatures : ",
        length(Seurat::VariableFeatures(orthomap_obj)), "\n")
  }
  # Scale the data
  if (verbose) cat("[Data scaling] started at ", date(), "\n")
  orthomap_obj <- Seurat::ScaleData(
    object = orthomap_obj,
    split.by = config_routine$scale.features,
    verbose = verbose
  )

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
    cat("Number of PCA dimensions : ", length(pca_dim), "\n")
  }

  if (step == 1 || step == 2) {
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
      cat("Number of harmony dimensions : ", length(harmony_dim), "\n")
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

    # Step 2.6 Find cluster
    if (verbose) cat("[Build Cluster] started at ", date(), "\n")
    orthomap_obj <- Seurat::FindClusters(
      object = orthomap_obj,
      resolution = config_routine$cluster.resolution,
      random.seed = config_routine$cluster.seed,
      algorithm = config_routine$cluster.algorithm,
      verbose = verbose
    )
  }
  orthomap_obj
}