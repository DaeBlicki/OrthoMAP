#' ----------------------------------------------------------------------------
#' @title create_orthofinder_expression.R
#'
#' @description
#' Constructs an orthogroups (HOG or OG) based expression matrix from multiple
#' a single-species Seurat object using Orthofinder orthology data. For
#' each OG group, the function identifies the most highly expressed gene per
#' cell and generates a matrix with OG groups as features and cells as columns.
#' The output is a modified Seurat object containing the new matrix.
#'
#' @details
#' This function integrates intra-species single-cell expression data with
#' orthology mappings from OrthoFinder to build a cross-compatible, gene
#' expression matrix. Ideal for cross-species comparisons or collapsing
#' redundant gene families into representative orthologous features.
#'
#' @param orthologous_groups Data frame or file path to `N0.tsv` (OrthoFinder)
#'                           stored in `Phylogenetic_Hierarchical_Orthogroups`
#'                           mapping OG groups to species-specific gene IDs.
#' @param datatable read csv file from `data/scsRNA_metadata.csv`
#' @param single_cell_directory Input directory where the seurat data are stored
#' @param annotation_table_directory Path to the directory where the annotation
#'                                   tables are stored
#' @param output_directory Output directory of OG expression matrix.
#' @param verbose Logical; if TRUE, prints progress messages.
#'
#' @return A Seurat object with a new expression matrix (OG groups × cells)
#'         stored in the output directory
#'
#' @importFrom Seurat CreateSeuratObject
#' @importFrom Seurat GetAssayData
#'
#' @seealso \code{\link[Seurat]{CreateSeuratObject}}
#' @references \url{https://cran.r-project.org/web/packages/Seurat/index.html}
#' @references \url{https://github.com/davidemms/OrthoFinder}
#'
#' @examples
#' \dontrun{
#' og_matrix <- create_orthofinder_expression(
#'     orthologous_groups = "N0.tsv",
#'     datatable = read.csv("data/scsRNA_metadata.csv"),
#'     single_cell_directory = "data/sc_objects",
#'     annotation_table_directory = "data/annotation_table"
#'     output_directory = "results/Orthogroup_expression_matrices/OrthoFinder",
#'     ouput_name = "Alldata.Robj",
#'     verbose = TRUE
#' )
#' }
#'
#' @author  David Blickenstorfer, Technau Group (2025)
#' @created 24/06/2025
#'
#' @copyright Technau Group, University of Vienna, 2025
#' @copyright Boeva Lab, ETH Zurich, 2025
#'
#' @export
create_orthofinder_expression <- function(
  orthologous_groups,
  datatable,
  single_cell_directory,
  annotation_table_directory,
  output_directory = getwd(),
  output_name = "Alldata.Robj",
  verbose = TRUE
) {
  # Import Seurat library in case
  if (!"Seurat" %in% loadedNamespaces()) {
    message("Load Seurat library in Environment\n")
    library(Seurat)
  }

  ## create OG tool variables
  # Only read the OG lines after comments (readLines and grep)
  #       input[1] = OMA000001\tAc:scaffold4\tHv:HVAEP1
  og_input <- read.delim(orthologous_groups,
                         header = TRUE,
                         sep = "\t",
                         stringsAsFactors = FALSE)

  og_names <- og_input[[2]]
  og_counts_rownames <- unique(og_names)
  og_to_index <- setNames(seq_along(og_counts_rownames),
                          og_counts_rownames)

  if (verbose) {
    cat("Found OG Orthogroups : ", length(og_counts_rownames), "\n\n")
  }

  # -------------------------------------------------------------------------
  #' @title orthofinder_expression_helper
  #' @brief generate OG expression matrix per species
  #' @param file_id String, species file name
  #' @param sc_version String, single-cell data file version
  #' @param orig_ident String, source of the data (laboratory)
  #' @param cell_type String, cell type of cell
  #' @return List of three:
  #'        - sp_oma_matrix: dgTMatrix, OMA vs. cells expression for species
  #'        - sp_metadata  : Meta.data, meta data of the species
  #'        - sp_cells     : name of the cells
  # -------------------------------------------------------------------------
  orthofinder_expression_helper <- function(file_id,
                                            sc_version,
                                            orig_ident,
                                            cell_type) {
    # Inform the user: Which species is getting processed
    if (verbose) {
      cat("Create orthogroup expression matrix for species: ", file_id, "\n")
      cat("[Pre-Processing] Begin at ", date(), "\n")
    }
    # load the seurat object
    seurat_obj <- get_single_cell_data(file_id, sc_version,
                                       single_cell_directory)

    # generate gene_expression matrix and get tool variables
    sp_counts <- Seurat::GetAssayData(seurat_obj, layer = "counts")
    genes_seurat <- rownames(sp_counts)
    sp_genes <- og_input[[file_id]]
    annotation_table <- get_annotation_table(file_id,
                                             annotation_table_directory)
    genes_fasta <- unlist(annotation_table[genes_seurat])

    # create lookup (hash-table) from genes to orthologous group
    lookup_pairs <- lapply(seq_along(sp_genes), function(i) {
      parts <- strsplit(sp_genes[i], ",")[[1]]
      # For each gene in parts, assign og_names[i]
      sapply(parts, function(gene) {
        og_names[i]
      })
    })
    genes_to_og_table <- unlist(lookup_pairs)
    genes_to_og_table <- collapse_transcript(file_id,
                                             genes_to_og_table,
                                             verbose)

    ## create Orthogroup vs. cells matrix
    # 1. For each cell i, get the expression value for gene j
    # 2. For each gene j, get the Orthogroup k (lookup_table)
    # 3. Get expression value for cell i in orthogroup k in og_matrix
    # 4. Compare expression values
    # convert into dgRMatrix to focus on genes
    sp_counts <- as(sp_counts, "RsparseMatrix")

    if (verbose) {
      cat("[Matrix generation] Begin at ", date(), "\n")
    }
    sp_og_matrix <- convert_gene_to_og_expression(sp_counts,
                                                  genes_fasta,
                                                  og_counts_rownames,
                                                  genes_to_og_table,
                                                  og_to_index,
                                                  verbose)
    sp_og_matrix <- as(sp_og_matrix, "CsparseMatrix")
    sp_metadata <- generate_species_metadata(seurat_obj,
                                             file_id,
                                             sc_version,
                                             orig_ident,
                                             cell_type)
    sp_cell_names <- colnames(sp_counts)

    list(sp_og_matrix = sp_og_matrix,
         sp_metadata = sp_metadata,
         sp_cell_names = sp_cell_names)
  }
  # end function

  # Apply matrix multiplication
  file_ids <- datatable$file_id
  sc_versions <- datatable$sc_version
  orig_idents <- datatable$orig.ident
  cell_types <- datatable$ID.separate
  orthofinder_results <- mapply(orthofinder_expression_helper,
                                file_ids,
                                sc_versions,
                                orig_idents,
                                cell_types,
                                SIMPLIFY = FALSE)

  # concartenate everything
  if (verbose) {
    cat("[Create Seurat Object] Begin at ", date(), "\n")
  }
  # create seurat obj and store in output_directory
  alldata_orthofinder <- merge_into_seurat_object(orthofinder_results,
                                                  og_counts_rownames,
                                                  verbose = verbose)
  if (verbose) {
    cat("Create and store seurat object in", output_directory, "\n\n")
  }

  filename <- file.path(output_directory, "Alldata.Robj")
  save(alldata_orthofinder, file = filename)
}