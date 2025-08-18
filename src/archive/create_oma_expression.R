#' ----------------------------------------------------------------------------
#' @title create_oma_expression.R
#'
#' @description
#' Constructs an orthologous matrix alignment (OMA)-based expression matrix
#' from multiple a single-species Seurat object using OMA orthology data. For
#' each OMA group, the function identifies the most highly expressed gene per
#' cell and generates a matrix with OMA groups as features and cells as columns.
#' The output is a modified Seurat object containing the new matrix.
#'
#' @details
#' This function integrates intra-species single-cell expression data with
#' orthology mappings from OMA to build a cross-compatible, gene group–centric
#' expression matrix. Ideal for cross-species comparisons or collapsing
#' redundant gene families into representative orthologous features.
#'
#' @param orthologous_groups Data frame or file path to `OrthologousGroups.txt`
#'                           mapping OMA groups to species-specific gene IDs.
#' @param datatable read csv file from `data/scsRNA_metadata.csv`
#' @param single_cell_directory Path to the directory where the seurat data
#'                              are stored
#' @param annotation_table_directory Path to the directory where the annotation
#'                                   tables are stored
#' @param output_directory Output directory of OMA expression matrix.
#' @param verbose Logical; if TRUE, prints progress messages.
#'
#' @return A Seurat object with a new expression matrix (OMA groups × cells)
#'         stored in in output directory
#'
#' @importFrom Seurat CreateSeuratObject
#' @importFrom Seurat GetAssayData
#'
#' @seealso \code{\link[Seurat]{CreateSeuratObject}}
#' @references \url{https://cran.r-project.org/web/packages/Seurat/index.html}
#'
#' @examples
#' \dontrun{
#' oma_matrix <- create_OMA_expression(
#'     orthologous_groups = "OrthologousGroups.txt",
#'     datatable = read.csv("data/scsRNA_metadata.csv"),
#'     single_cell_directory = "data/sc_objects",
#'     annotation_table_directory = "data/annotation_table"
#'     output_directory = "results/Alison/Orthogroup_expression_matrices/OMA",
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
create_oma_expression <- function(
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

  ## create OMA tool variables
  # Only read the OMA lines after comments (readLines and grep)
  #       input[1] = OMA000001\tAc:scaffold4\tHv:HVAEP1
  oma_input <- readLines(orthologous_groups)
  oma_input <- oma_input[!grepl("^#", oma_input)]
  orthogroups_df <- read.delim(text = oma_input,
                               header = FALSE,
                               stringsAsFactors = FALSE)
  # Ortholog groups of all species
  oma_counts_rownames <- orthogroups_df[[1]]
  oma_to_index <- setNames(seq_along(oma_counts_rownames),
                           oma_counts_rownames)
  if (verbose) {
    cat("Found OMA Orthogroups : ", length(oma_counts_rownames), "\n\n")
  }

  # -------------------------------------------------------------------------
  #' @title oma_expression_matrix_helper
  #' @brief generate OMA expression matrix per species
  #' @param file_id String, species file name
  #' @param sc_version String, single-cell data file version
  #' @param orig_ident String, source of the data (laboratory)
  #' @param cell_type String, cell type of cell
  #' @return List of three:
  #'        - sp_oma_matrix: dgTMatrix, OMA vs. cells expression for species
  #'        - sp_metadata  : Meta.data, meta data of the species
  #'        - sp_cells     : name of the cells
  # -------------------------------------------------------------------------
  oma_expression_helper <- function(file_id,
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
    annotation_table <- get_annotation_table(file_id,
                                             annotation_table_directory)
    genes_fasta <- unlist(annotation_table[genes_seurat])

    ## create lookup (hash-table) from genes to Ortholog group
    # Use of lapply and "lambda" function
    lookup_pairs <- lapply(oma_input, function(line) {
      # split the line by tabs
      parts <- strsplit(line, "\t")[[1]]
      oma_id <- parts[1]
      # Extract only species-specific fields
      sp_genes <- grep(paste0("^", file_id, ":"), parts, value = TRUE)
      # Extract gene IDs after 'Nv:'
      gene_ids <- sub(paste0("^", file_id, ":"), "", sp_genes)
      if (length(gene_ids) > 0) {
        setNames(rep(oma_id, length(gene_ids)), gene_ids)
      } else {
        setNames(character(0), character(0))
      }
    })
    genes_to_og_table <- unlist(lookup_pairs)
    genes_to_og_table <- collapse_transcript(file_id,
                                             genes_to_og_table,
                                             verbose)

    ## create Orthogroup vs. cells matrix
    # 1. For each cell i, get the expression value for gene j
    # 2. For each gene j, get the Orthogroup k (lookup_table)
    # 3. Get expression value for cell i in orthogroup k in oma_matrix
    # 4. Compare expression values
    # convert into dgRMatrix to focus on genes
    sp_counts <- as(sp_counts, "RsparseMatrix")
    if (verbose) {
      cat("[Matrix generation] Begin at ", date(), "\n")
    }
    sp_oma_matrix <- convert_gene_to_og_expression(sp_counts,
                                                   genes_fasta,
                                                   oma_counts_rownames,
                                                   genes_to_og_table,
                                                   oma_to_index,
                                                   verbose)
    sp_oma_matrix <- as(sp_oma_matrix, "CsparseMatrix")
    sp_metadata <- generate_species_metadata(seurat_obj,
                                             file_id,
                                             sc_version,
                                             orig_ident,
                                             cell_type)
    sp_cell_names <- colnames(sp_counts)

    list(sp_og_matrix = sp_oma_matrix,
         sp_metadata = sp_metadata,
         sp_cell_names = sp_cell_names)
  }
  # end function

  # Apply matrix multiplication
  file_ids <- datatable$file_id
  sc_versions <- datatable$sc_version
  orig_idents <- datatable$orig.ident
  cell_types <- datatable$ID.separate

  oma_results <- mapply(oma_expression_helper,
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
  alldata_oma <- merge_into_seurat_object(oma_results,
                                          oma_counts_rownames,
                                          verbose = verbose)
  if (verbose) {
    cat("Create and store seurat object in", output_directory, "\n\n")
  }

  filename <- file.path(output_directory, output_name)
  save(alldata_oma, file = filename)
}