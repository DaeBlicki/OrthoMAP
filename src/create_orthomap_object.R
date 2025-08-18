#' ----------------------------------------------------------------------------
#' @title create_orthomap_object.R
#'
#' @description
#' Constructs an orthogroups (HOG, OOG or OMA) based expression matrix from
#' multiple single-species Seurat object using processed results from the
#' orthology interference analysis. The function identifies the most highly
#' expressed gene per cell and generates a matrix with OG groups as features
#' and cells as columns. The output is a modified Seurat object containing
#' the merged orthologous cluster matrix and meta.data.
#'
#' @details
#' This function integrates intra-species single-cell expression data with
#' orthology mappings from orthology interference tools to build a Seurat
#' Object compatible for cross-species comparisons or collapsing redundant
#' gene families into representative orthologous features.
#'
#' @param orthologous_groups_path Path to processed orthology interference
#'                                analysis. Maps OG groups to gene IDs. Maps
#'                                Orthogroup to species-specific gene IDs.
#' @param datatable read csv file from `data/scsRNA_metadata.csv`
#' @param single_cell_directory Input directory where the seurat data are stored
#' @param annotation_table_directory Path to the directory where the annotation
#'                                   tables are stored
#' @param output_directory Output directory of OG expression matrix.
#' @param output_name Character Vector, Name of the OrthoMAP Object
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
#'
#' @examples
#' \dontrun{
#' og_matrix <- create_orthomap_object(
#'     orthologous_groups_path = "OrthologsTable_HOG.txt",
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
#' @created 21/07/2025
#'
#' @copyright Technau Group, University of Vienna, 2025
#' @copyright Boeva Lab, ETH Zurich, 2025
#'
#' @export
#' ----------------------------------------------------------------------------
create_orthomap_object <- function(
  orthologous_groups_path,
  datatable,
  single_cell_directory,
  annotation_table_directory,
  output_directory = getwd(),
  output_name = "Alldata.Robj",
  verbose = TRUE
) {
  # ---------------------------------------------------------------------------
  # 1. Read input variables and generate hashmaps
  # ---------------------------------------------------------------------------
  ## create OG tool variables
  og_input <- read.delim(orthologous_groups_path,
                         header = TRUE,
                         sep = "\t",
                         stringsAsFactors = FALSE)

  og_counts_rownames <- og_input[[1]]
  og_to_index <- setNames(seq_along(og_counts_rownames),
                          og_counts_rownames)

  if (verbose) {
    cat("Found OG Orthogroups : ", length(og_counts_rownames), "\n\n")
  }

  # -------------------------------------------------------------------------
  #' @title orthomap_expression_helper
  #' @brief generate OG expression matrix per species
  #' @param file_id String, species file name
  #' @param sc_version String, single-cell data file version
  #' @param orig_ident String, source of the data (laboratory)
  #' @param cell_type String, cell type of cell
  #' @param cell_type_separate String, more precise celltype
  #' @return List of three:
  #'        - sp_oma_matrix: dgTMatrix, OG vs. cells expression for species
  #'        - sp_metadata  : Meta.data, meta data of the species
  #'        - sp_cells     : name of the cells
  # -------------------------------------------------------------------------
  orthomap_expression_helper <- function(file_id,
                                         sc_version,
                                         orig_ident,
                                         cell_type,
                                         cell_type_separate) {
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

    # create lookup (Named vectors) from genes to orthologous group
    lookup_pairs <- lapply(seq_along(sp_genes), function(i) {
      parts <- strsplit(sp_genes[i], ",")[[1]]
      # For each gene in parts, assign og_names[i]
      sapply(parts, function(gene) {
        og_counts_rownames[i]
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

    if (verbose) cat("[Matrix generation] Begin at ", date(), "\n")
    sp_og_matrix <- convert_gene_to_og_expression(sp_counts,
                                                  genes_fasta,
                                                  og_counts_rownames,
                                                  genes_to_og_table,
                                                  og_to_index,
                                                  verbose)
    sp_og_matrix <- as(sp_og_matrix, "CsparseMatrix")
    sp_metadata <- generate_species_metadata(
      sc_object = seurat_obj,
      file_id = file_id,
      orig_ident = orig_ident,
      cell_type = cell_type,
      cell_type_separate = cell_type_separate
    )
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
  cell_types <- datatable$IDs
  cell_type_separates <- datatable$ID.separate
  orthology_results <- mapply(orthomap_expression_helper,
                              file_id = file_ids,
                              sc_version = sc_versions,
                              orig_ident = orig_idents,
                              cell_type = cell_types,
                              cell_type_separate = cell_type_separates,
                              SIMPLIFY = FALSE)

  # concartenate everything
  if (verbose) cat("[Create Seurat Object] Begin at ", date(), "\n")
  orthomap_obj <- merge_into_seurat_object(orthology_results,
                                           og_counts_rownames,
                                           verbose = verbose)
  # make IDs, ID.separate, and species as factors with levels
  orthomap_obj$IDs <- factor(orthomap_obj$IDs,
                             levels = unique(orthomap_obj$IDs))
  orthomap_obj$ID.separate <- factor(orthomap_obj$ID.separate,
                                     levels = unique(orthomap_obj$ID.separate))
  orthomap_obj$species <- factor(orthomap_obj$species,
                                 levels = unique(orthomap_obj$species))

  # save orthomap object
  filename <- file.path(output_directory, output_name)
  if (verbose) cat("Store Seurat Object: ", filename, "\n\n")
  save(orthomap_obj, file = filename)
  orthomap_obj
}