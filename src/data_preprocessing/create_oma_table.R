#' ----------------------------------------------------------------------------
#' @title create_oma_table.R
#'
#' @description
#' Generate `OrthologsTable_OMA.tsv`
#'
#' The OMA groups are given as one per row, starting with a unique group
#' identifier, followed by all group members, all separated by tabs.
#'
#' @details:
#' This script creates an OMA group table using `OrthologousMatrix.txt` and
#' `OrthologousGroups.txt` located in the given directory ("results/OMA").
#' The result is stored in the same subdirectory.
#'
#' @author  David Blickenstorfer, Technau Group (2025)
#' @created 21/07/2025
#'
#' @copyright Technau Group, University of Vienna, 2025
#' @copyright Boeva Lab, ETH Zurich, 2025
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# Step 1: Read the OrthologousGroups.txt file
# -----------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1) {
  message(
    "❌ Invalid usage of create_oma_table.R [directory]\n",
    "[directory] : input path to `results/OMA/`\n"
  )
  stop("Usage: create_oma_table.R [directory]")
}
directory <- args[1]
# -----------------------------------------------------------------------------
#' Step 2: Read OMA input and create header
# -----------------------------------------------------------------------------
oma_input_path <- file.path(directory, "OrthologousGroups.txt")
oma_input <- read.delim(oma_input_path,
                        header = FALSE,
                        sep = "\t",
                        comment.char = "#",
                        stringsAsFactors = TRUE)
species_path <- file.path(directory, "OrthologousMatrix.txt")
species_names <- read.delim(species_path,
                            header = FALSE,
                            sep = "\t",
                            comment.char = "#",
                            stringsAsFactors = TRUE)
species_names <- as.vector(unlist(species_names[1, ]))
# get header (colnames) and the oma group names without header (cbind with data)
colnames <- c("OMA", species_names)
rownames <- oma_input[1]

# -----------------------------------------------------------------------------
#' Step 3: Create clean OMA group output in results/OMA
# -----------------------------------------------------------------------------
# take only the data
oma_table <- oma_input[, -c(1)]
# iterate row-wise (OMA group), get all genes and assign to species
oma_line_list <- apply(oma_table, 1, function(oma_group) {
  # Get OMA line string and then iterate over species names
  oma_sp_gene_list <- lapply(species_names, function(sp_name) {
    # take gene id (after "sp_name:" and before "\t")
    identifier <- paste0(sp_name, ":")
    sp_gene <- grep(paste0("^", identifier), oma_group, value = TRUE)
    # remove everything after \t, then everything before :
    sp_gene <- sub("\t.*", "", sp_gene)
    sp_gene <- sub(paste0("^", identifier), "", sp_gene)
    if (length(sp_gene) == 0) {
      sp_gene <- " "
    }
    sp_gene
  })
  oma_line <- paste(unlist(oma_sp_gene_list), collapse = "\t")
  strsplit(oma_line, "\t")
})
oma_vec_list <- lapply(oma_line_list, "[[", 1)
oma_table <- data.frame(do.call(rbind, oma_vec_list), stringsAsFactors = FALSE)
# add the col and rownames
oma_table <- cbind(rownames, oma_table)
colnames(oma_table) <- colnames

filename <- file.path(directory, "OrthologsTable_OMA.tsv")
write.table(oma_table, filename, quote = FALSE,
            row.names = FALSE, col.names = colnames, sep = "\t")

# Clean workspace
rm(list = ls())