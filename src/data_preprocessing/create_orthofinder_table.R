#' ----------------------------------------------------------------------------
#' @title create_orthofinder_table.R
#'
#' @description
#' Generate `OrthologsTable_HOGs.tsv` and `OrthologsTable_one_to_one.tsv`
#'
#' The groups of orthologs are given as one per row, starting with a unique
#' group identifier, followed by all group members, all separated by tabs.
#'
#' @details:
#' This script creates 1:1 and all-to-all orthologs table for OrthoFinder
#' using the results in `Phylogenetic_Hierarchical_Orthogroups/N0.tsv`. The
#' result is stored in the same subdirectory.
#'
#' @author  David Blickenstorfer, Technau Group (2025)
#' @created 17/07/2025
#'
#' @copyright Technau Group, University of Vienna, 2025
#' @copyright Boeva Lab, ETH Zurich, 2025
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# Step 1: Read the N0.tsv file and create clean data.frame
# -----------------------------------------------------------------------------
library(tidyverse)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1) {
  message(
    "❌ Invalid usage of create_orthofinder_table.R [directory]\n",
    "[directory] : input path to `Phylogenetic_Hierarchical_Orthogroups/`\n"
  )
  stop("Usage: create_orthofinder_table.R [directory]")
}
directory <- args[1]
# -----------------------------------------------------------------------------
#' Step 2: Create clean list of HOGs (all:all)
#'
#' OrthologsTable_HOGs.tsv
#' | HOG       | speciesA          | speciesB    | speciesC    |
#' *-----------*-------------------*-------------*-------------*
#' | HOG0001   | "g1.t1", "g2.t1"  | "g1"        | NA          |
#' | HOG0002   | "g3.t1", "g3.t2"  | "g2", "g3"  | "g1", "g2"  |
# -----------------------------------------------------------------------------
orthofinder_input_path <- file.path(directory, "N0.tsv")
orthofinder_input <- read.delim(orthofinder_input_path,
                                header = TRUE,
                                sep = "\t",
                                comment.char = "#",
                                stringsAsFactors = TRUE)
colnames <- strsplit(readLines(orthofinder_input_path, n = 1), "\t")[[1]]
# remove OG and Gene.Tree.Parent.Clade
colnames <- colnames[-c(2, 3)]
hog_df <- orthofinder_input[, -c(2, 3)]
non_empty_count <- rowSums(!is.na(hog_df[-1]) & hog_df[-1] != "")
hog_df <- hog_df[non_empty_count > 1, ]
# rename the HOGs using dynamic padding
num_hogs <- length(hog_df[[1]])
hog_sequence <- 0:(num_hogs - 1)
width <- nchar(as.character(num_hogs))
format <- sprintf("HOG%%0%dd", width)
hog_names <- unlist(lapply(hog_sequence, function(id) {
  sprintf(format, id)
}))
hog_df[, 1] <- hog_names
# store the result as `.tsv`
filename <- file.path(directory, "OrthologsTable_HOG.tsv")
write.table(hog_df, filename, quote = FALSE,
            row.names = FALSE, col.names = colnames, sep = "\t")

# -----------------------------------------------------------------------------
#' Step 2.2: Create clean list of OGs (all:all)
#'
#' OrthologsTable_HOGs.tsv
#' | OG       | speciesA          | speciesB    | speciesC    |
#' *-----------*-------------------*-------------*-------------*
#' | OG0001   | "g1.t1", "g2.t1"  | "g1"        | NA          |
#' | OG0002   | "g3.t1", "g3.t2"  | "g2", "g3"  | "g1", "g2"  |
# -----------------------------------------------------------------------------
og_df <- orthofinder_input[, -c(1, 3)]
colnames <- strsplit(readLines(orthofinder_input_path, n = 1), "\t")[[1]]
# remove HOG and Gene.Tree.Parent.Clade
colnames <- colnames[-c(1, 3)]
og_df <- og_df %>% dplyr::group_by(OG) %>% dplyr::summarise(
  across(everything(),
  ~ paste(na.omit(.),
  collapse = " ")),
  .groups = "drop"
)
non_empty_count <- rowSums(!is.na(og_df[-1]) & og_df[-1] != "")
og_df <- og_df[non_empty_count > 1, ]
# rename the HOGs using dynamic padding
num_ogs <- length(og_df[[1]])
hog_sequence <- 0:(num_ogs - 1)
width <- nchar(as.character(num_ogs))
format <- sprintf("OG%%0%dd", width)
og_names <- unlist(lapply(hog_sequence, function(id) {
  sprintf(format, id)
}))
og_df[, 1] <- og_names
filename <- file.path(directory, "OrthologsTable_OG.tsv")
write.table(og_df, filename, quote = FALSE,
            row.names = FALSE, col.names = colnames, sep = "\t")

# -----------------------------------------------------------------------------
#' Step 3: Create One-to-One Orthologs table
#'    (1) Collapse isoforms into genes
#'    (2) Identify genes mapping to multiple OG
#'    (3.1) Remove HOG when multiple genes of species is in entry of dataframe
#'    (3.2) Remove HOG when single gene mapped to multiple OG before filtering
#'    (4) Remove HOG when #gene != #species
#'    (5) Rename HOG into OOG (One-to-One orthologous groups)
#'
#' OrthologsTable_one_to_one.tsv
#' | OG        | speciesA    | speciesB    | speciesC    |
#' *-----------*-------------*-------------*-------------*
#' | OG0001    | "g2"        | "g1"        | ""          |
#' | OG0002    | "g3"        | "g3"        | "g1"        |
# -----------------------------------------------------------------------------
# (1) Collapse isoforms into genes
oog_df <- sapply(hog_df, function(transcripts) {
  # get geneID by remove .t%d
  lapply(transcripts, function(isoforms_in_species) {
    # get isoform, remove `.t%d`, collapse dublicates, and return as string
    isoforms <- unlist(strsplit(as.character(isoforms_in_species), ",\\s*"))
    genes <- gsub("\\.t[0-9]+", "", isoforms)
    unique_genes <- unique(genes)
    paste(unique_genes, collapse = ",")
  })
})
oog_df <- as.data.frame(oog_df)
genes_df <- oog_df[, -c(1)]

# (2) Mark genes mapping to multiple OG
# get genes in dataframe as list and remove ""
all_genes_entries <- unlist(genes_df)
all_genes_entries <- all_genes_entries[all_genes_entries != ""]
all_genes <- unlist(strsplit(all_genes_entries, ",\\s*"))
genes_count <- table(all_genes)
genes_dubs_idx <- which(genes_count > 1)
genes_dubs <- names(genes_count[genes_dubs_idx])

# (3) Mark the dataframe when entry is in genes_dubs or contains multiple genes
marked_df <- sapply(genes_df, function(genes_entry) {
  sapply(genes_entry, function(genes) {
    gene <- as.character(genes)
    gene_list <- unlist(strsplit(gene, ",\\s*"))
    # Check (3.1) and then (3.2)
    if (length(gene_list) > 1 || gene %in% genes_dubs) {
      gene <- FALSE
    }
    gene
  })
})

# (4) Remove HOGs when marked with FALSE
oog_df <- oog_df[!apply(marked_df, 1, function(hog) any(hog == FALSE)), ]
oog_df <- sapply(oog_df, function(gene) {
  sapply(gene, function(x) as.character(x))
})
oog_df <- as.data.frame(oog_df)

# (5) Rename HOG into OOG (one-to-one orthologous group)
num_oogs <- length(oog_df[[1]])
oog_sequence <- 0:(num_oogs - 1)
width <- nchar(as.character(num_oogs))
format <- sprintf("OOG%%0%dd", width)
oog_names <- unlist(lapply(oog_sequence, function(id) {
  sprintf(format, id)
}))
oog_df[, 1] <- oog_names
colnames[1] <- "OOG"

# store the result as `.tsv`
filename <- file.path(directory, "OrthologsTable_OOG.tsv")
write.table(oog_df, filename, quote = FALSE,
            row.names = FALSE, col.names = colnames, sep = "\t")

# Clean workspace
rm(list = ls())