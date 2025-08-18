#' ----------------------------------------------------------------------------
#' @title:      Create OrthologousGroups for OMA pipeline
#' @description
#' Generate `OrthologousGroups_one_to_one.txt`
#'          `OrthologousGroups_many_to_many.txt`
#'
#' The groups of orthologs are given as one per row, starting with a unique
#' group identifier, followed by all group members, all separated by tabs.
#'
#' @details:
#' This script creates one-to-one and many-to-many orthogroup table for OMA
#' using the results in `OMA/PairwiseOrthologs`. The result is stored in
#' `OMA/Orthogroups`
#'
#' @author  David Blickenstorfer, Technau Group (2025)
#' @created 16/07/2025
#'
#' @copyright Technau Group, University of Vienna, 2025
#' @copyright Boeva Lab, ETH Zurich, 2025
# -----------------------------------------------------------------------------
library(dplyr)
library(tidyr)
library(Matrix)
library(hash)
library(igraph)
# -----------------------------------------------------------------------------
# Step 1: Create adcency matrix using results in PairwiseOrthologs/
# -----------------------------------------------------------------------------
# Load list of all genes and generate hashing tables
gene_list_path <- file.path("results", "OMA", "Map-SeqNum-ID.txt")
gene_list <- read.table(file = gene_list_path, header = FALSE, sep = "\t")
species <- gene_list[[1]]
genes <- gene_list[[3]]
indices <- seq_along(genes)
gene_to_id <- hash(keys = genes, values = indices)
id_to_gene <- hash(keys = indices, values = genes)
gene_to_species <- hash(keys = genes, values = species)
unique_species <- unique(species)

# Load PairwiseOrthologs and return edges
oma_input_path <- file.path("results", "OMA", "PairwiseOrthologs")
oma_input_list <- list.files(path = oma_input_path, full.names = TRUE)
row_list <- sapply(oma_input_list, function(input) {
  oma_input <- read.table(file = input, header = FALSE, sep = "\t")
  row_gene <- oma_input[[3]]
})
col_list <- sapply(oma_input_list, function(input) {
  oma_input <- read.table(file = input, header = FALSE, sep = "\t")
  col_gene <- oma_input[[4]]
})

# Use the edges row -> gene to generate adjency matrix
row_genes <- unlist(row_list, use.names = FALSE)
row_id <- sapply(row_genes, function(gene) {
  gene_to_id[[gene]]
})
row_id <- as.numeric(row_id)
col_genes <- unlist(col_list, use.names = FALSE)
col_id <- sapply(col_genes, function(gene) {
  gene_to_id[[gene]]
})
col_id <- as.numeric(col_id)

# Create dgTMatrix and get Clusters
# The dgTMatrix considered directed graph (Hv -> Ac but not Ac -> Hv)
# We need to use mode = "weak" for clustering, else no cluster is found
n <- length(genes)
adjency <- sparseMatrix(i = row_id, j = col_id, x = 1, dims = c(n, n),
                        dimnames = list(genes, genes), giveCsparse = FALSE)
graph <- graph_from_adjacency_matrix(adjency,
                                     mode = "directed",
                                     weighted = NULL)

# -----------------------------------------------------------------------------
# Step 2: Generate many-to-many
#   (1) Remove genes without any connections
#   (2) The cluster must contain genes from more than one species
#    => this was ensured with the PairwiseOrthologs
# -----------------------------------------------------------------------------
# remove vertices with degree = zero
graph <- delete_vertices(graph, degree(graph) == 0)
clusters <- components(graph, mode = "weak")

# Write the clusters as dataframe
create_gene_row <- function(cluster_id, transcript_name) {
  # write data frame with cluster ID, transcriptID, geneID, species
  gene_name <- sub("\\.t[0-9]+$", "", transcript_name)
  species <- gene_to_species[[transcript_name]]
  data.frame(clusterID = cluster_id,
             transcriptID = transcript_name,
             geneID = gene_name,
             species = species)
}
#' --------------------------------------------
#' clusterID | transcriptID | geneID | species
#' --------------------------------------------
df_genes <- mapply(create_gene_row,
                   clusters$membership,
                   names(clusters$membership),
                   SIMPLIFY = FALSE,
                   USE.NAMES = FALSE)
df_genes <- do.call(rbind, df_genes)
#'--------------------------------------
#' OG     | "Nv"    | "Ac"    | "Hv"   |
#'--------+---------+---------+--------+
#' OG0001 | n1.t1   | a1, a2  | NA     |
#' OG0002 | n1.t2   | a3      | h3     |
#' ...... | ...     | ...     | ...    |
#'--------------------------------------
num_clusters <- clusters$no
df_cluster_list <- lapply(1:num_clusters, function(cluster_id) {
  # [INSIDE THE CLUSTER ID: GET THE WHOLE GENES TO CLUSTER]
  df_cluster_genes <- df_genes[df_genes$clusterID == cluster_id, ]
  df_cluster <- split(df_cluster_genes$transcriptID, df_cluster_genes$species)
  list(clusterID = cluster_id, df_cluster)
})

# -----------------------------------------------------------------------------
#' @title collapse_species
#' @brief collapse the species vector into single string separated with ","
#' @param cluster Cluster in df_cluster_list, contains following structure
#'                | clusterID | speciesA    | speciesB    | speciesC    |
#'                |     1     | c(g1, g2)   | c(g1)       |             |
#'                |     2     | c(g3, g4)   | c(g2, g3)   | c(g1, g2)   |
#' @return collapsed dataframe, column does not have vectors
#'                | OG        | speciesA    | speciesB    | speciesC    |
#'                | OG0001    | "g1", "g2"  | "g1"        | NA          |
#'                | OG0002    | "g3", "g4"  | "g2", "g3"  | "g1", "g2"  |
# -----------------------------------------------------------------------------
width <- nchar(as.character(num_clusters))
format <- sprintf("OG0%%0%dd", width)
collapse_species <- function(cluster) {
  # get inner species data
  species_data <- cluster[[2]]
  # collapse the vector into one string for each species column
  collapsed <- lapply(species_data, function(genes) {
    paste(genes, collapse = ",")
  })
  cluster_id <- cluster[[1]]
  og_id <- sprintf(format, cluster_id)
  collapsed$OG <- og_id
  collapsed <- as.data.frame(collapsed, stringsAsFactors = FALSE)
  collapsed <- collapsed[, c("OG",
                             setdiff(names(collapsed), "OG"))]
}
df_cluster <- bind_rows(lapply(df_cluster_list, collapse_species))

# store many-to-many
result_path <- file.path("results", "OMA", "OrthologousGroups")
dir.create(result_path, recursive = TRUE, showWarnings = FALSE)
filename <- file.path(result_path, "Orthologs_many_to_many.csv")
write.table(df_cluster, filename, quote = FALSE,
            row.names = FALSE, col.names = FALSE, sep = "\t")

# -----------------------------------------------------------------------------
# Step 3: Generate One-to-one
#   (1) Collapse all transcript into one geneID
#   (2) Ensure that transcripts project to the same orthogroup
#   (3) Assure number of species = size of cluster
# -----------------------------------------------------------------------------
# Remove all transcripts when mapping to different cluster
df_orthologs_list <- lapply(1:num_clusters, function(cluster_id) {
  # [INSIDE THE CLUSTER ID: GET THE WHOLE GENES TO CLUSTER]
  df_cluster_genes <- df_genes[df_genes$clusterID == cluster_id, ]
  df_cluster <- split(df_cluster_genes$geneID, df_cluster_genes$species)
  # Remove all dublicated genes in species list
  df_orthologs <- lapply(df_cluster, function(sp_genes) {
    unique(sp_genes)
  })
  list(clusterID = cluster_id, df_orthologs)
})
df_orthologs <- bind_rows(lapply(df_orthologs_list, collapse_species))
# remove rows when they contain "," and thus, multiple gens
df_one_to_ones <- df_orthologs[!apply(df_orthologs, 1, function(row) {
  any(grepl(",", row, fixed = TRUE), na.rm = TRUE)
}), ]
