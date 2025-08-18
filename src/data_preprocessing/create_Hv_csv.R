#' ----------------------------------------------------------------------------
#' @title:      Create Annotation Table for Hv
#' @description Generate `Hv.csv`
#'
#' @details:
#' The seurat file `Nv2_Cole2024.Robj` used preprocessed the raw data
#' for the paper "Updated single cell reference atlas for the starlet anemone
#' Nematostella vectensis (2024)" [https://doi.org/10.1186/s12983-024-00529-z].
#' This Code will reverse the gene names used in the paper and transform into
#' the names used for the Orthogroup Analysis.
#'
#' @author  David Blickenstorfer, Technau Group (2025)
#' @created 07/07/2025
#'
#' @copyright Technau Group, University of Vienna, 2025
#' @copyright Boeva Lab, ETH Zurich, 2025
# -----------------------------------------------------------------------------
library(tidyverse)  # tidyverse (v.2.0.0), use dplyr, purrr

# load lookup table and `Nv2_Cole2024.Robj`
data <- read.csv(file.path("download", "HVAEP.csv"))
base_id <- data$Base.Gene.ID
base_id <- unique(base_id)
gene_name <- lapply(base_id, function(x) {
  name <- paste0("HVAEP1-", x)
  name
})
gene_id <- sapply(seq_along(base_id), function(i) {
  scaffold <- case_when(
    i >= 1     & i <= 2665  ~ 1,
    i >= 2666  & i <= 5082  ~ 2,
    i >= 5083  & i <= 7345  ~ 3,
    i >= 7346  & i <= 9051  ~ 4,
    i >= 9052  & i <= 10389 ~ 5,
    i >= 10390 & i <= 12371 ~ 6,
    i >= 12372 & i <= 13703 ~ 7,
    i >= 13704 & i <= 15525 ~ 8,
    i >= 15526 & i <= 17676 ~ 9,
    i >= 17677 & i <= 19351 ~ 10,
    i >= 19352 & i <= 21260 ~ 11,
    i >= 21261 & i <= 24043 ~ 12,
    i >= 24044 & i <= 25825 ~ 13,
    i >= 25826 & i <= 27142 ~ 14,
    i >= 27143 & i <= 28917 ~ 15,
    TRUE                   ~ NA_integer_
  )
  paste0("HVAEP", scaffold, ".", base_id[i])
})

df <- data.frame(baseID = as.character(base_id),
                 geneID = as.character(gene_id),
                 gene.name = as.character(gene_name),
                 stringsAsFactors = FALSE)

filename <- file.path("data", "annotation_table", "Hv.csv")
write.csv(df, filename, row.names = FALSE)

lines <- readLines(filename)
clean_lines <- gsub('"', '', lines)
writeLines(clean_lines, filename)

# Clean workspace
rm(list = ls())