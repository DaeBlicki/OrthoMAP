#'-----------------------------------------------------------------------------
#' @title   Estimate best leiden resolution for SAMap
#'
#' @concept
#' Estimates the best SAMap resolution based on the four bio conservations:
#' ARI for IDs, ARI for ID.separate, NMI for IDs, and NMI for ID.separate
#'
#' This script is not optimized and hardcoded
#'
#' @author  David Blickenstorfer, Technau Group (2025)
#' @created 11/08/2025
#'
#' Version: v1.0.0
#' Changelog:
#'  - v1.0.0: Initial release
#'
#' @copyright Technau Group, University of Vienna, 2025
#' @copyright Boeva Lab, ETH Zurich, 2025
#' ----------------------------------------------------------------------------
library(tidyverse)
library(aricode)
# -----------------------------------------------------------------------------
# Step 1: Load Data and Preparation
# -----------------------------------------------------------------------------
md <- read.csv(file.path("results", "raw_md.csv"), row.names = 1)
md$Cell_ID <- rownames(md)
resolutions <- c("0.5", "1.0", "1.5", "2.0", "2.5",
                 "3.0", "3.5", "4.0", "4.5", "5.0")
filenames <- file.path("results", "leiden_clusters",
                       paste0("resolution_", resolutions, ".csv"))
data_list <- lapply(filenames, function(x) {
  data <- read.csv(x, row.names = 1)
  data$Cell_ID <- rownames(data)
  data
})
data_list <- lapply(data_list, function(x) {
  merge(x, md, by = "Cell_ID", all.x = TRUE)
})

# Tool variables
md_colnames <- c("Cell_ID", "IDs", "ID.separate", "species", "orig.ident")

# -----------------------------------------------------------------------------
# Step 2: Estimation of ARI and NMI
# -----------------------------------------------------------------------------
# ARI for IDs
median_ari_bio <- unlist(lapply(data_list, function(x) {
  data <- apply(x[, !names(x) %in% md_colnames], 2, function(y) {
    aricode::ARI(y, x[, "IDs"])
  })
  median(data)
}))
# ARI for orig.ident
median_ari_batch <- unlist(lapply(data_list, function(x) {
  data <- apply(x[, !names(x) %in% md_colnames], 2, function(y) {
    aricode::ARI(y, x[, "orig.ident"])
  })
  median(data)
}))
# NMI for IDs
median_nmi_bio <- unlist(lapply(data_list, function(x) {
  data <- apply(x[, !names(x) %in% md_colnames], 2, function(y) {
    aricode::NMI(y, x[, "IDs"])
  })
  median(data)
}))
# NMI for orig.ident
median_nmi_batch <- unlist(lapply(data_list, function(x) {
  data <- apply(x[, !names(x) %in% md_colnames], 2, function(y) {
    aricode::ARI(y, x[, "orig.ident"])
  })
  median(data)
}))
# Create data frame
list_of_medians <- list(median_ari_bio, median_ari_batch,
                        median_nmi_bio, median_nmi_batch)
df <- do.call(rbind, list_of_medians)
colsums <- apply(df, 2, function(x) sum(unlist(x)))

cat("Best Resolution: ", resolutions[which.max(colsums)], "\n")

df <- rbind(df, colsums)
colnames(df) <- resolutions
rownames(df) <- c("cARI", "iARI", "cNMI", "iNMI", "Sum")
df <- as.data.frame(df)

filename <- file.path("results", "conservation_table.csv")
write.csv(df, filename)