#'-----------------------------------------------------------------------------
#' @title   Run Benchmark for OrthoMAP and SAMap of the best measurements
#'
#' @concept
#' Estimates the global clustering in form of Adjusted Random Index (ARI) and
#' Normalized Mutual Information (NMI), and local clustering in form of harm.
#' Shannon entropy of each cluster. The results are in form of dotplots.
#'
#' This script is not optimized and hardcoded
#'
#' @author  David Blickenstorfer, Technau Group (2025)
#' @created 15/08/2025
#'
#' Version: v1.0.0
#' Changelog:
#'  - v1.0.0: Initial release
#'
#' @copyright Technau Group, University of Vienna, 2025
#' @copyright Boeva Lab, ETH Zurich, 2025
#' ----------------------------------------------------------------------------
library(tidyverse)
# -----------------------------------------------------------------------------
# Step 1: Load the entropies of the best results
# -----------------------------------------------------------------------------
cat("[Step 1] Loading, started at ", date(), "\n")
entropy_path <- file.path("results", "OrthoMAP_Statistical_Results")
entropy_files <- c(
  file.path(entropy_path, "oma", "9", "local_statistic.csv"),
  file.path(entropy_path, "mmseqs2", "HOG", "9", "local_statistic.csv"),
  file.path(entropy_path, "diamond", "HOG", "9", "local_statistic.csv"),
  file.path(entropy_path, "mmseqs2", "OOG", "9", "local_statistic.csv"),
  file.path(entropy_path, "diamond", "OOG", "9", "local_statistic.csv"),
  file.path("results", "benchmark", "metrics", "samap_entropy.csv")
)
# -----------------------------------------------------------------------------
# Step 2: Extract the entropy
# -----------------------------------------------------------------------------
cat("[Step 2] Extracting, started at ", date(), "\n")
entropy_list <- lapply(entropy_files, function(file) {
  df <- read.csv(file)
  entropy <- df$Shannon_Entropy
  if (is.null(entropy)) entropy <- df$entropy
  entropy
})
# -----------------------------------------------------------------------------
# Step 3: Estimation of normalized Shannon Entropy
# -----------------------------------------------------------------------------
cat("[Step 3] Estimation of harmonized Shannon Entropy")
input_col <- c(rep("OMA", length(unlist(entropy_list[[1]]))),
               rep("HOG [mmseqs2]", length(unlist(entropy_list[[2]]))),
               rep("HOG [diamond]", length(unlist(entropy_list[[3]]))),
               rep("OOG [mmseqs2]", length(unlist(entropy_list[[4]]))),
               rep("OOG [diamond]", length(unlist(entropy_list[[5]]))),
               rep("SAMap", length(unlist(entropy_list[[6]]))))
entropy_col <- as.numeric(c(unlist(entropy_list)))
df <- data.frame(
  input = input_col,
  entropy = entropy_col
)
df$input <- as.factor(df$input)

# boxplot for harmonized Shannon Entropy
plot <- ggplot2::ggplot(df, aes(x=input, y=entropy, color=input, fill=input)) +
  ggplot2::geom_boxplot(alpha=0.5, width = 0.6,
                        outlier.shape = 21, outlier.size = 2) +
  ggplot2::scale_color_manual(values = c("OMA" =  "#f06180",
                                         "HOG [mmseqs2]" = "#54bae9",
                                         "HOG [diamond]" = "#c588eb",
                                         "OOG [mmseqs2]" = "#7bd47b",
                                         "OOG [diamond]" = "#dbdb81",
                                         "SAMap" = "#dda86c")) +
  ggplot2::scale_fill_manual(values = c("OMA" =  "#f06180",
                                        "HOG [mmseqs2]" = "#54bae9",
                                        "HOG [diamond]" = "#c588eb",
                                        "OOG [mmseqs2]" = "#7bd47b",
                                        "OOG [diamond]" = "#dbdb81",
                                        "SAMap" = "#dda86c")) +
  ggplot2::labs(title = "Shannon Entropy by Clustering Method",
                x = "Clustering Method", y = "Shannon Entropy") +
  ggplot2::geom_jitter(shape=16, position=position_jitter(0.2)) +
  ggplot2::theme_bw()
plotname <- file.path("results", "benchmark", "entropy_best.png")
ggplot2::ggsave(plotname, plot, width = 12, height = 8)