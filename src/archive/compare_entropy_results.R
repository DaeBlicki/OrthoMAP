#'-----------------------------------------------------------------------------
#' @title   Compare Entropy for OMA and OrthoFinder
#'
#' @concept
#' Compares the harmonized shannon entropy of the estimated clusters generated
#' by the OrthoMAP (and SAMap) workflow. The results must be in form of `csv`
#' and must have column `Shannon_Entropy`.
#'
#' @description
#' It calculates the mean and standard deviation of the Shannon entropy.
#' Furthermore, it plots an histogram distribution and it approximated to normal
#' distribution under the assumption, that the histogram represents a binomial
#' distribution.
#'
#' @author  David Blickenstorfer, Technau Group (2025)
#' @created 04/08/2025
#'
#' Version: v1.0.0
#' Changelog:
#'  - v1.0.0: Initial release
#'
#' @copyright Technau Group, University of Vienna, 2025
#' @copyright Boeva Lab, ETH Zurich, 2025
#' ----------------------------------------------------------------------------
library(tidyverse)
source("src/data_analyse/randomized_entropy_result.R")

# Load manually as interactive session
f1 <- file.path("results", "OrthoMAP_Statistical_Results",
                "oma", "5", "Summary_statistics.csv")
f2 <- file.path("results", "OrthoMAP_Statistical_Results",
                "mmseqs2", "HOG", "5", "Summary_statistics.csv")
f3 <- file.path("results", "OrthoMAP_Statistical_Results",
                "diamond", "HOG", "5", "Summary_statistics.csv")
f4 <- file.path("results", "OrthoMAP_Statistical_Results",
                "mmseqs2", "OOG", "5", "Summary_statistics.csv")
f5 <- file.path("results", "OrthoMAP_Statistical_Results",
                "diamond", "OOG", "5", "Summary_statistics.csv")

# generate output directory
output_directory <- file.path("results", "thesis")
dir.create(output_directory,
           recursive = TRUE,
           showWarnings = FALSE)

## create randomized control
# 1: using exponential cell distribution and uniform clustering
if (FALSE) {
  control <- randomized_entropy_result(num_cells = 150000,
                                       num_tissues = 29,
                                       num_clusters = 30,
                                       rng_seed = 42)
  ncontrol <- length(control)
}
# 2: using orthomap object to get actual cell and cluster distribution
if (TRUE) {
  orthomap_obj_path <- file.path("results", "OrthoMAP_Seurat_Objects",
                                 "mmseqs2", "OOG", "result2.Robj")
  orthomap_obj <- get(load(orthomap_obj_path))
  control <- simulate_entropy_result(orthomap_obj, rng_seed = 42)
  ncontrol <- length(control)
}

# read as data frame
f1_entropy <- read.csv(file = f1)$Shannon_Entropy
f2_entropy <- read.csv(file = f2)$Shannon_Entropy
f3_entropy <- read.csv(file = f3)$Shannon_Entropy
f4_entropy <- read.csv(file = f4)$Shannon_Entropy
f5_entropy <- read.csv(file = f5)$Shannon_Entropy
n1 <- length(f1_entropy)
n2 <- length(f2_entropy)
n3 <- length(f3_entropy)
n4 <- length(f4_entropy)
n5 <- length(f5_entropy)
input_col <- c(rep("OMA", n1),
               rep("HOG [mmseqs2]", n2),
               rep("HOG [diamond]", n3),
               rep("OOG [mmseqs2]", n4),
               rep("OOG [diamond]", n5),
               rep("Randomized", ncontrol))

entropy_col <- c(f1_entropy, f2_entropy, f3_entropy,
                 f4_entropy, f5_entropy, control)
df <- data.frame(
  input = input_col,
  entropy = entropy_col
)
df$input <- as.factor(df$input)

# data analyse
mu <- plyr::ddply(df, "input", summarise, grp.mean = mean(entropy))

# store histogram
plot <- ggplot2::ggplot(df, aes(x=entropy, color=input, fill=input)) +
  ggplot2::geom_histogram(position="identity", alpha=0.5, binwidth = 0.01) +
  ggplot2::scale_color_manual(values = c("OMA" =  "#f06180",
                                         "HOG [mmseqs2]" = "#54bae9",
                                         "HOG [diamond]" = "#c588eb",
                                         "OOG [mmseqs2]" = "#7bd47b",
                                         "OOG [diamond]" = "#dbdb81",
                                         "Randomized" = "#acacac")) +
  ggplot2::scale_fill_manual(values = c("OMA" =  "#f06180",
                                        "HOG [mmseqs2]" = "#54bae9",
                                        "HOG [diamond]" = "#c588eb",
                                        "OOG [mmseqs2]" = "#7bd47b",
                                        "OOG [diamond]" = "#dbdb81",
                                        "Randomized" = "#acacac")) +
  ggplot2::geom_vline(aes(xintercept = mean(f1_entropy)),
                      color = "#c9051f", linetype = "dashed", size = 0.5) +
  ggplot2::geom_vline(aes(xintercept = mean(f2_entropy)),
                      color = "#3905c9", linetype = "dashed", size = 0.5) +
  ggplot2::theme_bw()

plotname <- file.path(output_directory, "histogram.png")
ggplot2::ggsave(plotname, plot, width = 12, height = 8)

# boxplot for mean and median
plot <- ggplot2::ggplot(df, aes(x=input, y=entropy, color=input, fill=input)) +
  ggplot2::geom_boxplot(alpha=0.5, width = 0.6,
                        outlier.shape = 21, outlier.size = 2) +
  ggplot2::scale_color_manual(values = c("OMA" =  "#f06180",
                                         "HOG [mmseqs2]" = "#54bae9",
                                         "HOG [diamond]" = "#c588eb",
                                         "OOG [mmseqs2]" = "#7bd47b",
                                         "OOG [diamond]" = "#dbdb81",
                                         "Randomized" = "#acacac")) +
  ggplot2::scale_fill_manual(values = c("OMA" =  "#f06180",
                                        "HOG [mmseqs2]" = "#54bae9",
                                        "HOG [diamond]" = "#c588eb",
                                        "OOG [mmseqs2]" = "#7bd47b",
                                        "OOG [diamond]" = "#dbdb81",
                                        "Randomized" = "#acacac")) +
  ggplot2::labs(title = "Shannon Entropy by Clustering Method",
                x = "Clustering Method", y = "Shannon Entropy") +
  ggplot2::geom_jitter(shape=16, position=position_jitter(0.2)) +
  ggplot2::theme_bw()
plotname <- file.path(output_directory, "boxplot.png")
ggplot2::ggsave(plotname, plot, width = 12, height = 8)