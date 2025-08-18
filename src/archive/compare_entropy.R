#'-----------------------------------------------------------------------------
#' @title   Compare Entropy
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

# -----------------------------------------------------------------------------
# Load shannon entropy of clusters generated from OrthoMAP or SAMap
#   1) Load as Rscript: Rscript src/data_analyse/compare_entropy.R [-f1] [-f2]
#     [-f1]: path to table 1
#     [-f2]: path to table 2
#   2) Load manually in interactive session
# -----------------------------------------------------------------------------
# Load library
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 3) {
  f1 <- args[1]
  f2 <- args[2]
  output_directory <- args[3]
} else if (length(args) == 1 || args[1] == "-h") {
  # apply help
  message("Usage as Rscript: Rscript compare_entropy.R [-f1] [-f2] [-o]\n",
          "   [-f1] : Path to entropy table of 1\n",
          "   [-f2] : Path to entropy table of 2\n",
          "   [-o]  : Path to output directory\n\n")
} else {
  # Load manually as interactive session
  f1 <- file.path("results", "OrthoMAP_Statistical_Results",
                  "oma", "4", "Summary_statistics.csv")
  f2 <- file.path("results", "OrthoMAP_Statistical_Results",
                  "oma", "5", "Summary_statistics.csv")
  output_directory <- "test"
}

# read as data frame
f1_entropy <- read.csv(file = f1)$Shannon_Entropy
f2_entropy <- read.csv(file = f2)$Shannon_Entropy
n1 <- length(f1_entropy)
n2 <- length(f2_entropy)
input_col <- c(rep(1, n1), rep(2, n2))
input_col <- sapply(input_col, function(x) {
  if (x == 1) res <- "euclidean"
  if (x == 2) res <- "cosine"
  res
})
entropy_col <- c(f1_entropy, f2_entropy)
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
  ggplot2::scale_color_manual(values = c("euclidean" =  "#f06180",
                                         "cosine" = "#54bae9")) +
  ggplot2::scale_fill_manual(values = c("euclidean" = "#f06180",
                                        "cosine" = "#54bae9")) +
  ggplot2::geom_vline(aes(xintercept = mean(f1_entropy)),
                      color = "#c9051f", linetype = "dashed", size = 0.5) +
  ggplot2::geom_vline(aes(xintercept = mean(f2_entropy)),
                      color = "#3905c9", linetype = "dashed", size = 0.5) +
  ggplot2::theme_bw()

plotname <- file.path(output_directory, "entropy.png")
ggplot2::ggsave(plotname, plot, width = 6, height = 4)

# boxplot for mean and median
plot <- ggplot2::ggplot(df, aes(x=input, y=entropy, color=input, fill=input)) +
  ggplot2::geom_boxplot(alpha=0.5, width = 0.6,
                        outlier.shape = 21, outlier.size = 2) +
  ggplot2::scale_color_manual(values = c("euclidean" =  "#f06180",
                                         "cosine" = "#54bae9")) +
  ggplot2::scale_fill_manual(values = c("euclidean" = "#f06180",
                                        "cosine" = "#54bae9")) +
  ggplot2::labs(title = "Shannon Entropy by Clustering Method",
                x = "Clustering Method", y = "Shannon Entropy") +
  ggplot2::geom_jitter(shape=16, position=position_jitter(0.2))
  ggplot2::theme_bw()
plotname <- file.path(output_directory, "mean_and_median.png")
ggplot2::ggsave(plotname, plot, width = 6, height = 4)