#'-----------------------------------------------------------------------------
#' @title   Create parameter plot
#'
#' @description
#' This script visualize the overview of the clustering parameter
#'
#' @author  David Blickenstorfer, Technau Group (2025)
#' @created 17/08/2025
#'
#' Version: v1.0.0
#' Changelog:
#'  - v1.0.0: Initial release
#'
#' @copyright Technau Group, University of Vienna, 2025
#' @copyright Boeva Lab, ETH Zurich, 2025
#' ----------------------------------------------------------------------------
# Load library
library(tidyverse)
library(patchwork)
sessionInfo()

color_vector <- c(
  "OMA" = "#f06180",
  "HOG [mmseqs2]" = "#54bae9",
  "HOG [diamond]" = "#c588eb",
  "OOG [mmseqs2]" = "#7bd47b",
  "OOG [diamond]" = "#dbdb81"
)
methodnames <- c("OMA", "HOG [mmseqs2]", "HOG [diamond]",
                 "OOG [mmseqs2]", "OOG [diamond]")
# -----------------------------------------------------------------------------
# Step 1: Load Data and Preparation
# -----------------------------------------------------------------------------
result_path <- file.path("results", "OrthoMAP_Statistical_Results")
experiment3 <- c(
  file.path(result_path, "oma", "3", "local_statistic.csv"),
  file.path(result_path, "mmseqs2", "HOG", "3", "local_statistic.csv"),
  file.path(result_path, "diamond", "HOG", "3", "local_statistic.csv"),
  file.path(result_path, "mmseqs2", "OOG", "3", "local_statistic.csv"),
  file.path(result_path, "diamond", "OOG", "3", "local_statistic.csv")
)
experiment4 <- c(
  file.path(result_path, "oma", "4", "local_statistic.csv"),
  file.path(result_path, "mmseqs2", "HOG", "4", "local_statistic.csv"),
  file.path(result_path, "diamond", "HOG", "4", "local_statistic.csv"),
  file.path(result_path, "mmseqs2", "OOG", "4", "local_statistic.csv"),
  file.path(result_path, "diamond", "OOG", "4", "local_statistic.csv")
)
experiment5 <- c(
  file.path(result_path, "oma", "5", "local_statistic.csv"),
  file.path(result_path, "mmseqs2", "HOG", "5", "local_statistic.csv"),
  file.path(result_path, "diamond", "HOG", "5", "local_statistic.csv"),
  file.path(result_path, "mmseqs2", "OOG", "5", "local_statistic.csv"),
  file.path(result_path, "diamond", "OOG", "5", "local_statistic.csv")
)
experiment6 <- c(
  file.path(result_path, "oma", "6", "local_statistic.csv"),
  file.path(result_path, "mmseqs2", "HOG", "6", "local_statistic.csv"),
  file.path(result_path, "diamond", "HOG", "6", "local_statistic.csv"),
  file.path(result_path, "mmseqs2", "OOG", "6", "local_statistic.csv"),
  file.path(result_path, "diamond", "OOG", "6", "local_statistic.csv")
)

experiment7 <- c(
  file.path(result_path, "oma", "7", "local_statistic.csv"),
  file.path(result_path, "mmseqs2", "HOG", "7", "local_statistic.csv"),
  file.path(result_path, "diamond", "HOG", "7", "local_statistic.csv"),
  file.path(result_path, "mmseqs2", "OOG", "7", "local_statistic.csv"),
  file.path(result_path, "diamond", "OOG", "7", "local_statistic.csv")
)

experiment8 <- c(
  file.path(result_path, "oma", "8", "local_statistic.csv"),
  file.path(result_path, "mmseqs2", "HOG", "8", "local_statistic.csv"),
  file.path(result_path, "diamond", "HOG", "8", "local_statistic.csv"),
  file.path(result_path, "mmseqs2", "OOG", "8", "local_statistic.csv"),
  file.path(result_path, "diamond", "OOG", "8", "local_statistic.csv")
)
# -----------------------------------------------------------------------------
# Step 2: Create Entropy Table
# -----------------------------------------------------------------------------
pca_entropy_list <- mapply(function(file, method) {
  entropy <- read.csv(file = file)$Shannon_Entropy
  method <- rep(method, length(entropy))
  data.frame(method = method, entropy = entropy)
}, experiment3, methodnames, SIMPLIFY = FALSE)
pca_table <- do.call(rbind, pca_entropy_list)
rownames(pca_table) <- NULL
pca_table$annoy.metric <- "cosine"
pca_table$reduction <- "pca"

harmony_entropy_list <- mapply(function(file, method) {
  entropy <- read.csv(file = file)$Shannon_Entropy
  method <- rep(method, length(entropy))
  data.frame(method = method, entropy = entropy)
}, experiment4, methodnames, SIMPLIFY = FALSE)
harmony_table <- do.call(rbind, harmony_entropy_list)
rownames(harmony_table) <- NULL
harmony_table$annoy.metric <- "cosine"
harmony_table$reduction <- "harmony"

euclidean_entropy_list <- mapply(function(file, method) {
  entropy <- read.csv(file = file)$Shannon_Entropy
  method <- rep(method, length(entropy))
  data.frame(method = method, entropy = entropy)
}, experiment5, methodnames, SIMPLIFY = FALSE)
euclidean_table <- do.call(rbind, euclidean_entropy_list)
rownames(euclidean_table) <- NULL
euclidean_table$annoy.metric <- "euclidean"
euclidean_table$reduction <- "pca"

res_01_entropy_list <- mapply(function(file, method) {
  entropy <- read.csv(file = file)$Shannon_Entropy
  method <- rep(method, length(entropy))
  data.frame(method = method, entropy = entropy)
}, experiment5, methodnames, SIMPLIFY = FALSE)
res_01_table <- do.call(rbind, res_01_entropy_list)
rownames(res_01_table) <- NULL
res_01_table$coarse_resolution <- "0.1"
res_01_table$fine_resolution <- "0.1"

res_02_entropy_list <- mapply(function(file, method) {
  entropy <- read.csv(file = file)$Shannon_Entropy
  method <- rep(method, length(entropy))
  data.frame(method = method, entropy = entropy)
}, experiment6, methodnames, SIMPLIFY = FALSE)
res_02_table <- do.call(rbind, res_02_entropy_list)
rownames(res_02_table) <- NULL
res_02_table$coarse_resolution <- "0.1"
res_02_table$fine_resolution <- "0.2"

res_03_entropy_list <- mapply(function(file, method) {
  entropy <- read.csv(file = file)$Shannon_Entropy
  method <- rep(method, length(entropy))
  data.frame(method = method, entropy = entropy)
}, experiment7, methodnames, SIMPLIFY = FALSE)
res_03_table <- do.call(rbind, res_03_entropy_list)
rownames(res_03_table) <- NULL
res_03_table$coarse_resolution <- "0.2"
res_03_table$fine_resolution <- "0.1"

res_04_entropy_list <- mapply(function(file, method) {
  entropy <- read.csv(file = file)$Shannon_Entropy
  method <- rep(method, length(entropy))
  data.frame(method = method, entropy = entropy)
}, experiment8, methodnames, SIMPLIFY = FALSE)
res_04_table <- do.call(rbind, res_04_entropy_list)
rownames(res_04_table) <- NULL
res_04_table$coarse_resolution <- "0.2"
res_04_table$fine_resolution <- "0.2"

df_res <- rbind(res_01_table, res_02_table, res_03_table, res_04_table)
df_res$condition <- interaction(df_res$coarse_resolution,
                                df_res$fine_resolution, sep = "-")
df_res$method <- factor(df_res$method,
                        levels = unique(df_res$method))
df_res$fine_resolution <- factor(df_res$fine_resolution,
                                 levels = unique(df_res$fine_resolution))
df_res$coarse_resolution <- factor(df_res$coarse_resolution,
                                   levels = unique(df_res$coarse_resolution))
df_res$condition <- factor(df_res$condition,
                           levels = unique(df_res$condition))

# boxplot for mean and median
plot_res <- ggplot2::ggplot(df_res, aes(x = condition, y = entropy,
                                        color = method, fill = method)) +
  ggplot2::geom_boxplot(alpha=0.5, width = 0.6,
                        outlier.shape = 21, outlier.size = 2) +
  ggplot2::scale_color_manual(values = color_vector) +
  ggplot2::scale_fill_manual(values = color_vector) +
  ggplot2::labs(y = "Harmonized Shannon Entropy") +
  ggplot2::geom_jitter(shape = 16, position = position_jitter(0.2)) +
  ggplot2::theme_bw() + ggplot2::theme(axis.title.x = ggplot2::element_blank())

df_con <- rbind(pca_table, harmony_table, euclidean_table)
df_con$condition <- interaction(df_con$reduction,
                                df_con$annoy.metric, sep = "-")
df_con$method <- factor(df_con$method,
                        levels = unique(df_con$method))
df_con$reduction <- factor(df_con$reduction,
                           levels = unique(df_con$reduction))
df_con$annoy.metric <- factor(df_con$annoy.metric,
                              levels = unique(df_con$annoy.metric))
df_con$condition <- factor(df_con$condition,
                           levels = unique(df_con$condition))

# boxplot for mean and median
plot_con <- ggplot2::ggplot(df_con, aes(x = condition, y = entropy,
                                        color = method, fill = method)) +
  ggplot2::geom_boxplot(alpha = 0.5, width = 0.6,
                        outlier.shape = 21, outlier.size = 2) +
  ggplot2::scale_color_manual(values = color_vector) +
  ggplot2::scale_fill_manual(values = color_vector) +
  ggplot2::labs(y = "Harmonized Shannon Entropy") +
  ggplot2::geom_jitter(shape = 16, position = position_jitter(0.2)) +
  ggplot2::theme_bw() + ggplot2::theme(axis.title.x = ggplot2::element_blank())

plot_final <- plot_con / plot_res +
  patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(tag_levels = "A") &
  theme(
    legend.position = "right",
    title = element_blank(),
    plot.tag = element_text(size = 20, face = "bold")
  )
plotname <- file.path("results", "benchmark", "configuration.pdf")
ggplot2::ggsave(plotname, plot_final, width = 12, height = 10)
plotname <- file.path("results", "benchmark", "configuration_entropy.pdf")
ggplot2::ggsave(plotname, plot_con, width = 12, height = 8)
plotname <- file.path("results", "benchmark", "resolution_entropy.pdf")
ggplot2::ggsave(plotname, plot_res, width = 12, height = 8)