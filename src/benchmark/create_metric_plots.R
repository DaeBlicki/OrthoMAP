#'-----------------------------------------------------------------------------
#' @title   Run Benchmark for OrthoMAP and SAMap
#'
#' @concept
#' Estimates the global clustering in form of Adjusted Random Index (ARI) and
#' Normalized Mutual Information (NMI), and local clustering in form of harm.
#' Shannon entropy of each cluster. The results are in form of dotplots.
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
library(patchwork)
color_vector <- c(
  "OMA" = "#f06180",
  "HOG [mmseqs2]" = "#54bae9",
  "HOG [diamond]" = "#c588eb",
  "OOG [mmseqs2]" = "#7bd47b",
  "OOG [diamond]" = "#dbdb81",
  "SAMap" = "#dda86c"
)
methodnames <- c("OMA", "HOG [mmseqs2]", "HOG [diamond]",
                 "OOG [mmseqs2]", "OOG [diamond]", "SAMap")
methodnames_no_samap <- c("OMA", "HOG [mmseqs2]", "HOG [diamond]",
                          "OOG [mmseqs2]", "OOG [diamond]")
# -----------------------------------------------------------------------------
# Step 1: Load Data and Preparation
# -----------------------------------------------------------------------------
metric_path <- file.path("results", "benchmark", "metrics")
ari_nmi_files <- c(
  file.path(metric_path, "oma_ari_nmi.csv"),
  file.path(metric_path, "mmseqs2_hog_ari_nmi.csv"),
  file.path(metric_path, "diamond_hog_ari_nmi.csv"),
  file.path(metric_path, "mmseqs2_oog_ari_nmi.csv"),
  file.path(metric_path, "diamond_oog_ari_nmi.csv"),
  file.path(metric_path, "samap_ari_nmi.csv")
)
entropy_files <- c(
  file.path(metric_path, "oma_entropy.csv"),
  file.path(metric_path, "mmseqs2_hog_entropy.csv"),
  file.path(metric_path, "mmseqs2_oog_entropy.csv"),
  file.path(metric_path, "diamond_hog_entropy.csv"),
  file.path(metric_path, "diamond_oog_entropy.csv"),
  file.path(metric_path, "samap_entropy.csv")
)
asw_lisi_files <- c(
  file.path(metric_path, "oma_asw_lisi.csv"),
  file.path(metric_path, "mmseqs2_hog_asw_lisi.csv"),
  file.path(metric_path, "mmseqs2_oog_asw_lisi.csv"),
  file.path(metric_path, "diamond_hog_asw_lisi.csv"),
  file.path(metric_path, "diamond_oog_asw_lisi.csv")
)
# -----------------------------------------------------------------------------
# Step 2: Process batch and bio conservation to get statistic table
# -----------------------------------------------------------------------------
ari_nmi_list <- lapply(ari_nmi_files, function(file) {
  raw_table <- read.csv(file, row.names = 1)
  ari_bio_mean <- mean(as.numeric(raw_table["cARI", ]))
  ari_bio_sd <- sd(as.numeric(raw_table["cARI", ]))
  ari_batch_mean <- mean(as.numeric(raw_table["iARI", ]))
  ari_batch_sd <- sd(as.numeric(raw_table["iARI", ]))
  nmi_bio_mean <- mean(as.numeric(raw_table["cNMI", ]))
  nmi_bio_sd <- sd(as.numeric(raw_table["cNMI", ]))
  nmi_batch_mean <- mean(as.numeric(raw_table["iNMI", ]))
  nmi_batch_sd <- sd(as.numeric(raw_table["iNMI", ]))
  names <- c("cARI_mean", "cARI_sd", "iARI_mean", "iARI_sd",
             "cNMI_mean", "cNMI_sd", "iNMI_mean", "iNMI_sd")
  metrics <- c(ari_bio_mean, ari_bio_sd, ari_batch_mean, ari_batch_sd,
               nmi_bio_mean, nmi_bio_sd, nmi_batch_mean, nmi_batch_sd)
  names(metrics) <- names
  metrics
})
ari_nmi_table <- as.data.frame(do.call(rbind, ari_nmi_list))
rownames(ari_nmi_table) <- methodnames

asw_lisi_list <- lapply(asw_lisi_files, function(file) {
  raw_table <- read.csv(file, row.names = 1)
  asw_bio_mean <- mean(as.numeric(raw_table[, "cASW"]))
  asw_bio_sd <- sd(as.numeric(raw_table[, "cASW"]))
  asw_batch_mean <- mean(as.numeric(raw_table[, "iASW"]))
  asw_batch_sd <- sd(as.numeric(raw_table[, "iASW"]))
  lisi_bio_mean <- mean(as.numeric(raw_table[, "cLISI"]))
  lisi_bio_sd <- sd(as.numeric(raw_table[, "cLISI"]))
  lisi_batch_mean <- mean(as.numeric(raw_table[, "iLISI"]))
  lisi_batch_sd <- sd(as.numeric(raw_table[, "iLISI"]))
  names <- c("cASW_mean", "cASW_sd", "iASW_mean", "iASW_sd",
             "cLISI_mean", "cLISI_sd", "iLISI_mean", "iLISI_sd")
  metrics <- c(asw_bio_mean, asw_bio_sd, asw_batch_mean, asw_batch_sd,
               lisi_bio_mean, lisi_bio_sd, lisi_batch_mean, lisi_batch_sd)
  names(metrics) <- names
  metrics
})
asw_lisi_table <- as.data.frame(do.call(rbind, asw_lisi_list))
rownames(asw_lisi_table) <- methodnames_no_samap

metric_table <- merge(ari_nmi_table, asw_lisi_table,
                      by = "row.names", all = TRUE)
rownames(metric_table) <- metric_table[, "Row.names"]
metric_table$Row.names <- NULL
filename <- file.path("results", "benchmark", "conservation_table.csv")
write.csv(metric_table, filename, row.names = TRUE)

metric_table$methods <- rownames(metric_table)
# -----------------------------------------------------------------------------
# Step 3: Merge the results in the Entropy table
# -----------------------------------------------------------------------------
entropy_list <- lapply(entropy_files, function(file) {
  df <- read.csv(file)
  df$entropy
})

# -----------------------------------------------------------------------------
# Step 4: Create Lollipop plot
# -----------------------------------------------------------------------------
samap_data <- which(metric_table$methods == "SAMap")
ari_nmi_label <- ggplot(metric_table,
                        aes(y = methods, x = 1, label = methods)) +
  geom_text(size = 3, hjust = 0.5) +  # centered horizontally
  theme_void() +
  theme(
    plot.margin = margin(0, 0, 0, 0),    # remove plot margins
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

asw_lisi_label <- ggplot(metric_table[-samap_data, ],
                        aes(y = methods, x = 1, label = methods)) +
  geom_text(size = 3, hjust = 0.5) +  # centered horizontally
  theme_void() +
  theme(
    plot.margin = margin(0, 0, 0, 0),    # remove plot margins
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

# create batch plot
create_batch_plot <- function(metric_name,
                              segment_color = "black",
                              point_color = "black") {
  clean_metric_name <- sub("_mean$", "", metric_name)
  df <- metric_table
  if (as.character(metric_name) == "iASW_mean" ||
      as.character(metric_name) == "iLISI_mean") {
    df <- metric_table[-samap_data, ]
  }
  mean <- mean(as.numeric(df[, metric_name]))
  plot <- ggplot2::ggplot(df,
                          aes_string(x = "methods", y = metric_name)) +
    ggplot2::geom_segment(aes_string(xend = "methods", yend = 0),
                          color = segment_color) +
    ggplot2::geom_point(color = point_color, size = 4) +
    ggplot2::scale_x_discrete(position = "top") +
    ggplot2::theme_bw() +
    ggplot2::coord_flip() +
    ggplot2::theme(
      plot.margin = ggplot2::margin(0, 0, 0, 0),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(colour = "black")
    ) +
    ggplot2::scale_y_reverse(limits = c(1, 0),
                             breaks = seq(0, 1, by = 0.1)) +
    ggplot2::labs(y = clean_metric_name) +
    ggplot2::geom_hline(yintercept = mean,
                        linetype = "dashed",
                        color = point_color,
                        size = 0.3)
  plot
}


create_bio_plot <- function(metric_name,
                            segment_color = "black",
                            point_color = "black") {
  clean_metric_name <- sub("_mean$", "", metric_name)
  df <- metric_table
  if (as.character(metric_name) == "cASW_mean" ||
      as.character(metric_name) == "cLISI_mean") {
    df <- metric_table[-samap_data, ]
  }
  mean <- mean(as.numeric(df[, metric_name]))
  plot <- ggplot2::ggplot(df,
                          aes_string(x = "methods", y = metric_name)) +
    ggplot2::geom_segment(aes_string(xend = "methods", yend = 0),
                          color = segment_color) +
    ggplot2::geom_point(color = point_color, size = 4) +
    ggplot2::theme_bw() +
    ggplot2::coord_flip() +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.margin = ggplot2::margin(0, 0, 0, 0),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(colour = "black")
    ) +
    ggplot2::scale_y_continuous(limits = c(0, 1),
                                breaks = seq(0, 1, by = 0.1)) +
    ggplot2::labs(y = clean_metric_name) +
    ggplot2::geom_hline(yintercept = mean,
                        linetype = "dashed",
                        color = point_color,
                        size = 0.3)
  plot
}
# Create plots for ARI and NMI
ari_nmi_plot <- create_batch_plot("iNMI_mean", point_color = "#00c510") +
  create_batch_plot("iARI_mean", point_color = "#af0a2e") +
  ari_nmi_label +
  create_bio_plot("cARI_mean", point_color = "#af0a2e") +
  create_bio_plot("cNMI_mean", point_color = "#00c510") +
  patchwork::plot_layout(ncol = 5, width = c(2, 2, 1, 2, 2)) +
  ggplot2::theme(plot.margin = margin(0, 0, 0, 0))

plotname <- file.path("results", "benchmark", "ari_nmi.pdf")
ggplot2::ggsave(plotname, ari_nmi_plot, width = 15, height = 4)

# Create plots for ASW and LISI
asw_lisi_plot <- create_batch_plot("iLISI_mean", point_color = "#ecb500") +
  create_batch_plot("iASW_mean", point_color = "#0724ca") +
  asw_lisi_label +
  create_bio_plot("cASW_mean", point_color = "#0724ca") +
  create_bio_plot("cLISI_mean", point_color = "#ecb500") +
  patchwork::plot_layout(ncol = 5, width = c(2, 2, 1, 2, 2)) +
  ggplot2::theme(plot.margin = margin(0, 0, 0, 0))

plotname <- file.path("results", "benchmark", "asw_lisi.pdf")
ggplot2::ggsave(plotname, asw_lisi_plot, width = 15, height = 4)

# -----------------------------------------------------------------------------
# Step 5: Estimation of normalized Shannon Entropy
# -----------------------------------------------------------------------------
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
  ggplot2::scale_color_manual(values = color_vector) +
  ggplot2::scale_fill_manual(values = color_vector) +
  ggplot2::labs(y = "Harmonized Shannon Entropy") +
  ggplot2::geom_jitter(shape=16, position=position_jitter(0.2)) +
  ggplot2::theme_bw() + ggplot2::theme(axis.title.x = ggplot2::element_blank())


plotname <- file.path("results", "benchmark", "entropy.pdf")
ggplot2::ggsave(plotname, plot, width = 12, height = 6)