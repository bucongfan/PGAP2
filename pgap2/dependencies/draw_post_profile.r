#!/usr/bin/env Rscript
# install.packages(c("ggplot2", "dplyr", "optparse", "readr", "scales", "patchwork"))
library(ggpubr)

library(dplyr)
library(optparse)
library(patchwork)

option_list <- list(
  make_option(c("-a", "--pan_group_stat"), type = "character", help = "postprocess.pan_group_stat.tsv"),
  make_option(c("-b", "--clust_strain_freq"), type = "character", help = "postprocess.clust_strain_freq.tsv"),
  make_option(c("-c", "--rarefaction"), type = "character", help = "postprocess.rarefaction.tsv"),
  make_option(c("-d", "--new_clusters"), type = "character", help = "postprocess.new_clusters.tsv"),
  make_option(c("-e", "--para_stat"), type = "character", help = "postprocess.para_stat.tsv"),
  make_option(c("-s", "--single_file"), type = "logical", action = "store_true", default = FALSE, help = "Generate each plot to the single file"),
  make_option(c("-o", "--output_dir"), type = "character", help = "Output directory")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (!dir.exists(opt$output_dir)) {
  dir.create(opt$output_dir)
}


pan_group_stat_data <- read.csv(opt$pan_group_stat, header = TRUE, sep = "\t")

clust_strain_freq_data <- read.csv(opt$clust_strain_freq, header = TRUE, sep = "\t")

rarefaction_data <- read.csv(opt$rarefaction, header = TRUE, sep = "\t")

new_clusters_data <- read.csv(opt$new_clusters, header = TRUE, sep = "\t")

para_stat_data <- read.csv(opt$para_stat, header = TRUE, sep = "\t")

#-----------------------------------------------------------------------------#
save_para_plots <- function(A, B, single_file, output_dir) {
  if (single_file) {
    ggsave(file.path(output_dir, "postprocess.para_stat.pdf"), A + ggtitle(""))
    ggsave(file.path(output_dir, "postprocess.para_stat_facet.pdf"), B + ggtitle(""))
  } else {
    combined_plot <- A + B +
      plot_layout(
        ncol = 2,
        nrow = 1,
        widths = c(1, 1),
        heights = c(1, 1)
      )
    ggsave(file.path(output_dir, "pgap2.postprocess_stat_para.pdf"), combined_plot, width = 10.52, height = 5.8)
  }
}

save_basic_plots <- function(A, B, C, D, single_file, output_dir) {
  if (single_file) {
    ggsave(file.path(output_dir, "postprocess.pan_group_stat.pdf"), A)
    ggsave(file.path(output_dir, "postprocess.clust_strain_freq.pdf"), B)
    ggsave(file.path(output_dir, "postprocess.rarefaction.pdf"), C)
    ggsave(file.path(output_dir, "postprocess.new_clusters.pdf"), D)
  } else {
    combined_plot <- A + B + C + D +
      plot_layout(
        ncol = 2,
        nrow = 2,
        widths = c(1, 1),
        heights = c(1, 1)
      )
    ggsave(file.path(output_dir, "pgap2.postprocess_profile.pdf"), combined_plot, width = 9.52, height = 7.25)
  }
}

draw_pan_group_stat <- function(pan_group_stat_data) {
  labs <- paste0(pan_group_stat_data$Group, ": ", pan_group_stat_data$Count, " (", pan_group_stat_data$Proportion, "%)")

  pan_group_stat <- ggpie(pan_group_stat_data, "Count",
    label = labs,
    fill = "Group", color = "white", lab.pos = "in", size = 1,
    palette = c("#B8DBB3", "#72B063", "#719AAC", "#E29135", "#94C6CD"),
    lab.font = c(4, "plain", "black"),
  ) + theme(
    legend.position = "none",
    plot.margin = margin(0, 0, 0, 0),
    axis.text = element_blank(),
    axis.title = element_blank()
  ) + ggtitle("A")
  return(pan_group_stat)
}

draw_clust_strain_freq <- function(clust_strain_freq_data) {
  clust_strain_freq <- gghistogram(clust_strain_freq_data,
    x = "Freq", bins = 50, fill = "#B8DBB3",
    xlab = "Genome Frequency", ylab = "Gene Cluster Number"
  ) + scale_y_log10() + ggtitle("B")
  return(clust_strain_freq)
}


draw_rarefaction <- function(rarefaction_data) {
  rarefaction <- ggscatter(rarefaction_data,
    x = "Strain", y = "Sampling", color = "Type", fill = NA, alpha = 0.5, size = 1, palette = c("#72B063", "#719AAC"),
    xlab = "Genome Number", ylab = "Gene Cluster Number"
  ) + theme(
    axis.line = element_line(linewidth = 0),
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid = element_blank(),
    plot.margin = margin(10, 10, 10, 10),
    legend.position = c(0.15, 0.85),
    legend.box = "horizontal",
    legend.background = element_blank(),
    legend.key.size = unit(0.5, "cm"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  ) + guides(color = guide_legend(override.aes = list(alpha = 1))) + ggtitle("C")
  # ggboxplot(rarefaction.data, color = 'Type', x = 'Strain', y = 'Sampling',fill=NA) +
  #  theme(
  #    axis.line = element_line(linewidth = 0),
  #    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.75),
  #    panel.grid = element_blank(),
  #    plot.margin = margin(10, 10, 10, 10)
  #  )+
  #  scale_x_discrete(breaks = unique(rarefaction.data$Strain)[seq(1, length(unique(rarefaction.data$Strain)), by = 5)])

  return(rarefaction)
}

draw_new_clusters <- function(new_clusters_data) {
  new_clusters <- ggbarplot(new_clusters_data,
    x = "Strain", y = "Sampling", add = "mean_se", fill = "#94C6CD",
    xlab = "Genome Number", ylab = "New Gene Cluster Number"
  ) + theme(
    axis.line = element_line(linewidth = 0),
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid = element_blank(),
    plot.margin = margin(10, 10, 10, 10),
    legend.box = "horizontal",
    legend.key.size = unit(0.5, "cm"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  ) + guides(color = guide_legend(override.aes = list(alpha = 1))) + ggtitle("D")
  return(new_clusters)
}

draw_para_stat <- function(para_stat_data) {
  para_stat_data$duplication_degree <- log(para_stat_data$Para_gene / para_stat_data$Para_strain)
  para_stat_a <- ggscatter(para_stat_data,
    x = "Para_strain", y = "Para_gene", size = "duplication_degree",
    color = "black", shape = 21, fill = "Group",
    conf.int = TRUE,
    palette = c("#72B063", "#719AAC", "#B8DBB3", "#94C6CD", "#E29135"),
    xlab = "Paralogous Strain",
    ylab = "Paralogous Gene"
  ) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    guides(
      size = "none"
    ) + scale_fill_manual(
      values = c(
        "Cloud" = "#72B063",
        "Shell" = "#719AAC",
        "Soft_core" = "#B8DBB3",
        "Core" = "#94C6CD",
        "Strict_core" = "#E29135"
      )
    ) + ggtitle("A")
  return(para_stat_a)
}

draw_para_stat_facet <- function(para_stat_data) {
  para_stat_data$duplication_degree <- log(para_stat_data$Para_gene / para_stat_data$Para_strain)
  para_stat_b <- ggscatter(para_stat_data,
    x = "Para_strain", y = "Para_gene", size = "duplication_degree", color = "duplication_degree",
    conf.int = TRUE, facet.by = "Group",
    xlab = "Paralogous Strain",
    ylab = "Paralogous Gene"
  ) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    scale_y_log10() + scale_x_log10() +
    scale_color_gradient(low = "#dddfe6", high = "#d90429") +
    labs(color = "Dup degree") +
    guides(
      size = "none"
    ) + ggtitle("B")
  return(para_stat_b)
}

A <- draw_pan_group_stat(pan_group_stat_data)
B <- draw_clust_strain_freq(clust_strain_freq_data)
C <- draw_rarefaction(rarefaction_data)
D <- draw_new_clusters(new_clusters_data)

save_basic_plots(A, B, C, D, opt$single_file, opt$output_dir)
# continue if para_stat_data is not empty
if (nrow(para_stat_data) != 0) {
  E <- draw_para_stat(para_stat_data)
  F <- draw_para_stat_facet(para_stat_data)
  save_para_plots(E, F, opt$single_file, opt$output_dir)
}
