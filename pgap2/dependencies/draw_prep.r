#!/usr/bin/env Rscript

library(ggpubr)
library(ggrepel)

library(dplyr)
library(tidyr)
library(optparse)
library(patchwork)

option_list <- list(
  make_option(c("-a", "--stat_attrs"), type = "character", help = "preprocess.stat.tsv"),
  make_option(c("-b", "--gene_code"), type = "character", help = "preprocess.gene_code.csv"),
  make_option(c("-c", "--ani_thre"), type = "character", help = "ANI threshold"),
  make_option(c("-s", "--single_file"), action = "store_true", default = FALSE, help = "Generate each plot to the single file"),
  make_option(c("-o", "--output_dir"), type = "character", help = "Output directory")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (!dir.exists(opt$output_dir)) {
  dir.create(opt$output_dir)
}

input_prep <- read.csv(opt$stat_attrs, header = TRUE, sep = "\t")
input_code <- read.csv(opt$gene_code, header = TRUE, sep = ",", row.names = 1)


#-----------------------------------------------------------------------------#

save_basic_plots <- function(A, B, C, D, E, single_file, output_dir) {
  if (single_file) {
    ggsave(file.path(output_dir, "preprocess.ANI.pdf"), A)
    ggsave(file.path(output_dir, "preprocess.gene_number.pdf"), B)
    ggsave(file.path(output_dir, "preprocess.half_core.pdf"), C)
    ggsave(file.path(output_dir, "preprocess.proportion.pdf"), D)
    ggsave(file.path(output_dir, "preprocess.gene_code.pdf"), E)
  } else {
    p1 <- A + ggtitle("A") +
      B + ggtitle("B") +
      C + ggtitle("C") +
      plot_layout(guides = "collect", ncol = 3)
    p2 <- D + ggtitle("D") +
      E + ggtitle("E") +
      plot_layout(guides = "collect", ncol = 2) & 
      theme(legend.position = "right")
    combined_plot <- p1 / p2 + plot_layout(heights = c(1, 2))
    ggsave(file.path(output_dir, "pgap2.preprocess.pdf"), combined_plot, width = 9.8, height = 8.4)
  }
}

draw_ani <- function(input_prep, ani_thre) {
  p <- ggviolin(input_prep,
    y = "ani", add = "jitter",
    add.params = list(size = 3, color = "black", fill = "#B8DBB3", shape = 21),
    xlab = "Genome", ylab = "ANI"
  ) +
    geom_hline(yintercept = ani_thre, linetype = "dotted", col = "red") +
    theme(legend.position = "none") +
    scale_x_discrete(labels = c("")) +
    geom_text_repel(
      data = subset(input_prep, ani == 100 | is_outlier_ani == 1), 
      aes(x = as.factor(1), y = ani, label = strain), 
      arrow = arrow(type = "closed", length = unit(0.1, "inches")), 
      vjust = -0.5, size = 3, color = "black", 
      box.padding = 0.2, min.segment.length = 1 
    )
  return(p)
}

draw_gene_number <- function(input_prep) {
  p <- ggviolin(input_prep,
    y = "total_gene_num", add = "jitter", xlab = "Genome", ylab = "Number of genes",
    add.params = list(size = 3, color = "black", fill = "#B8DBB3", shape = 21)
  ) + theme(legend.position = "None") + scale_x_discrete(labels = c(""))
  return(p)
}

draw_half_core <- function(input_prep) {
  p <- ggscatter(input_prep,
    x = "half_core", y = "single_cloud",
    size = "total_gene_num", color = "total_gene_num", xlab = "Half core", ylab = "Single cloud"
  ) +
    guides(size = "none", color = guide_colorbar(barheight = 5)) +
    scale_color_gradient(high = "#990033", low = "white") +
    theme(legend.position = "right") +
    labs(color = "Count")
  return(p)
}

draw_proportion <- function(input_prep) {
  input_prep$others <- input_prep$total_gene_num - input_prep$half_core - input_prep$single_cloud
  used_prep <- input_prep %>% select(strain, total_gene_num, half_core, single_cloud, others)
  # Transform the data into a long format
  long_data <- used_prep %>%
    pivot_longer(cols = c(half_core, single_cloud, others), names_to = "Category", values_to = "Count") %>%
    mutate(
      Category = factor(Category, levels = c("single_cloud", "others", "half_core")), # Reorder the categories
      strain = factor(strain, levels = used_prep$strain)
    ) # Maintain the order of strains

  p <- ggplot(long_data, aes(x = strain, y = Count, fill = Category)) +
    geom_bar(position = "fill", stat = "identity") +
    coord_flip() +
    theme_minimal() +
    labs(x = "Strains", y = "Proportion", fill = "Category") +
    scale_fill_manual(values = c("single_cloud" = "#72B063", "others" = "#719AAC", "half_core" = "#E29135"), labels = c("Single", "Others", "Half")) +
    scale_x_discrete(
      breaks = function(x) {
        y_levels <- levels(factor(long_data$strain))
        if (length(y_levels) > 25) {
          return(y_levels[seq(1, length(y_levels), length.out = 25)]) 
        } else {
          return(y_levels) 
        }
      }
    )
  return(p)
}

draw_gene_code <- function(input_code) {
  df <- as.data.frame(as.table(as.matrix(input_code)))
  df <- df %>%
    rename(Count = Freq)
  df$Count <- as.numeric(as.character(df$Count))
  p <- ggplot(df, aes(x = Var1, y = Var2, fill = Count)) +
    geom_tile() +
    scale_fill_gradient2(high = "#990033", low = "white") +
    theme_bw() +
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5), ) +
    xlab("Start|Stop codon") +
    ylab(NULL) +
    scale_y_discrete(
      breaks = function(x) {
        y_levels <- levels(factor(df$Var2))
        if (length(y_levels) > 25) {
          return(y_levels[seq(1, length(y_levels), length.out = 25)])
        } else {
          return(y_levels) 
        }
      }
    )
  return(p)
}

A <- draw_ani(input_prep, as.numeric(opt$ani_thre))
B <- draw_gene_number(input_prep)
C <- draw_half_core(input_prep)
D <- draw_proportion(input_prep)
E <- draw_gene_code(input_code)
save_basic_plots(A, B, C, D, E, opt$single_file, opt$output_dir)
