#!/usr/bin/env Rscript
library(ggpubr)

library(dplyr)
library(optparse)
library(patchwork)

# 设置命令行参数解析
option_list <- list(
  make_option(c("-a", "--stat_attrs"), type = "character", help = "postprocess.stat_attrs.tsv"),
  make_option(c("-s", "--single_file"), action = "store_true", default = FALSE, help = "Generate each plot to the single file"),
  make_option(c("-o", "--output_dir"), type = "character", help = "Output directory")
)

# 解析命令行参数
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# 创建输出目录
if (!dir.exists(opt$output_dir)) {
  dir.create(opt$output_dir)
}

stat_attrs_data <- read.csv(opt$stat_attrs, header = TRUE, sep = "\t")


#-----------------------------------------------------------------------------#
# save figure

save_basic_plots <- function(A, B, C, D, single_file, output_dir) {
  if (single_file) {
    ggsave(file.path(output_dir, "postprocess.stat_attrs_mean.pdf"), A)
    ggsave(file.path(output_dir, "postprocess.stat_attrs_min.pdf"), B)
    ggsave(file.path(output_dir, "postprocess.stat_attrs_var.pdf"), C)
    ggsave(file.path(output_dir, "postprocess.stat_attrs_uni.pdf"), D)
  } else {
    combined_plot <- A + B + C + D +
      plot_layout(
        guides = "collect",
        ncol = 2, 
        nrow = 2,
        widths = c(1, 1), 
        heights = c(1, 1)
      ) +
      theme(legend.position = "bottom")
    ggsave(file.path(output_dir, "pgap2.postprocess_stat.pdf"), combined_plot, width = 8.6, height = 7.4)
  }
}

draw_stat_attr <- function(stat_attrs_data, attr, xlab_name) {
  if (attr == "mean" || attr == "min") {
    left_position <- 0.15
  } else {
    left_position <- 0.85
  }


  attr_plot <- ggline(subset(stat_attrs_data, Attr == attr),
    x = "Edge", y = "Prop",
    color = "Group", scales = "free", size = 1,
    xlab = xlab_name, ylab = "Gene Cluster Proportion",
    # palette = c("#72B063", "#719AAC", "#B8DBB3", "#94C6CD", "#E29135")
  ) +
    scale_y_log10() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.line = element_line(linewidth = 0),
      panel.background = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), 
      panel.grid = element_blank(),
      plot.margin = margin(10, 10, 10, 10), 
      legend.box = "vertical", 
      legend.position = c(left_position, 0.85),
      legend.background = element_blank(),
      legend.direction = "vertical",
      legend.key.size = unit(0.5, "cm"), 
      legend.title = element_text(size = 10), 
      legend.text = element_text(size = 8)
    ) + scale_color_manual(
      values = c(
        "Cloud" = "#72B063",
        "Shell" = "#719AAC",
        "Soft_core" = "#B8DBB3",
        "Core" = "#94C6CD",
        "Strict_core" = "#E29135"
      ))+
    guides(color = guide_legend(override.aes = list(alpha = 1))) + ggtitle(toupper(attr))

  return(attr_plot)
}


stat_attrs_data <- stat_attrs_data %>%
  group_by(Attr, Group) %>% 
  mutate(Prop = Count / sum(Count)) %>% 
  ungroup() 
stat_attrs_data$Group <- factor(stat_attrs_data$Group, levels = c("Strict_core", "Core", "Soft_core", "Shell", "Cloud"))

A <- draw_stat_attr(stat_attrs_data, "mean", xlab_name = "Gene Identity")
B <- draw_stat_attr(stat_attrs_data, "min", xlab_name = "Gene Identity")
C <- draw_stat_attr(stat_attrs_data, "var", xlab_name = "Gene Cluster Variance")
D <- draw_stat_attr(stat_attrs_data, "uni", xlab_name = "Gene Identity")

save_basic_plots(A, B, C, D, opt$single_file, opt$output_dir)
