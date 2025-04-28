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
# 保存图形

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
        ncol = 2, # 2列
        nrow = 2, # 2行
        widths = c(1, 1), # 第一列宽度为第二列的两倍
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
    palette = c("#B8DBB3", "#72B063", "#719AAC", "#E29135", "#94C6CD")
  ) +
    scale_y_log10() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.line = element_line(linewidth = 0),
      panel.background = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # 添加边框
      panel.grid = element_blank(), # 去除网格线
      plot.margin = margin(10, 10, 10, 10), # 增加图形的外边距，防止文字被裁剪
      legend.box = "vertical", # 使图例在框内水平排列
      legend.position = c(left_position, 0.85),
      legend.background = element_blank(),
      legend.direction = "vertical",
      legend.key.size = unit(0.5, "cm"), # 控制图例块的大小
      legend.title = element_text(size = 10), # 设置图例标题字体大小
      legend.text = element_text(size = 8)
    ) +
    guides(color = guide_legend(override.aes = list(alpha = 1))) + ggtitle(toupper(attr))

  return(attr_plot)
}


stat_attrs_data <- stat_attrs_data %>%
  group_by(Attr, Group) %>% # 按 Attr 分组
  mutate(Prop = Count / sum(Count)) %>% # 计算每行的 Prop 值
  ungroup() # 取消分组
stat_attrs_data$Group <- factor(stat_attrs_data$Group, levels = c("Strict_core", "Core", "Soft_core", "Shell", "Cloud"))

A <- draw_stat_attr(stat_attrs_data, "mean", xlab_name = "Gene Identity")
B <- draw_stat_attr(stat_attrs_data, "min", xlab_name = "Gene Identity")
C <- draw_stat_attr(stat_attrs_data, "var", xlab_name = "Gene Cluster Variance")
D <- draw_stat_attr(stat_attrs_data, "uni", xlab_name = "Gene Identity")

save_basic_plots(A, B, C, D, opt$single_file, opt$output_dir)
