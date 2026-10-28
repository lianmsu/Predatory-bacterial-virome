# =============================================================================
# Predatory-bacterial-virome -- phage annotation / Gene dot-plot, simple variant
# =============================================================================
# Simplified variant of the gene dot-plot, driven by embedded example data;
# useful as a template for panel styling.
#
# Manuscript      : Fig. S6 (variant)
# Main outputs    : virus_protein_dotplot.pdf
# Original script : 1mobilome_predatory_bacteria/7phage_gene_annotation/phage_gene_anno_plot/virus_dotplot_simple.R
# -----------------------------------------------------------------------------
# Input data and output files are resolved from GTR_DATA (current analysis root),
# GTR_WORK (earlier working copy) and GTR_DEF (shared defensome data). All three
# default to the current working directory; set them instead of editing files:
#
#     Sys.setenv(GTR_DATA = "/path/to/06_predatory_bacteria_virome")
#
# =============================================================================
GTR_DATA <- normalizePath(Sys.getenv("GTR_DATA", unset = "."), mustWork = FALSE)
GTR_DEF <- normalizePath(Sys.getenv("GTR_DEF", unset = "."), mustWork = FALSE)

# 古菌病毒蛋白质特征点图 - 简化版本
library(ggplot2)

# 设置工作目录
setwd(file.path(GTR_DEF, "代码调试"))

# 1. 创建病毒名称
viruses <- c(
  "Helarchaeota virus Nidhogg Meg22_1012",
  "Helarchaeota virus Nidhogg Meg22_1214", 
  "Helarchaeota virus Ratatoskr Meg22_1012",
  "Lokiarchaeota virus Fenrir Meg22_1012",
  "Lokiarchaeota virus Fenrir Meg22_1214",
  "Lokiarchaeota virus Sköll Meg22_1214"
)

# 2. 创建蛋白质特征名称
features <- c(
  # 病毒结构和形态发生 (12个)
  "Major capsid protein", "Tail fibre", "Baseplate J", "Terminase large subunit",
  "Baseplate wedge", "Minor head protein", "Tail completion protein S", 
  "Portal protein", "Tail sheath protein", "Aminopeptidase", 
  "HNH endonuclease", "I7/C57 endopeptidase-like",
  
  # DNA复制、重组和修复 (9个)
  "AAA ATPase repC-like", "DNA ligase with TFIIS-like Zn finger", 
  "3' flap repair endonuclease", "ERCC4 XPF DNA repair endonuclease", 
  "PCNA", "NucS mismatch repair endonuclease", "Archaeal-eukaryote primase", 
  "Primase-polymerase", "DNA polymerase family B",
  
  # RNA和转录 (7个)
  "RNA ligase", "Homeobox-like", "RNaseH1-like", "Ro60/TROVE2", 
  "Elf1-like glycosyltransferase", "Elp3/MiaB/NifB", 
  "Ribonucleoside diphosphate reductase",
  
  # 核苷酸代谢 (4个)
  "Deoxynucleo(tide/side) monophosphate kinase", "Cytidine deaminase", 
  "5'(3')-deoxyribonucleotidase", "Thioredoxin",
  
  # 宿主-病毒相互作用 (3个)
  "Von willebrand factor A", "Ubiquitin-activating enzyme", 
  "Homing endonuclease-intein"
)

# 3. 创建数据矩阵
create_data <- function() {
  data_list <- list()
  
  for (i in 1:length(viruses)) {
    for (j in 1:length(features)) {
      association <- "None"
      
      # 根据图像描述设置关联类型
      if (i <= 5) {  # 前5个病毒
        if (j <= 12) {  # 病毒结构和形态发生特征
          if (j %in% c(1, 3, 4, 7, 8)) {
            association <- "Caudovirales"
          } else if (j %in% c(11, 12)) {
            association <- "NCVOG"
          }
        } else if (j == 18) {  # NucS mismatch repair endonuclease
          association <- "Other"
        } else if (j == 21) {  # DNA polymerase family B
          association <- "NCVOG"
        } else if (i == 2 && j == 34) {  # 第二个病毒的Ubiquitin-activating enzyme
          association <- "Other"
        }
      }
      
      data_list[[length(data_list) + 1]] <- data.frame(
        virus = viruses[i],
        feature = features[j],
        association = association,
        feature_num = j,
        virus_num = i
      )
    }
  }
  
  return(do.call(rbind, data_list))
}

# 创建数据
plot_data <- create_data()

# 4. 设置因子顺序
plot_data$virus <- factor(plot_data$virus, levels = rev(viruses))
plot_data$feature <- factor(plot_data$feature, levels = features)
plot_data$association <- factor(plot_data$association, 
                               levels = c("Caudovirales", "NCVOG", "Other", "None"))

# 5. 创建分类背景数据
category_ranges <- list(
  c(1, 12),    # 病毒结构和形态发生
  c(13, 21),   # DNA复制、重组和修复
  c(22, 28),   # RNA和转录
  c(29, 32),   # 核苷酸代谢
  c(33, 35)    # 宿主-病毒相互作用
)

category_colors <- c("#E6E6FA", "#E0F6FF", "#E0F6FF", "#E0F6FF", "#FFE6E6")
category_names <- c("Virion structure and morphogenesis", 
                   "DNA replication, recombination and repair",
                   "RNA and transcription", 
                   "Nucleotide metabolism",
                   "Host-virus interaction")

# 6. 创建背景矩形数据
background_data <- data.frame()
for (i in 1:length(category_ranges)) {
  range <- category_ranges[[i]]
  background_data <- rbind(background_data, data.frame(
    xmin = range[1] - 0.5,
    xmax = range[2] + 0.5,
    ymin = 0.5,
    ymax = length(viruses) + 0.5,
    category = category_names[i],
    fill_color = category_colors[i]
  ))
}

# 7. 创建颜色映射
association_colors <- c(
  "Caudovirales" = "#8B4B8C",  # 紫色
  "NCVOG" = "#20B2AA",         # 青绿色
  "Other" = "#696969",         # 深灰色
  "None" = "white"             # 白色
)

# 8. 绘制图形
p <- ggplot() +
  
  # 添加分类背景
  geom_rect(data = background_data,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = background_data$fill_color,
            alpha = 0.3) +
  
  # 添加点
  geom_point(data = plot_data,
             aes(x = feature_num, y = virus_num, 
                 color = association, fill = association),
             size = 4, shape = 21, stroke = 0.5) +
  
  # 设置颜色
  scale_color_manual(values = association_colors, name = "Association Type") +
  scale_fill_manual(values = association_colors, name = "Association Type") +
  
  # 设置坐标轴
  scale_x_continuous(breaks = 1:length(features), labels = features) +
  scale_y_continuous(breaks = 1:length(viruses), labels = rev(viruses)) +
  
  # 设置主题
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 9),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 10)
  ) +
  
  # 添加标题
  labs(
    title = "Archaeal Virus Protein Feature Association",
    subtitle = "Presence and association type of proteins across six archaeal viruses"
  )

# 9. 添加分类标签
for (i in 1:length(category_ranges)) {
  range <- category_ranges[[i]]
  x_pos <- mean(range)
  p <- p + 
    annotate("text", x = x_pos, y = length(viruses) + 0.8, 
             label = category_names[i], angle = 45, hjust = 0, 
             size = 3, fontface = "bold")
}

# 10. 显示图形
print(p)

# 11. 保存图形
ggsave("virus_protein_dotplot.png", plot = p, width = 16, height = 8, dpi = 300)
ggsave("virus_protein_dotplot.pdf", plot = p, width = 16, height = 8)

# 12. 输出数据摘要
cat("数据摘要:\n")
cat("病毒数量:", length(viruses), "\n")
cat("蛋白质特征数量:", length(features), "\n")
cat("功能分类数量:", length(category_ranges), "\n")

# 显示关联类型统计
association_summary <- table(plot_data$association)
cat("\n关联类型统计:\n")
print(association_summary)

cat("\n图形已保存为 virus_protein_dotplot.png 和 virus_protein_dotplot.pdf\n")
