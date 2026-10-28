# =============================================================================
# Predatory-bacterial-virome -- phage annotation / Gene dot-plot, full variant
# =============================================================================
# Full variant of the gene dot-plot with the complete styling; reads and
# writes inside an external prototyping folder (GTR_DEF).
#
# Manuscript      : Fig. S6 (variant)
# Main outputs    : virus_protein_dotplot.pdf
# Original script : 1mobilome_predatory_bacteria/7phage_gene_annotation/phage_gene_anno_plot/virus_protein_dotplot.R
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

# 古菌病毒蛋白质特征点图
# 基于图像描述创建的数据可视化

library(ggplot2)
library(dplyr)
library(tidyr)

# 设置工作目录和输出
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

# 2. 创建蛋白质特征名称（按功能分类）
# 病毒结构和形态发生 (12个特征)
virion_features <- c(
  "Major capsid protein",
  "Tail fibre", 
  "Baseplate J",
  "Terminase large subunit",
  "Baseplate wedge",
  "Minor head protein",
  "Tail completion protein S",
  "Portal protein",
  "Tail sheath protein",
  "Aminopeptidase",
  "HNH endonuclease",
  "I7/C57 endopeptidase-like"
)

# DNA复制、重组和修复 (9个特征)
dna_features <- c(
  "AAA ATPase repC-like",
  "DNA ligase with TFIIS-like Zn finger",
  "3' flap repair endonuclease", 
  "ERCC4 XPF DNA repair endonuclease",
  "PCNA",
  "NucS mismatch repair endonuclease",
  "Archaeal-eukaryote primase",
  "Primase-polymerase",
  "DNA polymerase family B"
)

# RNA和转录 (7个特征)
rna_features <- c(
  "RNA ligase",
  "Homeobox-like",
  "RNaseH1-like",
  "Ro60/TROVE2",
  "Elf1-like glycosyltransferase",
  "Elp3/MiaB/NifB",
  "Ribonucleoside diphosphate reductase"
)

# 核苷酸代谢 (4个特征)
nucleotide_features <- c(
  "Deoxynucleo(tide/side) monophosphate kinase",
  "Cytidine deaminase",
  "5'(3')-deoxyribonucleotidase",
  "Thioredoxin"
)

# 宿主-病毒相互作用 (3个特征)
host_virus_features <- c(
  "Von willebrand factor A",
  "Ubiquitin-activating enzyme", 
  "Homing endonuclease-intein"
)

# 合并所有特征
all_features <- c(virion_features, dna_features, rna_features, nucleotide_features, host_virus_features)

# 3. 创建功能分类信息
feature_categories <- data.frame(
  feature = all_features,
  category = c(
    rep("Virion structure and morphogenesis", length(virion_features)),
    rep("DNA replication, recombination and repair", length(dna_features)),
    rep("RNA and transcription", length(rna_features)),
    rep("Nucleotide metabolism", length(nucleotide_features)),
    rep("Host-virus interaction", length(host_virus_features))
  ),
  category_order = c(
    rep(1, length(virion_features)),
    rep(2, length(dna_features)),
    rep(3, length(rna_features)),
    rep(4, length(nucleotide_features)),
    rep(5, length(host_virus_features))
  )
)

# 4. 根据图像描述创建数据矩阵
# 创建关联类型数据
create_association_data <- function() {
  # 初始化数据框
  data_list <- list()
  
  # 为每个病毒创建数据
  for (i in 1:length(viruses)) {
    virus_data <- data.frame(
      virus = viruses[i],
      feature = all_features,
      association = "None",  # 默认为无关联
      stringsAsFactors = FALSE
    )
    
    # 根据图像描述设置关联类型
    if (i <= 5) {  # 前5个病毒
      # 病毒结构和形态发生特征 (前12个)
      virus_data$association[1:12] <- c("Caudovirales", "None", "Caudovirales", "Caudovirales", 
                                       "None", "None", "Caudovirales", "Caudovirales", 
                                       "None", "None", "NCVOG", "NCVOG")
      
      # DNA复制、重组和修复特征 (13-21个)
      virus_data$association[13:21] <- c("None", "None", "None", "None", "None", 
                                        "Other", "None", "None", "NCVOG")
      
      # 其他特征保持为None
      virus_data$association[22:35] <- "None"
      
      # 特殊处理：第二个病毒的宿主-病毒相互作用特征
      if (i == 2) {
        virus_data$association[34] <- "Other"  # Ubiquitin-activating enzyme
      }
    } else {  # 第6个病毒 (Sköll)
      # 所有特征都无关联
      virus_data$association <- "None"
    }
    
    data_list[[i]] <- virus_data
  }
  
  return(do.call(rbind, data_list))
}

# 创建数据
plot_data <- create_association_data()

# 合并功能分类信息
plot_data <- merge(plot_data, feature_categories, by = "feature", all.x = TRUE)

# 5. 设置因子顺序
plot_data$virus <- factor(plot_data$virus, levels = rev(viruses))
plot_data$feature <- factor(plot_data$feature, levels = all_features)
plot_data$association <- factor(plot_data$association, 
                               levels = c("Caudovirales", "NCVOG", "Other", "None"))

# 6. 创建分类背景数据
category_background <- feature_categories %>%
  group_by(category, category_order) %>%
  summarise(
    start = min(as.numeric(feature)),
    end = max(as.numeric(feature)),
    .groups = 'drop'
  ) %>%
  mutate(
    y_min = 0.5,
    y_max = length(viruses) + 0.5
  )

# 7. 创建颜色映射
association_colors <- c(
  "Caudovirales" = "#8B4B8C",  # 紫色
  "NCVOG" = "#20B2AA",         # 青绿色
  "Other" = "#696969",         # 深灰色
  "None" = "white"             # 白色
)

# 8. 绘制图形
p <- ggplot(plot_data, aes(x = feature, y = virus)) +
  
  # 添加分类背景
  geom_rect(data = category_background, 
            aes(xmin = start - 0.5, xmax = end + 0.5, ymin = y_min, ymax = y_max),
            fill = c("#E6E6FA", "#E0F6FF", "#E0F6FF", "#E0F6FF", "#FFE6E6"),
            alpha = 0.3, inherit.aes = FALSE) +
  
  # 添加点
  geom_point(aes(color = association, fill = association), 
             size = 4, shape = 21, stroke = 0.5) +
  
  # 设置颜色
  scale_color_manual(values = association_colors, name = "Association Type") +
  scale_fill_manual(values = association_colors, name = "Association Type") +
  
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
  ) +
  
  # 反转y轴使第一个病毒在顶部
  scale_y_discrete(limits = rev(levels(plot_data$virus)))

# 9. 添加分类标签
# 计算每个分类的中间位置
category_labels <- feature_categories %>%
  group_by(category, category_order) %>%
  summarise(
    x_pos = mean(as.numeric(feature)),
    .groups = 'drop'
  ) %>%
  mutate(
    y_pos = length(viruses) + 0.8
  )

# 添加分类标签到图形
p <- p + 
  geom_text(data = category_labels,
            aes(x = x_pos, y = y_pos, label = category),
            angle = 45, hjust = 0, size = 3, fontface = "bold")

# 10. 显示图形
print(p)

# 11. 保存图形
ggsave("virus_protein_dotplot.png", plot = p, width = 16, height = 8, dpi = 300)
ggsave("virus_protein_dotplot.pdf", plot = p, width = 16, height = 8)

# 12. 输出数据摘要
cat("数据摘要:\n")
cat("病毒数量:", length(viruses), "\n")
cat("蛋白质特征数量:", length(all_features), "\n")
cat("功能分类数量:", length(unique(feature_categories$category)), "\n")

# 显示关联类型统计
association_summary <- table(plot_data$association)
cat("\n关联类型统计:\n")
print(association_summary)
