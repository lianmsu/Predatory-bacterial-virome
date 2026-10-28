# =============================================================================
# Predatory-bacterial-virome -- phage annotation / Gene presence/absence dot-plot
# =============================================================================
# Dot-plot of the annotated genes across a set of viruses (filled circle =
# orthologue present); the script installs its packages and builds the plot
# from a prepared table.
#
# Manuscript      : Fig. S6
# Main outputs    : p_phage_gene_anno2.pdf
# Original script : 1mobilome_predatory_bacteria/7phage_gene_annotation/phage_gene_anno_plot/create_virus_dotplot.R
# -----------------------------------------------------------------------------
# Input data and output files are resolved from GTR_DATA (current analysis root),
# GTR_WORK (earlier working copy) and GTR_DEF (shared defensome data). All three
# default to the current working directory; set them instead of editing files:
#
#     Sys.setenv(GTR_DATA = "/path/to/06_predatory_bacteria_virome")
#
# =============================================================================
GTR_DATA <- normalizePath(Sys.getenv("GTR_DATA", unset = "."), mustWork = FALSE)

# 古菌病毒蛋白质特征点图
# 完整的R脚本，包含包安装和错误处理

# 检查并安装必要的包
required_packages <- c("ggplot2", "dplyr", "tidyr")

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# 设置工作目录
current_dir <- getwd()
cat("当前工作目录:", current_dir, "\n")

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
cat("创建数据矩阵...\n")
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
        virus_num = i,
        stringsAsFactors = FALSE
      )
    }
  }
  
  return(do.call(rbind, data_list))
}

# 创建数据
plot_data <- create_data()
cat("数据矩阵创建完成，共", nrow(plot_data), "行数据\n")

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

# 创建背景矩形数据
background_data <- data.frame()
for (i in 1:length(category_ranges)) {
  range <- category_ranges[[i]]
  background_data <- rbind(background_data, data.frame(
    xmin = range[1] - 0.5,
    xmax = range[2] + 0.5,
    ymin = 0.5,
    ymax = length(viruses) + 0.5,
    category = category_names[i],
    fill_color = category_colors[i],
    stringsAsFactors = FALSE
  ))
}

# 6. 创建颜色映射
association_colors <- c(
  "Caudovirales" = "#8B4B8C",  # 紫色
  "NCVOG" = "#20B2AA",         # 青绿色
  "Other" = "#696969",         # 深灰色
  "None" = "white"             # 白色
)

# 7. 绘制图形
cat("开始绘制图形...\n")
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

# 8. 添加分类标签
for (i in 1:length(category_ranges)) {
  range <- category_ranges[[i]]
  x_pos <- mean(range)
  p <- p + 
    annotate("text", x = x_pos, y = length(viruses) + 0.8, 
             label = category_names[i], angle = 45, hjust = 0, 
             size = 3, fontface = "bold")
}

################

my_phage <- readxl::read_xlsx("all_17_phage_gff_withtype2_simplifiedname.xlsx")
phage_14 <- readxl::read_xlsx("Candidate_14_new_phages.xlsx")

my_phage <- my_phage[my_phage$contig %in% phage_14$clean_seq,]
my_phage <- my_phage %>%
  left_join(phage_14 %>% select(clean_seq, Name), 
            by = c("contig" = "clean_seq")) %>%
  rename(virus_name = Name)

my_phage <- my_phage[,c(11,2,7)]

head(my_phage)

# 可自定义过滤条件的病毒蛋白质点图绘制脚本
library(ggplot2)
library(dplyr)


# 检查my_phage数据框是否存在
if (!exists("my_phage")) {
  cat("错误：my_phage数据框不存在！\n")
  stop("数据框不存在")
}

# 过滤参数设置
MIN_VIRUS_COUNT <- 2  # 蛋白质至少需要在多少个病毒中出现
# 可以修改这个值：
# MIN_VIRUS_COUNT <- 1  # 不过滤，保留所有蛋白质
# MIN_VIRUS_COUNT <- 3  # 只保留在3个或更多病毒中出现的蛋白质

cat("过滤参数: 蛋白质至少需要在", MIN_VIRUS_COUNT, "个病毒中出现\n")

# 数据预处理
cat("开始数据预处理...\n")

# 1. 获取唯一值
viruses <- unique(my_phage$virus_name)
proteins <- unique(my_phage$name[!is.na(my_phage$name)])
categories <- unique(my_phage$type)

cat("原始数据:\n")
cat("病毒数量:", length(viruses), "\n")
cat("蛋白质数量:", length(proteins), "\n")
cat("功能分类:", paste(categories, collapse = ", "), "\n")

# 2. 过滤蛋白质
cat("\n过滤蛋白质...\n")

# 统计每个蛋白质出现的病毒数量
protein_counts <- my_phage %>%
  group_by(name) %>%
  summarise(virus_count = n(), .groups = 'drop')

# 显示蛋白质出现频率分布
cat("蛋白质出现频率分布:\n")
freq_table <- table(protein_counts$virus_count)
print(freq_table)

# 应用过滤条件
proteins_filtered <- protein_counts %>%
  filter(virus_count >= MIN_VIRUS_COUNT) %>%
  pull(name)

cat("过滤结果:\n")
cat("原始蛋白质数量:", length(proteins), "\n")
cat("过滤后蛋白质数量:", length(proteins_filtered), "\n")
cat("删除的蛋白质数量:", length(proteins) - length(proteins_filtered), "\n")

# 显示被删除的蛋白质
deleted_proteins <- setdiff(proteins, proteins_filtered)
if (length(deleted_proteins) > 0) {
  cat("被删除的蛋白质:\n")
  print(deleted_proteins)
}

# 如果过滤后没有蛋白质，停止执行
if (length(proteins_filtered) == 0) {
  cat("警告：过滤后没有蛋白质剩余！请调整过滤条件。\n")
  stop("没有蛋白质可用于绘图")
}

# 3. 创建绘图数据
cat("\n创建绘图数据...\n")
plot_data <- expand.grid(
  virus = viruses,
  protein = proteins_filtered,
  stringsAsFactors = FALSE
)

# 合并原始数据获取type信息
plot_data <- merge(plot_data, my_phage, 
                   by.x = c("virus", "protein"), 
                   by.y = c("virus_name", "name"), 
                   all.x = TRUE)

# 标记蛋白质是否存在
plot_data$present <- !is.na(plot_data$type)

# 4. 设置因子顺序
plot_data$virus <- factor(plot_data$virus, levels = viruses)
plot_data$virus <- factor(plot_data$virus, levels = c("Xuanmingvirus","Dizangvirus",
                                                      "Shennongvirus","Taotievirus",
                                                      "Aobingvirus","Qiongqivirus",
                                                      "Hairuovirus","Qinyuanvirus",
                                                      "Zhulongvirus","Yinglongvirus",
                                                      "Xiangliuvirus","Rushouvirus",
                                                      "Bingyivirus","Xirangvirus"))
#plot_data$protein <- factor(plot_data$protein, levels = proteins_filtered)

# 根据type分类和首字母对蛋白质进行排序
protein_order <- plot_data %>%
  select(protein, type) %>%
  distinct() %>%
  arrange(type, protein) %>%
  pull(protein)
protein_order <- protein_order[1:54]

plot_data$protein <- factor(plot_data$protein, levels = protein_order)

# 5. 创建颜色映射
# n_colors <- length(categories)
# if (n_colors <= 8) {
#   colors <- RColorBrewer::brewer.pal(max(3, n_colors), "Set3")[1:n_colors]
# } else {
#   colors <- rainbow(n_colors)
# }

colors <- c("#8020f3","#01cabf","#1f78b4","#cab2d6","#d9d9d9","#ff7f00","#e5049c")
names(colors) <- categories
# 将NA值转换为字符串"NA"
plot_data$type[is.na(plot_data$type)] <- "NA"
# 为NA值添加白色
colors <- c(colors, "white")
names(colors)[length(colors)] <- "NA"

# 6. 绘制图形
cat("开始绘制图形...\n")

p <- ggplot(plot_data, aes(x = protein, y = virus)) +
  
  # 添加点
  geom_point(aes( fill = type), color = "black",
             size = 4, shape = 21, stroke = 0.2) +
  
  # 设置颜色
#  scale_color_manual(values = colors, name = "Protein Type") +
  scale_fill_manual(values = colors, name = "Protein Type",na.value = "white") +
  
  # 设置主题
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 10)
  ) 
  
  # # 添加标题
  # labs(
  #   #title = paste("Virus protein feature association (filtered, min =", MIN_VIRUS_COUNT, ")"),
  #   subtitle = paste("Protein feature across", length(viruses), "viruses", "(filtered, min = 2)")
  # )

# 7. 显示图形
print(p)

ggsave("p_phage_gene_anno2.pdf",plot = p, units = 'in', width = 9, height = 4.5)

# 9. 数据摘要
cat("\n=== 过滤后数据摘要 ===\n")
cat("过滤条件: 蛋白质至少需要在", MIN_VIRUS_COUNT, "个病毒中出现\n")
cat("病毒数量:", length(viruses), "\n")
cat("过滤后蛋白质数量:", length(proteins_filtered), "\n")

# 各功能分类的统计（基于过滤后的数据）
filtered_data <- my_phage[my_phage$name %in% proteins_filtered, ]
type_summary <- table(filtered_data$type)
cat("\n过滤后各功能分类的蛋白质数量:\n")
print(type_summary)

# 各病毒的蛋白质数量（基于过滤后的数据）
virus_summary <- table(filtered_data$contig)
cat("\n过滤后各病毒的蛋白质数量:\n")
print(virus_summary)

# 创建过滤后的交叉表
summary_table <- filtered_data %>%
  group_by(contig, type) %>%
  summarise(count = n(), .groups = 'drop') %>%
  pivot_wider(names_from = type, values_from = count, values_fill = 0)

cat("\n过滤后病毒-功能分类交叉表:\n")
print(summary_table)

cat("\n=== 脚本执行完成 ===\n")
cat("过滤后的图形已保存为", filename_base, ".png 和", filename_base, ".pdf\n")






