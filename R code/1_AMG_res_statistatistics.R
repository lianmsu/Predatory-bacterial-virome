setwd("D:/02_pku/01_research/06_predatory_bacteria_virome/06_allphage_AMGs")

library(dplyr)
library(tidyr)

phage_AMG_anno <- 
  read.delim2("output_annotated_all_AMG.ORF2GENE_with_meta.tsv",header = F)
colnames(phage_AMG_anno) <- 
  c("gene_id","gene_name","level1","level2","level3")
# 定义要排除的基因名
exclude_genes <- c("dcm", "queC", "queD", "queE", "queF")

# 过滤数据，保留不在排除列表中的行
phage_AMG_anno <- phage_AMG_anno[!phage_AMG_anno$gene_name %in% exclude_genes, ]

table(phage_AMG_anno$gene_name)

##
library(ggplot2)
library(dplyr)
library(viridis)

# 统计每个 gene_name 的数量
gene_counts <- phage_AMG_anno %>%
  group_by(level3, gene_name) %>%
  summarise(count = n(), .groups = "drop")

# 按 level3 分组重新排序 gene_name
gene_counts <- gene_counts %>%
  arrange(level3, gene_name) %>%
  mutate(gene_name = factor(gene_name, levels = gene_name))

# 自动生成颜色
n_colors <- length(unique(gene_counts$level3))
colors <- viridis(n_colors, option = "turbo")

# 绘图，添加柱子数量标签
ggplot(gene_counts, aes(x = gene_name, y = count, fill = level3)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = count), vjust = -0.3, size = 3) +  # 数字标签在柱子上方
  scale_fill_manual(values = colors) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank()
  ) +
  labs(y = "Count", fill = "Level3")

#基因数量####

library(ggplot2)
library(dplyr)
library(viridis)

# 统计每个 gene_name 的数量
gene_counts <- phage_AMG_anno %>%
  group_by(level2, gene_name) %>%
  summarise(count = n(), .groups = "drop")

# 按 level2 分组重新排序 gene_name
gene_counts <- gene_counts %>%
  arrange(level2, gene_name) %>%
  mutate(gene_name = factor(gene_name, levels = gene_name))

# 自动生成颜色，按 level2
n_colors <- length(unique(gene_counts$level2))
colors <- viridis(n_colors, option = "turbo")

# 绘图，横向柱状图，添加数量标签
p1_AMG_number <- 
ggplot(gene_counts, aes(x = gene_name, y = count, fill = level2)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = count), hjust = 0.5, vjust = -1,size = 3) +  # 数字标签在柱子末端
  scale_fill_manual(values = colors) +
#  coord_flip() +  # 横向翻转
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 10),
    axis.title.y = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  ) +

  labs(y = "Count", fill = "Level2")
p1_AMG_number

p1_AMG_number_legend <- 
  ggplot(gene_counts, aes(x = gene_name, y = count, fill = level2)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = count), hjust = 0.5, vjust = -1,size = 3) +  # 数字标签在柱子末端
  scale_fill_manual(values = colors) +
  #  coord_flip() +  # 横向翻转
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 10),
    axis.title.y = element_blank(),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  ) +
  labs(y = "Count", fill = "Level2")

p1_AMG_number_legend

ggsave(plot = p1_AMG_number_legend,
       '1_p1_AMG_number_legend.pdf',
       width = 7.2,height = 4,
       units = 'in')

#lifestyle####
bacphilip_pre <- read.delim2("Myxococcota_all_phage_info_filterd.fa.bacphlip")

bacphilip_pre <- bacphilip_pre %>%
  mutate(lifestyle = case_when(
    grepl("prophage", Phage) ~ "Temperate",
    Virulent >= 0.9        ~ "Virulent",
    Temperate >= 0.9       ~ "Temperate",
    TRUE                    ~ NA_character_
  ))
bacphilip_pre$lifestyle[is.na(bacphilip_pre$lifestyle)] <- "Unknown"

all_phage_lifestyle <- table(bacphilip_pre$lifestyle)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

# 准备数据
df_lifestyle <- as.data.frame(all_phage_lifestyle)
colnames(df_lifestyle) <- c("Lifestyle", "Count")

# 计算百分比和标签位置
df_lifestyle <- df_lifestyle %>%
  mutate(
    Percentage = Count / sum(Count) * 100,
    Label = paste0(Count, " (", round(Percentage, 1), "%)"),
    ymax = cumsum(Percentage),
    ymin = c(0, head(ymax, n = -1))
  )

# 创建环形图
p_lifestyle_pie <- ggplot(df_lifestyle, aes(ymax = ymax, ymin = ymin, 
                                            xmax = 4, xmin = 3, 
                                            fill = Lifestyle)) +
  geom_rect() +
  coord_polar(theta = "y") +
  xlim(c(2, 4)) +  # 控制环形的宽度
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "right"
  ) +
  scale_fill_brewer(palette = "Set2") +  # 使用漂亮的颜色
  geom_text(aes(x = 3.5, y = (ymin + ymax) / 2, 
                label = paste0(Lifestyle, "\n", Label)),
            size = 3.5, color = "black")

print(p_lifestyle_pie)

ggsave(plot = p_lifestyle_pie,'1_p_lifestyle_pie.pdf',width = 4.4,height = 4.2,units = 'in')


# library(Cairo)
# library(cowplot)
# 
# g1_lifestyle_AMG <- plot_grid(p_lifestyle_pie,p1_AMG_number_legend,
#                 labels = LETTERS[1:2], ncol = 2,
#                 rel_widths = c(1, 2)
#                 # label_x = 0.09,            # 标签x坐标（1为最右侧）
#                 #  label_y = 0.9             # 标签y坐标（1为最上方）
# )
# g1_lifestyle_AMG

ggsave(plot = g1_lifestyle_AMG,'1_MGEnum_antiDS2.pdf',width = 8,height = 2.5,units = 'in')







save.image("1_AMG_res_statistatistics.RData")
load("1_AMG_res_statistatistics.RData")




