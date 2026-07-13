
setwd("D:/02_pku/01_research/06_predatory_bacteria_virome/10_myxo_function_anno")

library(stringr)
library(dplyr)
library(tidyr)
#
myxo_genome_metadata <- read.delim2("myxococcota_mag_metadata.tsv")


#DS####
## data prepare####
myxo_DF_gene <- read.delim2("merged_myxo_mag_defense_finder_genes.tsv")
DF_gene_gff <- read.delim2("DF_genes.gff", header = F)
colnames(DF_gene_gff) <-
  c("contig","source","CDS","start","end","info","strand","info2","detail")
DF_gene_gff$gene_ID <- str_extract(DF_gene_gff$detail, "(?<=^ID=)[^;]+")
DF_gene_info <- left_join(myxo_DF_gene, DF_gene_gff, 
                          by = c("hit_id" = "gene_ID"))
DF_gene_info <- 
  select(DF_gene_info, replicon, hit_id, start, end, 
         strand, gene_name, type, subtype, activity)
DF_system_info <- 
  read.delim2("merged_myxo_mag_defense_finder_systems.tsv")

#用DF_defense的sys_beg去匹配DG数据框的hit_id列，匹配到的行的V4列的值作为DF_defense对应行的begin
#用DF_defense的sys_end去匹配DG数据框的hit_id列，匹配到的行的V5列的值作为DF_defense对应行的end
# 首先为DF_defense创建begin和end列（如果不存在）
if(!"DSbegin" %in% names(DF_system_info)) DF_system_info$DSbegin <- NA
if(!"DSend" %in% names(DF_system_info)) DF_system_info$DSend <- NA
# 创建从hit_id到V4/V5的映射
v4_mapping <- setNames(DF_gene_info$start, DF_gene_info$hit_id)
v5_mapping <- setNames(DF_gene_info$end, DF_gene_info$hit_id)
# 用sys_beg匹配并填充begin列
DF_system_info$DSbegin <- ifelse(
  DF_system_info$sys_beg %in% names(v4_mapping),
  v4_mapping[as.character(DF_system_info$sys_beg)],
  DF_system_info$DSbegin
)
# 用sys_end匹配并填充end列
DF_system_info$DSend <- ifelse(
  DF_system_info$sys_end %in% names(v5_mapping),
  v5_mapping[as.character(DF_system_info$sys_end)],
  DF_system_info$DSend
)

DF_system_info <- DF_system_info %>%
  mutate(
    contig = str_replace(sys_beg, "_\\d+$", "")  
  )
DF_system_info <- DF_system_info %>%
  mutate(
    genome = str_replace(contig, "-\\d+$", "")  
  )
DF_system_info <- 
  DF_system_info %>%
  left_join(myxo_genome_metadata, by = c("genome" = "Genome"))

##统计及绘图####
DS <- DF_system_info[DF_system_info$activity == "Defense",]

###各类DS总数####
library(ggplot2)

# 统计每个 type 的出现次数（直接按行计）
DS_type_counts <- DS %>%
  count(type, sort = TRUE) %>%
  top_n(20, n)              # 保留数量最多的前 50 个 type


plot_DS_num <-
ggplot(DS_type_counts, aes(x = reorder(type, n), y = n)) +
  geom_col(fill = "#4A90D9", width = 0.7) +
  coord_flip() +                       # 横向条形图，便于展示
  # labs(x = "Type", y = "Count", title = "Top 50 Most Abundant Defense System Types") +
  # theme_minimal(base_size = 12) +
  # theme(panel.grid.major.y = element_blank(),
  #       axis.text.y = element_text(size = 9))+
  theme(
    panel.background = element_blank(),  # 设置背景为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    panel.border = element_rect(color = "grey20", fill = NA, linewidth = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    axis.title = element_blank(),
    #axis.title = element_text( size = 12, color = "black"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 10, color = "black"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5),  # 横轴标签旋转
    axis.ticks = element_line(color = "grey20", linewidth = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    legend.position = 'none',
    #legend.position = c(0.2,0.75),  # 设置图例位置在右侧
    legend.title = element_blank(), #设置图例标题字体为粗体，大小为12，颜色为深灰色
    legend.text = element_blank(), #设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", linewidth = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )
plot_DS_num

###DS流行率####
library(dplyr)
# 统计每个 class 的 MAG 数量，取前 10
top10_classes <- myxo_genome_metadata %>%
  count(class, name = "n_mags") %>%
  slice_max(n_mags, n = 10) %>%
  pull(class)
# 2.1 计算每个 class 的 MAG 总数（依据 myxo_genome_metadata）
class_mag_counts <- myxo_genome_metadata %>%
  filter(class %in% top10_classes) %>%
  count(class, name = "total_mags")

# 2.2 计算每个 class-type 组合的独特 MAG 数
class_DStype_mag <- DS %>%
  filter(class %in% top10_classes) %>%
  distinct(class, genome, type) %>%    # 同一 MAG 同一 type 只计一次
  count(class, type, name = "n_mags")

# 2.3 合并并计算流行率
DS_prevalence <- class_DStype_mag %>%
  left_join(class_mag_counts, by = "class") %>%
  mutate(prevalence = 100 * n_mags / total_mags)

# 2.4 选择要展示的 Type（例如：至少在一个 class 中流行率 > 5% 的 type）
DS_type_to_show <- DS_prevalence %>%
  group_by(type) %>%
  summarise(max_prev = max(prevalence)) %>%
  filter(max_prev > 5) %>%
  pull(type)

DS_heatmap_data <- DS_prevalence %>%
  filter(type %in% DS_type_to_show)

# 2.5 绘制热图
plot_DS_prevalence_heatmap <-
ggplot(DS_heatmap_data, aes(x = type, y = class, fill = prevalence)) +
  geom_tile(color = "white", lwd = 0.5) +
  scale_fill_gradient(low = "#f7fbff", high = "#08519c", name = "Prevalence (%)") +
  labs(x = "Defense System Type", y = "Class (Myxococcota)",
       title = "Prevalence of Defense Systems Across Top 10 Classes") +
  # theme_minimal(base_size = 12) +
  # theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
  #       axis.text.y = element_text(size = 10),
  #       legend.position = "right")+
  theme(
    panel.background = element_blank(),  # 设置背景为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    panel.border = element_rect(color = "grey20", fill = NA, linewidth = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    axis.title = element_blank(),
    #axis.title = element_text( size = 12, color = "black"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 10, color = "black"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),  # 横轴标签旋转
    axis.ticks = element_line(color = "grey20", linewidth = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    #legend.position = 'none',
    legend.position = c(0.2,0.75),  # 设置图例位置在右侧
    legend.title = element_blank(), #设置图例标题字体为粗体，大小为12，颜色为深灰色
    #legend.text = element_blank(), #设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", linewidth = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )
plot_DS_prevalence_heatmap

ggsave(plot = plot_DS_prevalence_heatmap, 
       "1_plot_DS_prevalence_heatmap_legend.pdf",
       width = 8, height = 4, units = "in")

#AMP####
##data prepare####
AMP <- read.delim2("macrel.out.prediction")
AMP_gene_gff <- read.delim2("AMP_genes.gff",header = F)
colnames(AMP_gene_gff) <-
  c("contig","source","CDS","start","end","info","strand","info2","detail")
AMP_gene_gff$gene_ID <- str_extract(AMP_gene_gff$detail, "(?<=^ID=)[^;]+")
AMP <- left_join(AMP, AMP_gene_gff, 
                 by = c("Access" = "gene_ID"))
AMP <- AMP %>%
  mutate(
    genome = str_replace(contig, "-\\d+$", "")  
  )
AMP <- 
  AMP %>%
  left_join(myxo_genome_metadata, by = c("genome" = "Genome"))
head(AMP)

##统计及绘图####

###各类AMP总数####
library(ggplot2)

# 统计每个 type 的出现次数（直接按行计）
AMP_type_counts <- AMP %>%
  count(AMP_family, sort = TRUE) %>%
  top_n(50, n)              # 保留数量最多的前 50 个 type

library(ggplot2)
plot_AMP_num <-
  ggplot(AMP_type_counts, aes(x = reorder(AMP_family, n), y = n)) +
  geom_col(fill = "#9b5de5", width = 0.7) +
  coord_flip() +                       # 横向条形图，便于展示
  theme(
    panel.background = element_blank(),  # 设置背景为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    panel.border = element_rect(color = "grey20", fill = NA, linewidth = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    axis.title = element_blank(),
    #axis.title = element_text( size = 12, color = "black"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 10, color = "black"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5),  # 横轴标签旋转
    axis.ticks = element_line(color = "grey20", linewidth = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    legend.position = 'none',
    #legend.position = c(0.2,0.75),  # 设置图例位置在右侧
    legend.title = element_blank(), #设置图例标题字体为粗体，大小为12，颜色为深灰色
    legend.text = element_blank(), #设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", linewidth = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )
plot_AMP_num

###DS流行率####
library(dplyr)
# 统计每个 class 的 MAG 数量，取前 10
top10_classes <- myxo_genome_metadata %>%
  count(class, name = "n_mags") %>%
  slice_max(n_mags, n = 10) %>%
  pull(class)
# 2.1 计算每个 class 的 MAG 总数（依据 myxo_genome_metadata）
class_mag_counts <- myxo_genome_metadata %>%
  filter(class %in% top10_classes) %>%
  count(class, name = "total_mags")

# 2.2 计算每个 class-type 组合的独特 MAG 数
class_AMPtype_mag <- AMP %>%
  filter(class %in% top10_classes) %>%
  distinct(class, genome, AMP_family) %>%    # 同一 MAG 同一 type 只计一次
  count(class, AMP_family, name = "n_mags")

# 2.3 合并并计算流行率
AMP_prevalence <- class_AMPtype_mag %>%
  left_join(class_mag_counts, by = "class") %>%
  mutate(prevalence = 100 * n_mags / total_mags)

# 2.4 选择要展示的 Type（例如：至少在一个 class 中流行率 > 5% 的 type）
AMP_type_to_show <- AMP_prevalence %>%
  group_by(AMP_family) %>%
  summarise(max_prev = max(prevalence)) %>%
  filter(max_prev > 5) %>%
  pull(AMP_family)

AMP_heatmap_data <- AMP_prevalence %>%
  filter(AMP_family %in% AMP_type_to_show)

# 2.5 绘制热图
plot_AMP_prevalence_heatmap <-
  ggplot(AMP_heatmap_data, aes(x = AMP_family, y = class, fill = prevalence)) +
  geom_tile(color = "white", lwd = 0.5) +
  scale_fill_gradient(low = "#f7fbff", high = "#9b5de5", name = "Prevalence (%)") +
  labs(x = "Defense System Type", y = "Class (Myxococcota)",
       title = "Prevalence of Defense Systems Across Top 10 Classes") +
  # theme_minimal(base_size = 12) +
  # theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
  #       axis.text.y = element_text(size = 10),
  #       legend.position = "right")+
  theme(
    panel.background = element_blank(),  # 设置背景为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    panel.border = element_rect(color = "grey20", fill = NA, linewidth = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    axis.title = element_blank(),
    #axis.title = element_text( size = 12, color = "black"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 10, color = "black"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),  # 横轴标签旋转
    axis.ticks = element_line(color = "grey20", linewidth = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    #legend.position = 'none',
    legend.position = c(0.2,0.75),  # 设置图例位置在右侧
    legend.title = element_blank(), #设置图例标题字体为粗体，大小为12，颜色为深灰色
    #legend.text = element_blank(), #设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", linewidth = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )
plot_AMP_prevalence_heatmap

ggsave(plot = plot_AMP_prevalence_heatmap, 
       "1_plot_AMP_prevalence_heatmap.pdf",
       width = 8, height = 4, units = "in")


#TXSS####
##data prepare####

txss <- read.delim2("merged_best_solutions.tsv")
txss_gene_gff <- read.delim2("TXSS_genes.gff",header = F)
colnames(txss_gene_gff) <-
  c("contig","source","CDS","start","end","info","strand","info2","detail")
txss_gene_gff$gene_ID <- str_extract(txss_gene_gff$detail, "(?<=^ID=)[^;]+")
txss <- left_join(txss, txss_gene_gff, 
                 by = c("hit_id" = "gene_ID"))
txss <- txss %>%
  mutate(
    genome = str_replace(contig, "-\\d+$", "")  
  )
txss$type <- 
  sub(".*/", "", txss$model_fqn)

txss_systems <- txss %>%
  group_by(sys_id) %>%
  arrange(hit_pos, .by_group = TRUE) %>%   # 按基因顺序排序
  summarise(
    type = first(type),                    # 系统类型（每组相同，取第一个）
    subtype = "",                          # TXSS 结果通常无子类型，置空
    activity = "Secretion",                # 活性标记，这里统一为分泌
    sys_beg = first(hit_id),               # 系统起始蛋白ID
    sys_end = last(hit_id),                # 系统末尾蛋白ID
    protein_in_syst = paste(hit_id, collapse = ","),  # 所有蛋白ID
    genes_count = n(),                     # 蛋白数量
    name_of_profiles_in_sys = paste(gene_name, collapse = ","), # 所有profile名称
    .groups = "drop"
  ) %>%
  select(sys_id, type, subtype, activity, sys_beg, sys_end, 
         protein_in_syst, genes_count, name_of_profiles_in_sys)

# 如果 txss_systems 还没有坐标列，先创建
if(!"SSbegin" %in% names(txss_systems)) txss_systems$SSbegin <- NA
if(!"SSend" %in% names(txss_systems)) txss_systems$SSend <- NA

# 创建 hit_id -> start / end 的映射表
start_map <- setNames(txss$start, txss$hit_id)
end_map   <- setNames(txss$end,   txss$hit_id)
# 填充 sys_beg 对应的起始坐标
txss_systems$SSbegin <- ifelse(
  txss_systems$sys_beg %in% names(start_map),
  start_map[as.character(txss_systems$sys_beg)],
  txss_systems$SSbegin
)
# 填充 sys_end 对应的终止坐标
txss_systems$SSend <- ifelse(
  txss_systems$sys_end %in% names(end_map),
  end_map[as.character(txss_systems$sys_end)],
  txss_systems$SSend
)


txss_systems <- txss_systems %>%
  mutate(
    contig = str_replace(sys_beg, "_\\d+$", "")  
  )
txss_systems <- txss_systems %>%
  mutate(
    genome = str_replace(contig, "-\\d+$", "")  
  )

txss_systems <- 
  txss_systems %>%
  left_join(myxo_genome_metadata, by = c("genome" = "Genome"))
head(txss_systems)



##统计及绘图####

###各类txss总数####
library(ggplot2)

# 统计每个 type 的出现次数（直接按行计）
txss_type_counts <- txss_systems %>%
  count(type, sort = TRUE) %>%
  top_n(50, n)              # 保留数量最多的前 50 个 type

library(ggplot2)
plot_txss_num <-
  ggplot(txss_type_counts, aes(x = reorder(type, n), y = n)) +
  geom_col(fill = "#f9844a", width = 0.7) +
  coord_flip() +                       # 横向条形图，便于展示
  theme(
    panel.background = element_blank(),  # 设置背景为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    panel.border = element_rect(color = "grey20", fill = NA, linewidth = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    axis.title = element_blank(),
    #axis.title = element_text( size = 12, color = "black"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 10, color = "black"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5),  # 横轴标签旋转
    axis.ticks = element_line(color = "grey20", linewidth = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    legend.position = 'none',
    #legend.position = c(0.2,0.75),  # 设置图例位置在右侧
    legend.title = element_blank(), #设置图例标题字体为粗体，大小为12，颜色为深灰色
    legend.text = element_blank(), #设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", linewidth = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )
plot_txss_num

###txss流行率####
library(dplyr)
# 统计每个 class 的 MAG 数量，取前 10
top10_classes <- myxo_genome_metadata %>%
  count(class, name = "n_mags") %>%
  slice_max(n_mags, n = 10) %>%
  pull(class)
# 2.1 计算每个 class 的 MAG 总数（依据 myxo_genome_metadata）
class_mag_counts <- myxo_genome_metadata %>%
  filter(class %in% top10_classes) %>%
  count(class, name = "total_mags")

# 2.2 计算每个 class-type 组合的独特 MAG 数
class_txsstype_mag <- txss_systems %>%
  filter(class %in% top10_classes) %>%
  distinct(class, genome, type) %>%    # 同一 MAG 同一 type 只计一次
  count(class, type, name = "n_mags")

# 2.3 合并并计算流行率
txss_prevalence <- class_txsstype_mag %>%
  left_join(class_mag_counts, by = "class") %>%
  mutate(prevalence = 100 * n_mags / total_mags)

# 2.4 选择要展示的 Type（例如：至少在一个 class 中流行率 > 5% 的 type）
txss_type_to_show <- txss_prevalence %>%
  group_by(type) %>%
  summarise(max_prev = max(prevalence)) %>%
  filter(max_prev > 5) %>%
  pull(type)

txss_heatmap_data <- txss_prevalence %>%
  filter(type %in% txss_type_to_show)

# 2.5 绘制热图
plot_txss_prevalence_heatmap <-
  ggplot(txss_heatmap_data, aes(x = type, y = class, fill = prevalence)) +
  geom_tile(color = "white", lwd = 0.5) +
  scale_fill_gradient(low = "#f7fbff", high = "#f9844a", name = "Prevalence (%)") +
  labs(x = "Defense System Type", y = "Class (Myxococcota)",
       title = "Prevalence of Defense Systems Across Top 10 Classes") +
  # theme_minimal(base_size = 12) +
  # theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
  #       axis.text.y = element_text(size = 10),
  #       legend.position = "right")+
  theme(
    panel.background = element_blank(),  # 设置背景为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    panel.border = element_rect(color = "grey20", fill = NA, linewidth = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    axis.title = element_blank(),
    #axis.title = element_text( size = 12, color = "black"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 10, color = "black"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),  # 横轴标签旋转
    axis.ticks = element_line(color = "grey20", linewidth = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    #legend.position = 'none',
    legend.position = c(0.2,0.75),  # 设置图例位置在右侧
    legend.title = element_blank(), #设置图例标题字体为粗体，大小为12，颜色为深灰色
    #legend.text = element_blank(), #设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", linewidth = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )
plot_txss_prevalence_heatmap

ggsave(plot = plot_txss_prevalence_heatmap, 
       "1_plot_txss_prevalence_heatmap_legend.pdf",
       width = 8, height = 4, units = "in")


#BGC####
##data prepare####

bgc <- read.delim2("all_merged_BGC.tsv")
bgc <- bgc[!grepl("^allorf_", bgc$gene_ID), ]
bgc$gene_ID <- gsub(" ", "", bgc$gene_ID, fixed = TRUE)
#bgc$gene_ID <- gsub("\\s+", "", bgc$gene_ID, perl = TRUE)
#bgc$gene_ID <- gsub("[[:space:]]+", "", bgc$gene_ID)

bgc_gene_gff <- read.delim2("BGC_genes.gff",header = F)
colnames(bgc_gene_gff) <-
  c("contig","source","CDS","start","end","info","strand","info2","detail")
bgc_gene_gff$gene_ID <- str_extract(bgc_gene_gff$detail, "(?<=^ID=)[^;]+")

bgc <- left_join(bgc, bgc_gene_gff, 
                  by = c("gene_ID" = "gene_ID"))


bgc <- bgc %>%
  mutate(
    genome = str_replace(contig, "-\\d+$", "")  
  )
bgc$type <- bgc$category
bgc$hit_position <- sub(".*_", "", bgc$gene_ID)
bgc$hit_position <- as.numeric(bgc$hit_position)
library(dplyr)

bgc_systems <- bgc %>%
  # 1. 按分组字段分组，并按照 hit_position 排序
  group_by(contig, type, product) %>%
  arrange(hit_position, .by_group = TRUE) %>%
  
  # 2. 计算相邻基因 hit_position 的差值，并标记新簇开始
  mutate(
    gap = hit_position - lag(hit_position, default = first(hit_position)),
    new_cluster = row_number() == 1 | gap > 10,
    cluster_id = cumsum(new_cluster)
  ) %>%
  
  # 3. 将簇编号加入分组，重新排序（确保簇内按 hit_position 有序）
  group_by(contig, type, product, cluster_id) %>%
  arrange(hit_position, .by_group = TRUE) %>%
  
  # 4. 汇总每个簇的信息
  summarise(
    type                    = first(type),
    product                 = first(product),
    activity                = "BGC",
    sys_beg                 = first(gene_ID),
    sys_end                 = last(gene_ID),
    protein_in_syst         = paste(gene_ID, collapse = ","),
    genes_count             = n(),
    name_of_profiles_in_sys = paste(gene_functions, collapse = ","),
    .groups = "drop"
  ) %>%
  
  # 5. 选择最终输出列
  select(contig, type, product, activity, sys_beg, sys_end, 
         protein_in_syst, genes_count, name_of_profiles_in_sys)

bgc_systems <- bgc_systems[!is.na(bgc_systems$contig), ]


# 如果 bgc_systems 还没有坐标列，先创建
if(!"bgcbegin" %in% names(bgc_systems)) bgc_systems$bgcbegin <- NA
if(!"bgcend" %in% names(bgc_systems)) bgc_systems$bgcend <- NA

# 创建 hit_id -> start / end 的映射表
start_map <- setNames(bgc$start, bgc$hit_id)
end_map   <- setNames(bgc$end,   bgc$hit_id)
# 填充 sys_beg 对应的起始坐标
bgc_systems$bgcbegin <- ifelse(
  bgc_systems$sys_beg %in% names(start_map),
  start_map[as.character(bgc_systems$sys_beg)],
  bgc_systems$bgcbegin
)
# 填充 sys_end 对应的终止坐标
bgc_systems$bgcend <- ifelse(
  bgc_systems$sys_end %in% names(end_map),
  end_map[as.character(bgc_systems$sys_end)],
  bgc_systems$bgcend
)


# bgc_systems <- bgc_systems %>%
#   mutate(
#     contig = str_replace(sys_beg, "_\\d+$", "")  
#   )
bgc_systems <- bgc_systems %>%
  mutate(
    genome = str_replace(contig, "-\\d+$", "")  
  )

bgc_systems <- 
  bgc_systems %>%
  left_join(myxo_genome_metadata, by = c("genome" = "Genome"))
head(bgc_systems)



##统计及绘图####

###各类bgc总数####
library(ggplot2)

# 统计每个 type 的出现次数（直接按行计）
bgc_type_counts <- bgc_systems %>%
  count(product, sort = TRUE) %>%
  top_n(20, n)              # 保留数量最多的前 50 个 type

library(ggplot2)
plot_bgc_num <-
  ggplot(bgc_type_counts, aes(x = reorder(product, n), y = n)) +
  geom_col(fill = "#00afb9", width = 0.7) +
  coord_flip() +                       # 横向条形图，便于展示
  theme(
    panel.background = element_blank(),  # 设置背景为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    panel.border = element_rect(color = "grey20", fill = NA, linewidth = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    axis.title = element_blank(),
    #axis.title = element_text( size = 12, color = "black"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 10, color = "black"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5),  # 横轴标签旋转
    axis.ticks = element_line(color = "grey20", linewidth = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    legend.position = 'none',
    #legend.position = c(0.2,0.75),  # 设置图例位置在右侧
    legend.title = element_blank(), #设置图例标题字体为粗体，大小为12，颜色为深灰色
    legend.text = element_blank(), #设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", linewidth = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )
plot_bgc_num

###bgc流行率####
library(dplyr)
# 统计每个 class 的 MAG 数量，取前 10
top10_classes <- myxo_genome_metadata %>%
  count(class, name = "n_mags") %>%
  slice_max(n_mags, n = 10) %>%
  pull(class)
# 2.1 计算每个 class 的 MAG 总数（依据 myxo_genome_metadata）
class_mag_counts <- myxo_genome_metadata %>%
  filter(class %in% top10_classes) %>%
  count(class, name = "total_mags")

# 2.2 计算每个 class-type 组合的独特 MAG 数
class_bgctype_mag <- bgc_systems %>%
  filter(class %in% top10_classes) %>%
  distinct(class, genome, type) %>%    # 同一 MAG 同一 type 只计一次
  count(class, type, name = "n_mags")

# 2.3 合并并计算流行率
bgc_prevalence <- class_bgctype_mag %>%
  left_join(class_mag_counts, by = "class") %>%
  mutate(prevalence = 100 * n_mags / total_mags)

# 2.4 选择要展示的 Type（例如：至少在一个 class 中流行率 > 5% 的 type）
bgc_type_to_show <- bgc_prevalence %>%
  group_by(type) %>%
  summarise(max_prev = max(prevalence)) %>%
  filter(max_prev > 1) %>%
  pull(type)

bgc_heatmap_data <- bgc_prevalence %>%
  filter(type %in% bgc_type_to_show)

# 2.5 绘制热图
plot_bgc_prevalence_heatmap <-
  ggplot(bgc_heatmap_data, aes(x = type, y = class, 
                               fill = prevalence)) +
  geom_tile(color = "white", lwd = 0.5) +
  scale_fill_gradient(low = "#f7fbff", high = "#00afb9", 
                      name = "Prevalence (%)") +
  labs(x = "Defense System Type", y = "Class (Myxococcota)",
       title = "Prevalence of BGCs Across Top 10 Classes") +
  # theme_minimal(base_size = 12) +
  # theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
  #       axis.text.y = element_text(size = 10),
  #       legend.position = "right")+
  theme(
    panel.background = element_blank(),  # 设置背景为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    panel.border = element_rect(color = "grey20", fill = NA, linewidth = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    axis.title = element_blank(),
    #axis.title = element_text( size = 12, color = "black"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 10, color = "black"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),  # 横轴标签旋转
    axis.ticks = element_line(color = "grey20", linewidth = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    #legend.position = 'none',
    legend.position = c(0.2,0.75),  # 设置图例位置在右侧
    legend.title = element_blank(), #设置图例标题字体为粗体，大小为12，颜色为深灰色
    #legend.text = element_blank(), #设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", linewidth = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )
plot_bgc_prevalence_heatmap

ggsave(plot = plot_bgc_prevalence_heatmap, 
       "1_plot_bgc_prevalence_heatmap_legend.pdf",
       width = 8, height = 4, units = "in")


#MGE################

##prophage####

genomad_prophage <- read.delim2("merged_myxo_mag_virus_summary.tsv")
# 过滤不包含plasmid的行
genomad_prophage <- genomad_prophage %>%
  filter(!grepl("plasmid", seq_name, ignore.case = TRUE))
genomad_prophage <- genomad_prophage[genomad_prophage$topology == "Provirus",]

# 提取 | 前面的部分
genomad_prophage$contig <- sub("\\|.*$", "", genomad_prophage$seq_name)
genomad_prophage$genome <- gsub("-.*","",genomad_prophage$contig)
genomad_prophage <- genomad_prophage %>%
  mutate(
    position = str_extract(seq_name, "(?<=provirus_)\\d+_\\d+"),       # 提取"数字_数字"格式的位置信息
    begin = as.numeric(str_extract(position, "^\\d+")),    # 提取_前的数字作为起始位置
    end = as.numeric(str_extract(position, "\\d+$"))       # 提取_后的数字作为结束位置
  ) %>% 
  select(-position)

# genomad_prophage <- genomad_prophage %>%
#   mutate(
#     begin = ifelse(topology != "Provirus", 1, begin),
#     end = ifelse(topology != "Provirus", length, end)
#   )

genomad_prophage$MGE_type <- "Prophage"

genomad_prophage$start <- genomad_prophage$begin
Prophage_n <- select(genomad_prophage, genome, contig, start, end, length, MGE_type)


##Plasmid####
plasmid_info <- read.delim2("merged_myxo_mag_plasmid_summary.tsv")
plasmid_info <- plasmid_info %>%
  mutate(
    genome = sub("(.*)-.*", "\\1", plasmid_info$seq_name)
  )
plasmid_info$contig <- plasmid_info$seq_name
plasmid_info$MGE_type <- "Plasmid"
plasmid_info$start <- 1
plasmid_info$end <- plasmid_info$length
plasmid_n <- select(plasmid_info, genome, contig, start, end, length, MGE_type)


##IS####
IS <- read.delim2("3229MAGs_filtered_completeIS.tsv")
IS <- IS %>% 
  mutate(
    genome = sub("(.*)-.*", "\\1", IS$seqID),
    contig = seqID
  )

IS$MGE_type <- "IS"
IS$start <- IS$isBegin
IS$end <- IS$isEnd
IS$length <- IS$isLen
IS_n <- select(IS, genome, contig, start, end, length, MGE_type)

##integron####
integrons <- read.delim2("merged_myxo_mag.integrons")
complete_integrons <- integrons[integrons$type == "complete",]
# 处理数据
complete_integrons_clean <- complete_integrons %>%
  group_by(ID_replicon) %>%          # 按ID_replicon分组
  summarise(
    min_pos_beg = min(pos_beg),      # 最小起始位置
    max_pos_end = max(pos_end),      # 最大结束位置
    length = max_pos_end - min_pos_beg +1,
    .groups = "drop"                 # 取消分组
  ) %>%
  arrange(ID_replicon)               # 按ID排序

complete_integrons_clean <- complete_integrons_clean %>%
  mutate(
    genome = sub("(.*)-.*", "\\1", ID_replicon),
    contig = ID_replicon
  )

complete_integrons_clean$MGE_type <- "Integron"
complete_integrons_clean$start <- complete_integrons_clean$min_pos_beg
complete_integrons_clean$end <- complete_integrons_clean$max_pos_end
integron_n <- select(complete_integrons_clean, genome, contig, start, end, length, MGE_type)


##ICE#####
ice <- read.delim2("processed_ice_sequences_ID.txt",header = F)
colnames(ice) <- c("seq_name","range","length")

ice$seq_name <- gsub("\\.","_",ice$seq_name)
ice$seq_name <- gsub("_contig_","-",ice$seq_name)
ice <- separate(ice, range, into = c("num1", "num2"), sep = "-", convert = TRUE) %>%
  mutate(
    start = pmin(num1, num2),  # 取较小值作为 start
    end = pmax(num1, num2)     # 取较大值作为 end
  ) %>%
  select(-num1, -num2)  # 移除中间列

ice$length <- ice$end - ice$start + 1
ice$contig <- ice$seq_name
ice$genome <- sub("(.*)-.*", "\\1", ice$contig)

ice$MGE_type <- "ICEs"
ice$contig <- ice$seq_name
ice_n <- select(ice, genome, contig, start, end, length, MGE_type)

##composite_transposons####
# 按基因组和contig分组，排序
IS_sorted <- IS_n %>%
  arrange(genome, contig, start)

# 找出邻近的IS对
find_composite_transposons <- function(df, max_distance = 100000) {
  results <- data.frame()
  
  # 按contig分组
  for(ctg in unique(df$contig)) {
    ctg_data <- df %>% filter(contig == ctg) %>% arrange(start)
    
    if(nrow(ctg_data) >= 2) {
      for(i in 1:(nrow(ctg_data)-1)) {
        for(j in (i+1):nrow(ctg_data)) {
          # 计算两个IS之间的距离
          distance <- ctg_data$start[j] - ctg_data$end[i]
          
          # 如果距离在阈值内，可能是复合转座子
          if(distance > 0 & distance <= max_distance) {
            # 确定方向
            direction <- ifelse(
              ctg_data$start[i] < ctg_data$start[j] & 
                ctg_data$end[i] < ctg_data$end[j],
              "same", "opposite"
            )
            
            results <- rbind(results, data.frame(
              genome = ctg_data$genome[i],
              contig = ctg,
              IS1_start = ctg_data$start[i],
              IS1_end = ctg_data$end[i],
              IS2_start = ctg_data$start[j],
              IS2_end = ctg_data$end[j],
              distance = distance,
              direction = direction,
              total_length = (ctg_data$end[j] - ctg_data$start[i]),
              IS1_length = ctg_data$length[i],
              IS2_length = ctg_data$length[j]
            ))
          }
        }
      }
    }
  }
  return(results)
}

# 应用函数
composite_candidates <- find_composite_transposons(IS_n, max_distance = 100000)
composite_transposons_n <- 
  select(composite_candidates, genome, contig, IS1_start, IS2_end, total_length)
composite_transposons_n$MGE_type <- "Composite_transposon"
colnames(composite_transposons_n) <- c("genome","contig","start","end","length","MGE_type")

##merged_all_MGE####

all_MGE <- rbind(plasmid_n, Prophage_n, IS_n, integron_n, ice_n, composite_transposons_n)

# 假设两个数据框的匹配列都叫 "genome"
all_MGE <- left_join(all_MGE, 
                     select(myxo_genome_metadata, Genome, class),
                     by = c("genome"="Genome"))

mge_wide <- all_MGE %>%
  count(genome, MGE_type, name = "count") %>%
  pivot_wider(
    id_cols = genome,
    names_from = MGE_type,
    values_from = count,
    values_fill = 0  # 关键：将NA填充为0
  )
myxo_genome_metadata <- left_join(myxo_genome_metadata, 
                                mge_wide,
                                by = c("Genome"="genome"))
myxo_genome_metadata <- myxo_genome_metadata %>%
  mutate(across(17:22, ~ replace(., is.na(.), 0)))

write.table(myxo_genome_metadata,
            "myxo_genome_metadata_withMGE.tsv",
            sep = "\t",           # TSV使用制表符分隔
            row.names = FALSE,    # 不写入行名
            col.names = TRUE,     # 写入列名
            quote = FALSE        # 不添加引号
)            # NA值的表示方式

mge_matrix <- select(myxo_genome_metadata, Genome, Prophage, 
                     Plasmid, ICEs, Integron, IS,Composite_transposon, class)
mge_long <- mge_matrix %>%
  pivot_longer(
    cols = c(Prophage, Plasmid, ICEs, Integron, IS,Composite_transposon),  # 指定要转换的列
    names_to = "MGE_type",       # 新列名：存储原来的列名
    values_to = "count"          # 新列名：存储原来的值
  )

###绘图:MGE数量热图####
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(scales)

# 1. 按Tax和MGE_type汇总
tax_summary <- mge_long %>%
  group_by(class, MGE_type) %>%
  summarise(
    mean_count = mean(count),
    median_count = median(count),
    total_count = sum(count),
    n_genomes = n(),
    .groups = "drop"
  )
# tax_summary$class <- factor(tax_summary$class, 
#                             levels = c("Ef","Sa","Kp","Ab","Pa","En"))
tax_summary <- tax_summary[tax_summary$class %in% top10_classes,]
# 2. 绘制聚合热图
p_mge_num_heatmap <- 
  ggplot(tax_summary, aes(x = MGE_type, y = class)) +
  geom_tile(aes(fill = mean_count), color = "white", size = 1) +
  geom_text(aes(label = round(mean_count, 1)), 
            color = "black", size = 3, fontface = "plain") +
  # scale_fill_gradient2(
  #   low = "white",
  #   mid = "#1982c4",
  #   high = "#ff9f1c",
  #   midpoint = 2,
  #   name = "Mean Count"
  # ) +
  scale_fill_gradient2(
    low = "#F5F5F5",
    mid = "#6C91BF",
    high = "#C44E52",
    midpoint = 2,
    name = "Mean Count"
  )+
  theme(
    panel.background = element_blank(),  # 设置背景为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    panel.border = element_rect(color = "grey20", fill = NA, linewidth = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    axis.title = element_blank(),
    #axis.title = element_text( size = 12, color = "black"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 10, color = "black"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),  # 横轴标签旋转
    axis.ticks = element_line(color = "grey20", linewidth = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    legend.position = 'none',
    #legend.position = c(0.2,0.75),  # 设置图例位置在右侧
    legend.title = element_blank(), #设置图例标题字体为粗体，大小为12，颜色为深灰色
    legend.text = element_blank(), #设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", linewidth = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )
p_mge_num_heatmap

####绘图###
# library(rstatix)
# md_p_val2 <- mge_long %>%
#   # 先过滤掉全为0的分组组合
#   group_by(Tax, MGE_type) %>%
#   filter(!all(count == 0)) %>%
#   ungroup() %>%
#   # 再进行检验
#   group_by(MGE_type) %>%
#   wilcox_test(count ~ Tax) %>%
#   adjust_pvalue(p.col = "p", method = "bonferroni") %>%
#   add_significance(p.col = "p.adj") %>%
#   add_xy_position(x = "Tax", dodge = 0.5, group = "MGE_type")
# 
# library(ggplot2)
# library(ggsignif)#geom_signif
# library(ggpubr)#stat_pvalue_manual stat_compare_means
# library(gghalves)
# library(ggsci)
# library(ggbeeswarm)
# p_mge_tax_box <- ggplot(mge_long, 
#                         aes(x = Tax, y = count)) +
#   geom_quasirandom(
#     aes(
#       x = Tax,
#       y = count,
#       color = Tax
#     ),
#     dodge.width = 0,
#     width = 0.25,   # 控制横向展开
#     varwidth = FALSE,
#     size = 0.5
#   )+
#   geom_boxplot(
#     aes(fill = Tax),
#     color = "black",
#     linewidth = 0.4,
#     alpha = 0.5,          # ← 整体透明度
#     outlier.shape = NA,
#     width = 0.8
#   )+
#   
#   #  stat_pvalue_manual(md_p_val2,label = "p.adj.signif",label.size=5,hide.ns = T)+
#   facet_wrap(.~MGE_type,nrow=5, scales = "free_y")+ #分⾯操作
#   # 定义颜⾊
#   # scale_fill_npg()+
#   # scale_color_npg()+
#   scale_fill_manual(values = c("Ef"="#9f86c0","Sa"="#FFBE7D",
#                                "Kp"="#76c893","Ab"="#F0E442","Pa"="#ff99c8",
#                                "En"="#56B4E9")
#   )+
#   scale_color_manual(values = c("Ef"="#9f86c0","Sa"="#FFBE7D",
#                                 "Kp"="#76c893","Ab"="#F0E442","Pa"="#ff99c8",
#                                 "En"="#56B4E9")
#   )+
#   labs(x="River",y="log10 (DS Density +1) × 1e-3")+
#   theme(
#     panel.spacing = unit(0.1, "lines"),  # 将分面间距设为0
#     strip.background = element_blank(),      # 移除分面背景
#     strip.text = element_blank(),            # 移除分面标题文字
#     panel.background = element_blank(),  # 设置背景为透明
#     # panel.grid.major = element_blank(),  # 去掉主要网格线
#     panel.grid.major = element_line(color = "gray90", linewidth = 0.2),
#     panel.grid.minor = element_blank(),  # 去掉次要网格线
#     panel.border = element_rect(color = "grey20", fill = NA, size = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
#     axis.title = element_blank(),
#     #axis.title = element_text( size = 12, color = "black"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
#     axis.text = element_text(size = 10, color = "black"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
#     axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5),  # 横轴标签旋转
#     axis.ticks = element_line(color = "grey20", size = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
#     legend.position = 'none',
#     #legend.position = c(0.2,0.75),  # 设置图例位置在右侧
#     legend.title = element_blank(), #设置图例标题字体为粗体，大小为12，颜色为深灰色
#     legend.text = element_blank(), #设置图例文本字体大小为10，颜色为深灰色
#     legend.background = element_blank(),  # 设置图例背景为透明
#     legend.box.background = element_rect(color = "grey80", linewidth = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
#     legend.key = element_blank()  # 去掉图例键的背景
#   )
# p_mge_tax_box
# 
# ggsave(plot = p_mge_tax_box,'1_p_mge_tax_box.pdf',width = 2.28,height = 4.37,units = 'in')


#MGE携带的功能基因####
##MGE与DS关系####

###统计质粒上的DS####
head(plasmid_n)
head(DS)

plasmid_DS <- DS %>%
  inner_join(plasmid_n, 
             by = c("contig" = "contig")) %>%
  filter(
    DSbegin >= start,
    DSend <= end
  )
plasmid_DS <- select(plasmid_DS, genome.x, type, subtype, MGE_type)

plasmid_DS <- 
  left_join(plasmid_DS, 
  select(myxo_genome_metadata, Genome, class),
  by = c("genome.x"="Genome"))

###prophage上的DS####

prophage_DS <- DS %>%
  inner_join(Prophage_n, 
             by = c("contig" = "contig")) %>%
  filter(
    DSbegin >= start,
    DSend <= end
  )

prophage_DS <- select(prophage_DS, genome.x, type, subtype, MGE_type)

prophage_DS <- 
  left_join(prophage_DS, 
            select(myxo_genome_metadata, Genome, class),
            by = c("genome.x"="Genome"))
###is上的DS####

IS_DS <- DS %>%
  inner_join(IS_n, 
             by = c("contig" = "contig")) %>%
  filter(
    DSbegin >= start,
    DSend <= end
  )

IS_DS <- select(IS_DS, genome.x, type, subtype, MGE_type)

IS_DS <- 
  left_join(IS_DS, 
            select(myxo_genome_metadata, Genome, class),
            by = c("genome.x"="Genome"))
###composite_transposons上的DS####

composite_transposons_DS <- 
  DS %>%
  inner_join(composite_transposons_n, 
             by = c("contig" = "contig")) %>%
  filter(
    DSbegin >= start,
    DSend <= end
  )

composite_transposons_DS <- 
  select(composite_transposons_DS, genome.x, type, subtype, MGE_type)

composite_transposons_DS <- 
  left_join(composite_transposons_DS, 
            select(myxo_genome_metadata, Genome, class),
            by = c("genome.x"="Genome"))



###Integron上的DS####

integron_DS <- DS %>%
  inner_join(integron_n, 
             by = c("contig" = "contig")) %>%
  filter(
    DSbegin >= start,
    DSend <= end
  )

integron_DS <- select(integron_DS, genome.x, type, subtype, MGE_type)

integron_DS <- left_join(integron_DS, 
                         select(myxo_genome_metadata, Genome, class),
                         by = c("genome.x"="Genome"))
###ICE上的DS####

ice_DS <- DS %>%
  inner_join(integron_n, 
             by = c("contig" = "contig")) %>%
  filter(
    DSbegin >= start,
    DSend <= end
  )

ice_DS <- select(ice_DS, genome.x, type, subtype, MGE_type)

ice_DS <- 
  left_join(ice_DS, 
  select(myxo_genome_metadata, Genome, class),
  by = c("genome.x"="Genome"))

##MGE_DS 归一化的点热图####
all_mge_ds <- 
  rbind(plasmid_DS, prophage_DS, 
        composite_transposons_DS, ice_DS, integron_DS)
all_mge_ds <- 
  all_mge_ds[,-1]


# 1. 计算计数
mge_ds_dot_heatmap <- all_mge_ds %>%
  count(type, MGE_type, class, name = "count")

# mge_ds_dot_heatmap$class <- 
#   factor(dot_heatmap_ds$Tax, levels = c("Ef","Sa","Kp","Ab","Pa","En"))
mge_ds_dot_heatmap <- 
  mge_ds_dot_heatmap[mge_ds_dot_heatmap$class %in% top10_classes,]


mge_ds_type_totals <- mge_ds_dot_heatmap %>%
  group_by(type) %>%
  summarise(total_count = sum(count, na.rm = TRUE))
mge_ds_type_totals <- mge_ds_type_totals[mge_ds_type_totals$total_count >10,]

mge_ds_dot_heatmap <- 
  mge_ds_dot_heatmap[mge_ds_dot_heatmap$type %in% mge_ds_type_totals$type,]


# 添加归一化列
mge_ds_dot_heatmap_normalized <- mge_ds_dot_heatmap %>%
  left_join(tax_summary, by = c("MGE_type", "class")) %>%
  mutate(
    count_normalized = count / total_count,
    percent = (count / total_count) * 100  # 可选：转换为百分比
  ) %>%
  arrange(type, MGE_type)

mge_ds_dot_heatmap_normalized$MGE_type[mge_ds_dot_heatmap_normalized$MGE_type == "Composite_transposon"] <- "IS_Tn"


# 3. 创建分面热图
p_dot_heatmap_mge_ds_normalized <-
  ggplot(mge_ds_dot_heatmap_normalized, aes(x = class, y = type)) +
  geom_point(aes(size = count_normalized, color = MGE_type), alpha = 1) +
  scale_size_continuous(
    name = "Count",
    range = c(1, 10),  # 控制点的大小范围
    breaks = pretty_breaks(n = 5)
  ) +
  scale_color_manual(values = c("IS"="#9f86c0","Plasmid"="#FFBE7D",
                                "Integron"="#76c893","ICEs"="#F0E442",
                                "IS_Tn"="#ff99c8",
                                "Prophage"="#56B4E9")
  )+
  facet_grid(~ MGE_type, scales = "free", space = "free") +
  #facet_wrap(~ MGE_type, scales = "free", ncol = 5) +
  theme(
    panel.spacing = unit(0.1, "lines"),  # 将分面间距设为0
    panel.background = element_blank(),  # 设置背景为透明
    # panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.major = element_line(color = "gray90", linewidth = 0.2),
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    panel.border = element_rect(color = "grey20", fill = NA, linewidth = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    axis.title = element_blank(),
    #axis.title = element_text( size = 12, color = "black"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 10, color = "black"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),  # 横轴标签旋转
    axis.ticks = element_line(color = "grey20", linewidth = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    legend.position = 'none',
    #legend.position = c(0.2,0.75),  # 设置图例位置在右侧
    legend.title = element_blank(), #设置图例标题字体为粗体，大小为12，颜色为深灰色
    legend.text = element_blank(), #设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", linewidth = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )
p_dot_heatmap_mge_ds_normalized

p_dot_heatmap_mge_ds_normalized_legend <-
  ggplot(mge_ds_dot_heatmap_normalized, aes(x = class, y = type)) +
  geom_point(aes(size = count_normalized, color = MGE_type), alpha = 1) +
  scale_size_continuous(
    name = "Count",
    range = c(1, 10),  # 控制点的大小范围
    breaks = pretty_breaks(n = 5)
  ) +
  scale_color_manual(values = c("IS"="#9f86c0","Plasmid"="#FFBE7D",
                                "Integron"="#76c893","ICEs"="#F0E442",
                                "Composite_transposon"="#ff99c8",
                                "Prophage"="#56B4E9")
  )+
  facet_grid(~ MGE_type, scales = "free", space = "free") +
  #facet_wrap(~ MGE_type, scales = "free", ncol = 5) +
  theme(
    panel.spacing = unit(0.1, "lines"),  # 将分面间距设为0
    panel.background = element_blank(),  # 设置背景为透明
    # panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.major = element_line(color = "gray90", linewidth = 0.2),
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    panel.border = element_rect(color = "grey20", fill = NA, linewidth = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    axis.title = element_blank(),
    #axis.title = element_text( size = 12, color = "black"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 10, color = "black"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),  # 横轴标签旋转
    axis.ticks = element_line(color = "grey20", linewidth = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    #legend.position = 'none',
    legend.position = c(0.2,0.75),  # 设置图例位置在右侧
    #legend.title = element_blank(), #设置图例标题字体为粗体，大小为12，颜色为深灰色
    #legend.text = element_blank(), #设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", linewidth = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )
p_dot_heatmap_mge_ds_normalized_legend

ggsave(plot = p_dot_heatmap_mge_ds_normalized_legend,
       '1_p_dot_heatmap_mge_ds_normalized_legend.pdf',
       width = 6.9,height = 5.5,units = 'in')

##MGE与txss关系####

###统计质粒上的ss####
head(plasmid_n)
head(txss_systems)

plasmid_txss <- txss_systems %>%
  inner_join(plasmid_n, 
             by = c("contig" = "contig")) %>%
  filter(
    SSbegin >= start,
    SSend <= end
  )

plasmid_txss <- 
  select(plasmid_txss, genome.x, type, subtype, MGE_type)

plasmid_txss <- 
  left_join(plasmid_txss, 
            select(myxo_genome_metadata, Genome, class),
            by = c("genome.x"="Genome"))

###prophage上的txss####

prophage_txss <- txss_systems %>%
  inner_join(Prophage_n, 
             by = c("contig" = "contig")) %>%
  filter(
    SSbegin >= start,
    SSend <= end
  )

prophage_txss <- 
  select(prophage_txss, genome.x, type, subtype, MGE_type)

prophage_txss <- 
  left_join(prophage_txss, 
            select(myxo_genome_metadata, Genome, class),
            by = c("genome.x"="Genome"))
###is上的txss####

IS_txss <- txss_systems %>%
  inner_join(IS_n, 
             by = c("contig" = "contig")) %>%
  filter(
    SSbegin >= start,
    SSend <= end
  )

IS_txss <- 
  select(IS_txss, genome.x, type, subtype, MGE_type)

IS_txss <- 
  left_join(IS_txss, 
            select(myxo_genome_metadata, Genome, class),
            by = c("genome.x"="Genome"))
###composite_transposons上的txss####

composite_transposons_txss <- 
  txss_systems %>%
  inner_join(composite_transposons_n, 
             by = c("contig" = "contig")) %>%
  filter(
    SSbegin >= start,
    SSend <= end
  )

composite_transposons_txss <- 
  select(composite_transposons_txss, genome.x, type, subtype, MGE_type)

composite_transposons_txss <- 
  left_join(composite_transposons_txss, 
            select(myxo_genome_metadata, Genome, class),
            by = c("genome.x"="Genome"))



###Integron上的txss####

integron_txss <- txss_systems %>%
  inner_join(integron_n, 
             by = c("contig" = "contig")) %>%
  filter(
    SSbegin >= start,
    SSend <= end
  )

integron_txss <- select(integron_txss, genome.x, type, subtype, MGE_type)

integron_txss <- left_join(integron_txss, 
                         select(myxo_genome_metadata, Genome, class),
                         by = c("genome.x"="Genome"))
###ICE上的txss####

ice_txss <- txss_systems %>%
  inner_join(integron_n, 
             by = c("contig" = "contig")) %>%
  filter(
    SSbegin >= start,
    SSend <= end
  )

ice_txss <- select(ice_txss, genome.x, type, subtype, MGE_type)

ice_txss <- 
  left_join(ice_txss, 
            select(myxo_genome_metadata, Genome, class),
            by = c("genome.x"="Genome"))

###MGE_txss 归一化的点热图####
all_mge_txss <- 
  rbind(plasmid_txss, prophage_txss, 
        composite_transposons_txss, ice_txss, integron_txss)
all_mge_txss <- 
  all_mge_txss[,-1]


# 1. 计算计数
mge_txss_dot_heatmap <- all_mge_txss %>%
  count(type, MGE_type, class, name = "count")

# mge_ds_dot_heatmap$class <- 
#   factor(dot_heatmap_ds$Tax, levels = c("Ef","Sa","Kp","Ab","Pa","En"))
mge_txss_dot_heatmap <- 
  mge_txss_dot_heatmap[mge_txss_dot_heatmap$class %in% top10_classes,]


mge_txss_type_totals <- mge_txss_dot_heatmap %>%
  group_by(type) %>%
  summarise(total_count = sum(count, na.rm = TRUE))
mge_txss_type_totals <- 
  mge_txss_type_totals[mge_txss_type_totals$total_count >=1,]

mge_txss_dot_heatmap <- 
  mge_txss_dot_heatmap[mge_txss_dot_heatmap$type %in% mge_txss_type_totals$type,]


# 添加归一化列
mge_txss_dot_heatmap_normalized <- 
  mge_txss_dot_heatmap %>%
  left_join(tax_summary, by = c("MGE_type", "class")) %>%
  mutate(
    count_normalized = count / total_count,
    percent = (count / total_count) * 100  # 可选：转换为百分比
  ) %>%
  arrange(type, MGE_type)
mge_txss_dot_heatmap_normalized$MGE_type[mge_txss_dot_heatmap_normalized$MGE_type == "Composite_transposon"] <- "IS_Tn"

# 3. 创建分面热图
p_dot_heatmap_mge_txss_normalized <-
  ggplot(mge_txss_dot_heatmap_normalized, aes(x = class, y = type)) +
  geom_point(aes(size = count_normalized, color = MGE_type), alpha = 1) +
  scale_size_continuous(
    name = "Count",
    range = c(1, 10),  # 控制点的大小范围
    breaks = pretty_breaks(n = 5)
  ) +
  scale_color_manual(values = c("IS"="#9f86c0","Plasmid"="#FFBE7D",
                                "Integron"="#76c893","ICEs"="#F0E442",
                                "IS_Tn"="#ff99c8",
                                "Prophage"="#56B4E9")
  )+
  facet_grid(~ MGE_type, scales = "free", space = "free") +
  #facet_wrap(~ MGE_type, scales = "free", ncol = 5) +
  theme(
    panel.spacing = unit(0.1, "lines"),  # 将分面间距设为0
    panel.background = element_blank(),  # 设置背景为透明
    # panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.major = element_line(color = "gray90", linewidth = 0.2),
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    panel.border = element_rect(color = "grey20", fill = NA, linewidth = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    axis.title = element_blank(),
    #axis.title = element_text( size = 12, color = "black"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 10, color = "black"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),  # 横轴标签旋转
    axis.ticks = element_line(color = "grey20", linewidth = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    legend.position = 'none',
    #legend.position = c(0.2,0.75),  # 设置图例位置在右侧
    legend.title = element_blank(), #设置图例标题字体为粗体，大小为12，颜色为深灰色
    legend.text = element_blank(), #设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", linewidth = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )
p_dot_heatmap_mge_txss_normalized

p_dot_heatmap_mge_txss_normalized_legend <-
  ggplot(mge_txss_dot_heatmap_normalized, aes(x = class, y = type)) +
  geom_point(aes(size = count_normalized, color = MGE_type), alpha = 1) +
  scale_size_continuous(
    name = "Count",
    range = c(1, 10),  # 控制点的大小范围
    breaks = pretty_breaks(n = 5)
  ) +
  scale_color_manual(values = c("IS"="#9f86c0","Plasmid"="#FFBE7D",
                                "Integron"="#76c893","ICEs"="#F0E442",
                                "Composite_transposon"="#ff99c8",
                                "Prophage"="#56B4E9")
  )+
  facet_grid(~ MGE_type, scales = "free", space = "free") +
  #facet_wrap(~ MGE_type, scales = "free", ncol = 5) +
  theme(
    panel.spacing = unit(0.1, "lines"),  # 将分面间距设为0
    panel.background = element_blank(),  # 设置背景为透明
    # panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.major = element_line(color = "gray90", linewidth = 0.2),
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    panel.border = element_rect(color = "grey20", fill = NA, linewidth = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    axis.title = element_blank(),
    #axis.title = element_text( size = 12, color = "black"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 10, color = "black"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),  # 横轴标签旋转
    axis.ticks = element_line(color = "grey20", linewidth = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    #legend.position = 'none',
    legend.position = c(0.2,0.75),  # 设置图例位置在右侧
    #legend.title = element_blank(), #设置图例标题字体为粗体，大小为12，颜色为深灰色
    #legend.text = element_blank(), #设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", linewidth = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )
p_dot_heatmap_mge_txss_normalized_legend

ggsave(plot = p_dot_heatmap_mge_txss_normalized_legend,
       '1_p_dot_heatmap_mge_txss_normalized_legend.pdf',
       width = 6.9,height = 5.5,units = 'in')

##MGE与AMP关系####

###统计质粒上的amp####
head(plasmid_n)
head(AMP)
AMP <- AMP %>%
  rename(ampbegin = start,
         ampend = end)

plasmid_amp <- 
  AMP %>%
  inner_join(plasmid_n, 
             by = c("contig" = "contig")) %>%
  filter(
    ampbegin >= start,
    ampend <= end
  )

plasmid_amp <- 
  select(plasmid_amp, genome.x, AMP_family, Hemolytic, MGE_type)

plasmid_amp <- 
  left_join(plasmid_amp, 
            select(myxo_genome_metadata, Genome, class),
            by = c("genome.x"="Genome"))

###prophage上的amp####

prophage_amp <- AMP %>%
  inner_join(Prophage_n, 
             by = c("contig" = "contig")) %>%
  filter(
    ampbegin >= start,
    ampend <= end
  )

prophage_amp <- 
  select(prophage_amp, genome.x, AMP_family, Hemolytic, MGE_type)

prophage_amp <- 
  left_join(prophage_amp, 
            select(myxo_genome_metadata, Genome, class),
            by = c("genome.x"="Genome"))
###is上的amp####

IS_amp <- AMP %>%
  inner_join(IS_n, 
             by = c("contig" = "contig")) %>%
  filter(
    ampbegin >= start,
    ampend <= end
  )

IS_amp <- 
  select(IS_amp, genome.x, AMP_family, Hemolytic, MGE_type)

IS_amp <- 
  left_join(IS_amp, 
            select(myxo_genome_metadata, Genome, class),
            by = c("genome.x"="Genome"))
###composite_transposons上的amp####

composite_transposons_amp <- 
  AMP %>%
  inner_join(composite_transposons_n, 
             by = c("contig" = "contig")) %>%
  filter(
    ampbegin >= start,
    ampend <= end
  )

composite_transposons_amp <- 
  select(composite_transposons_amp, genome.x, AMP_family, Hemolytic, MGE_type)

composite_transposons_amp <- 
  left_join(composite_transposons_amp, 
            select(myxo_genome_metadata, Genome, class),
            by = c("genome.x"="Genome"))



###Integron上的amp####

integron_amp <- AMP %>%
  inner_join(integron_n, 
             by = c("contig" = "contig")) %>%
  filter(
    ampbegin >= start,
    ampend <= end
  )

integron_amp <- select(integron_amp, genome.x, AMP_family, Hemolytic, MGE_type)

integron_amp <- left_join(integron_amp, 
                           select(myxo_genome_metadata, Genome, class),
                           by = c("genome.x"="Genome"))
###ICE上的amp####

ice_amp <- AMP %>%
  inner_join(integron_n, 
             by = c("contig" = "contig")) %>%
  filter(
    ampbegin >= start,
    ampend <= end
  )

ice_amp <- select(ice_amp, genome.x, AMP_family, Hemolytic, MGE_type)

ice_amp <- 
  left_join(ice_amp, 
            select(myxo_genome_metadata, Genome, class),
            by = c("genome.x"="Genome"))

###MGE_amp 归一化的点热图####
all_mge_amp <- 
  rbind(plasmid_amp, prophage_amp, 
        composite_transposons_amp, ice_amp, integron_amp)
all_mge_amp <- 
  all_mge_amp[,-1]


# 1. 计算计数
mge_amp_dot_heatmap <- all_mge_amp %>%
  count(AMP_family, MGE_type, class, name = "count")

# mge_ds_dot_heatmap$class <- 
#   factor(dot_heatmap_ds$Tax, levels = c("Ef","Sa","Kp","Ab","Pa","En"))
mge_amp_dot_heatmap <- 
  mge_amp_dot_heatmap[mge_amp_dot_heatmap$class %in% top10_classes,]


mge_amp_type_totals <- mge_amp_dot_heatmap %>%
  group_by(AMP_family) %>%
  summarise(total_count = sum(count, na.rm = TRUE))
mge_amp_type_totals <- 
  mge_amp_type_totals[mge_amp_type_totals$total_count >=1,]

mge_amp_dot_heatmap <- 
  mge_amp_dot_heatmap[mge_amp_dot_heatmap$AMP_family %in% mge_amp_type_totals$AMP_family,]


# 添加归一化列
mge_amp_dot_heatmap_normalized <- 
  mge_amp_dot_heatmap %>%
  left_join(tax_summary, by = c("MGE_type", "class")) %>%
  mutate(
    count_normalized = count / total_count,
    percent = (count / total_count) * 100  # 可选：转换为百分比
  ) %>%
  arrange(AMP_family, MGE_type)
mge_amp_dot_heatmap_normalized$MGE_type[mge_amp_dot_heatmap_normalized$MGE_type == "Composite_transposon"] <- "IS_Tn"

# 3. 创建分面热图
p_dot_heatmap_mge_amp_normalized <-
  ggplot(mge_amp_dot_heatmap_normalized, aes(x = class, y = AMP_family)) +
  geom_point(aes(size = count_normalized, color = MGE_type), alpha = 1) +
  scale_size_continuous(
    name = "Count",
    range = c(1, 10),  # 控制点的大小范围
    breaks = pretty_breaks(n = 5)
  ) +
  scale_color_manual(values = c("IS"="#9f86c0","Plasmid"="#FFBE7D",
                                "Integron"="#76c893","ICEs"="#F0E442",
                                "IS_Tn"="#ff99c8",
                                "Prophage"="#56B4E9")
  )+
  facet_grid(~ MGE_type, scales = "free", space = "free") +
  #facet_wrap(~ MGE_type, scales = "free", ncol = 5) +
  theme(
    panel.spacing = unit(0.1, "lines"),  # 将分面间距设为0
    panel.background = element_blank(),  # 设置背景为透明
    # panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.major = element_line(color = "gray90", linewidth = 0.2),
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    panel.border = element_rect(color = "grey20", fill = NA, linewidth = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    axis.title = element_blank(),
    #axis.title = element_text( size = 12, color = "black"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 10, color = "black"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),  # 横轴标签旋转
    axis.ticks = element_line(color = "grey20", linewidth = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    legend.position = 'none',
    #legend.position = c(0.2,0.75),  # 设置图例位置在右侧
    legend.title = element_blank(), #设置图例标题字体为粗体，大小为12，颜色为深灰色
    legend.text = element_blank(), #设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", linewidth = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )
p_dot_heatmap_mge_amp_normalized

p_dot_heatmap_mge_amp_normalized_legend <-
  ggplot(mge_amp_dot_heatmap_normalized, aes(x = class, y = AMP_family)) +
  geom_point(aes(size = count_normalized, color = MGE_type), alpha = 1) +
  scale_size_continuous(
    name = "Count",
    range = c(1, 10),  # 控制点的大小范围
    breaks = pretty_breaks(n = 5)
  ) +
  scale_color_manual(values = c("IS"="#9f86c0","Plasmid"="#FFBE7D",
                                "Integron"="#76c893","ICEs"="#F0E442",
                                "Composite_transposon"="#ff99c8",
                                "Prophage"="#56B4E9")
  )+
  facet_grid(~ MGE_type, scales = "free", space = "free") +
  #facet_wrap(~ MGE_type, scales = "free", ncol = 5) +
  theme(
    panel.spacing = unit(0.1, "lines"),  # 将分面间距设为0
    panel.background = element_blank(),  # 设置背景为透明
    # panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.major = element_line(color = "gray90", linewidth = 0.2),
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    panel.border = element_rect(color = "grey20", fill = NA, linewidth = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    axis.title = element_blank(),
    #axis.title = element_text( size = 12, color = "black"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 10, color = "black"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),  # 横轴标签旋转
    axis.ticks = element_line(color = "grey20", linewidth = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    #legend.position = 'none',
    legend.position = c(0.2,0.75),  # 设置图例位置在右侧
    #legend.title = element_blank(), #设置图例标题字体为粗体，大小为12，颜色为深灰色
    #legend.text = element_blank(), #设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", linewidth = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )
p_dot_heatmap_mge_amp_normalized_legend

ggsave(plot = p_dot_heatmap_mge_amp_normalized_legend,
       '1_p_dot_heatmap_mge_amp_normalized_legend.pdf',
       width = 6.9,height = 5.5,units = 'in')

##MGE与bgc关系####

###统计质粒上的bgc####
head(plasmid_n)
head(bgc_systems)

plasmid_bgc <- bgc_systems %>%
  inner_join(plasmid_n, 
             by = c("contig" = "contig")) %>%
  filter(
    bgcbegin >= start,
    bgcend <= end
  )

plasmid_bgc <- 
  select(plasmid_bgc, 
         genome.x, type, product, MGE_type)

plasmid_bgc <- 
  left_join(plasmid_bgc, 
            select(myxo_genome_metadata, 
                   Genome, class),
            by = c("genome.x"="Genome"))

###prophage上的bgc####

prophage_bgc <- bgc_systems %>%
  inner_join(Prophage_n, 
             by = c("contig" = "contig")) %>%
  filter(
    bgcbegin >= start,
    bgcend <= end
  )

prophage_bgc <- 
  select(prophage_bgc, 
         genome.x, type, product, MGE_type)

prophage_bgc <- 
  left_join(prophage_bgc, 
            select(myxo_genome_metadata, Genome, class),
            by = c("genome.x"="Genome"))
###is上的bgc####

IS_bgc <- bgc_systems %>%
  inner_join(IS_n, 
             by = c("contig" = "contig")) %>%
  filter(
    bgcbegin >= start,
    bgcend <= end
  )

IS_bgc <- 
  select(IS_bgc, genome.x, type, product, MGE_type)

IS_bgc <- 
  left_join(IS_bgc, 
            select(myxo_genome_metadata,
                   Genome, class),
            by = c("genome.x"="Genome"))
###composite_transposons上的bgc####

composite_transposons_bgc <- 
  bgc_systems %>%
  inner_join(composite_transposons_n, 
             by = c("contig" = "contig")) %>%
  filter(
    bgcbegin >= start,
    bgcend <= end
  )

composite_transposons_bgc <- 
  select(composite_transposons_bgc, genome.x, 
         type, product, MGE_type)

composite_transposons_bgc <- 
  left_join(composite_transposons_bgc, 
            select(myxo_genome_metadata, Genome, class),
            by = c("genome.x"="Genome"))



###Integron上的bgc####

integron_bgc <- bgc_systems %>%
  inner_join(integron_n, 
             by = c("contig" = "contig")) %>%
  filter(
    bgcbegin >= start,
    bgcend <= end
  )

integron_bgc <- 
  select(integron_txss, genome.x, 
         type, product, MGE_type)

integron_bgc <- left_join(integron_bgc, 
                           select(myxo_genome_metadata, Genome, class),
                           by = c("genome.x"="Genome"))
###ICE上的bgc####

ice_bgc <- bgc_systems %>%
  inner_join(integron_n, 
             by = c("contig" = "contig")) %>%
  filter(
    bgcbegin >= start,
    bgcend <= end
  )

ice_bgc <- select(ice_bgc, genome.x, type, product,
                  MGE_type)

ice_bgc <- 
  left_join(ice_bgc, 
            select(myxo_genome_metadata, Genome, class),
            by = c("genome.x"="Genome"))

###MGE_bgc 归一化的点热图####
# all_mge_bgc <- 
#   rbind(plasmid_bgc, prophage_bgc, 
#         composite_transposons_bgc, ice_bgc, 
#         integron_bgc)
# all_mge_bgc <- 
#   all_mge_bgc[,-1]
# 
# 
# # 1. 计算计数
# mge_bgc_dot_heatmap <- all_mge_bgc %>%
#   count(type, MGE_type, class, name = "count")
# 
# # mge_ds_dot_heatmap$class <- 
# #   factor(dot_heatmap_ds$Tax, levels = c("Ef","Sa","Kp","Ab","Pa","En"))
# mge_txss_dot_heatmap <- 
#   mge_txss_dot_heatmap[mge_txss_dot_heatmap$class %in% top10_classes,]
# 
# 
# mge_txss_type_totals <- mge_txss_dot_heatmap %>%
#   group_by(type) %>%
#   summarise(total_count = sum(count, na.rm = TRUE))
# mge_txss_type_totals <- 
#   mge_txss_type_totals[mge_txss_type_totals$total_count >=1,]
# 
# mge_txss_dot_heatmap <- 
#   mge_txss_dot_heatmap[mge_txss_dot_heatmap$type %in% mge_txss_type_totals$type,]
# 
# 
# # 添加归一化列
# mge_txss_dot_heatmap_normalized <- 
#   mge_txss_dot_heatmap %>%
#   left_join(tax_summary, by = c("MGE_type", "class")) %>%
#   mutate(
#     count_normalized = count / total_count,
#     percent = (count / total_count) * 100  # 可选：转换为百分比
#   ) %>%
#   arrange(type, MGE_type)
# mge_txss_dot_heatmap_normalized$MGE_type[mge_txss_dot_heatmap_normalized$MGE_type == "Composite_transposon"] <- "IS_Tn"
# 
# # 3. 创建分面热图
# p_dot_heatmap_mge_txss_normalized <-
#   ggplot(mge_txss_dot_heatmap_normalized, aes(x = class, y = type)) +
#   geom_point(aes(size = count_normalized, color = MGE_type), alpha = 1) +
#   scale_size_continuous(
#     name = "Count",
#     range = c(1, 10),  # 控制点的大小范围
#     breaks = pretty_breaks(n = 5)
#   ) +
#   scale_color_manual(values = c("IS"="#9f86c0","Plasmid"="#FFBE7D",
#                                 "Integron"="#76c893","ICEs"="#F0E442",
#                                 "IS_Tn"="#ff99c8",
#                                 "Prophage"="#56B4E9")
#   )+
#   facet_grid(~ MGE_type, scales = "free", space = "free") +
#   #facet_wrap(~ MGE_type, scales = "free", ncol = 5) +
#   theme(
#     panel.spacing = unit(0.1, "lines"),  # 将分面间距设为0
#     panel.background = element_blank(),  # 设置背景为透明
#     # panel.grid.major = element_blank(),  # 去掉主要网格线
#     panel.grid.major = element_line(color = "gray90", linewidth = 0.2),
#     panel.grid.minor = element_blank(),  # 去掉次要网格线
#     panel.border = element_rect(color = "grey20", fill = NA, linewidth = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
#     axis.title = element_blank(),
#     #axis.title = element_text( size = 12, color = "black"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
#     axis.text = element_text(size = 10, color = "black"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
#     axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),  # 横轴标签旋转
#     axis.ticks = element_line(color = "grey20", linewidth = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
#     legend.position = 'none',
#     #legend.position = c(0.2,0.75),  # 设置图例位置在右侧
#     legend.title = element_blank(), #设置图例标题字体为粗体，大小为12，颜色为深灰色
#     legend.text = element_blank(), #设置图例文本字体大小为10，颜色为深灰色
#     legend.background = element_blank(),  # 设置图例背景为透明
#     legend.box.background = element_rect(color = "grey80", linewidth = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
#     legend.key = element_blank()  # 去掉图例键的背景
#   )
# p_dot_heatmap_mge_txss_normalized
# 
# p_dot_heatmap_mge_txss_normalized_legend <-
#   ggplot(mge_txss_dot_heatmap_normalized, aes(x = class, y = type)) +
#   geom_point(aes(size = count_normalized, color = MGE_type), alpha = 1) +
#   scale_size_continuous(
#     name = "Count",
#     range = c(1, 10),  # 控制点的大小范围
#     breaks = pretty_breaks(n = 5)
#   ) +
#   scale_color_manual(values = c("IS"="#9f86c0","Plasmid"="#FFBE7D",
#                                 "Integron"="#76c893","ICEs"="#F0E442",
#                                 "Composite_transposon"="#ff99c8",
#                                 "Prophage"="#56B4E9")
#   )+
#   facet_grid(~ MGE_type, scales = "free", space = "free") +
#   #facet_wrap(~ MGE_type, scales = "free", ncol = 5) +
#   theme(
#     panel.spacing = unit(0.1, "lines"),  # 将分面间距设为0
#     panel.background = element_blank(),  # 设置背景为透明
#     # panel.grid.major = element_blank(),  # 去掉主要网格线
#     panel.grid.major = element_line(color = "gray90", linewidth = 0.2),
#     panel.grid.minor = element_blank(),  # 去掉次要网格线
#     panel.border = element_rect(color = "grey20", fill = NA, linewidth = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
#     axis.title = element_blank(),
#     #axis.title = element_text( size = 12, color = "black"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
#     axis.text = element_text(size = 10, color = "black"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
#     axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),  # 横轴标签旋转
#     axis.ticks = element_line(color = "grey20", linewidth = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
#     #legend.position = 'none',
#     legend.position = c(0.2,0.75),  # 设置图例位置在右侧
#     #legend.title = element_blank(), #设置图例标题字体为粗体，大小为12，颜色为深灰色
#     #legend.text = element_blank(), #设置图例文本字体大小为10，颜色为深灰色
#     legend.background = element_blank(),  # 设置图例背景为透明
#     legend.box.background = element_rect(color = "grey80", linewidth = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
#     legend.key = element_blank()  # 去掉图例键的背景
#   )
# p_dot_heatmap_mge_txss_normalized_legend
# 
# ggsave(plot = p_dot_heatmap_mge_ds_normalized_legend,
#        '1_p_dot_heatmap_mge_ds_normalized_legend.pdf',
#        width = 6.9,height = 5.5,units = 'in')





#组图####
library(Cairo)
library(cowplot)

g1 <- plot_grid(plot_txss_num,plot_AMP_num,
              #  labels = LETTERS[1:2], 
                ncol = 1,
              #  rel_widths = c(1, 1),
              rel_heights = c(2,1),
                align = "v",        # 垂直方向对齐绘图区域
                axis = "lr"         # 对齐左右两侧的边界（即 Y 轴）
                # label_x = 0.09,            # 标签x坐标（1为最右侧）
                #  label_y = 0.9             # 标签y坐标（1为最上方）
)+theme(plot.margin = margin(0, 0, 0, 0, "pt"))
g1

g1_2 <- plot_grid(plot_DS_num,g1,plot_bgc_num,
               # labels = LETTERS[1:2], 
                ncol = 3,
                rel_widths = c(1, 1,1),
           #     align = "v",        # 垂直方向对齐绘图区域
                axis = "lr"         # 对齐左右两侧的边界（即 Y 轴）
                # label_x = 0.09,            # 标签x坐标（1为最右侧）
                #  label_y = 0.9             # 标签y坐标（1为最上方）
)
g1_2



plot_DS_prevalence_heatmap

plot_bgc_num
plot_bgc_prevalence_heatmap
g2 <- 
  plot_grid(
            plot_txss_prevalence_heatmap,
            plot_AMP_prevalence_heatmap,
            plot_bgc_prevalence_heatmap,
            # labels = LETTERS[1:2], 
            ncol = 3,
            rel_widths = c(2, 1,2),
          #  align = "v",        # 垂直方向对齐绘图区域
          #  axis = "lr" ,        # 对齐左右两侧的边界（即 Y 轴）
            align = "h",
            axis = "tb"
            # label_x = 0.09,            # 标签x坐标（1为最右侧）
            #  label_y = 0.9             # 标签y坐标（1为最上方）
  )
g2

g2_1 <- 
  plot_grid(plot_DS_prevalence_heatmap,
            p_mge_num_heatmap, 
            # labels = LETTERS[1:2], 
            ncol = 2,
              rel_widths = c(3, 1.3),
           # rel_heights = c(10, 1),
            #  align = "v",        # 垂直方向对齐绘图区域
            #  axis = "lr" ,        # 对齐左右两侧的边界（即 Y 轴）
               align = "vh",
            axis = "tb"
            # label_x = 0.09,            # 标签x坐标（1为最右侧）
            #  label_y = 0.9             # 标签y坐标（1为最上方）
  )
g2_1


g2_2 <- 
  plot_grid(g2,
            g2_1,
            # labels = LETTERS[1:2], 
            ncol = 1,
          #  rel_widths = c(1, 1),
            rel_heights = c(1, 1.2),
            #  align = "v",        # 垂直方向对齐绘图区域
            #  axis = "lr" ,        # 对齐左右两侧的边界（即 Y 轴）
         #   align = "h",
            axis = "tb"
            # label_x = 0.09,            # 标签x坐标（1为最右侧）
            #  label_y = 0.9             # 标签y坐标（1为最上方）
  )
g2_2

g3 <- 
  plot_grid(p_dot_heatmap_mge_ds_normalized,
            p_dot_heatmap_mge_txss_normalized,
            p_dot_heatmap_mge_amp_normalized,
            # labels = LETTERS[1:2], 
            ncol = 3,
            #  rel_widths = c(1, 1),
            rel_heights = c(1, 1.2),
            #  align = "v",        # 垂直方向对齐绘图区域
            #  axis = "lr" ,        # 对齐左右两侧的边界（即 Y 轴）
               align = "vh",
            axis = "tb"
            # label_x = 0.09,            # 标签x坐标（1为最右侧）
            #  label_y = 0.9             # 标签y坐标（1为最上方）
  )
g3


ggsave(plot = g1_2,
       '1_g1_2_function_num.pdf',
       width = 10,height = 3,
       units = 'in')
ggsave(plot = g2_2,
       '1_g2_2_function_prevalence.pdf',
       width = 10,height = 6,
       units = 'in')
ggsave(plot = g3,
       '1_g3_MGE_function.pdf',
       width = 10,height = 3,
       units = 'in')



#MGE携带某类基因####

##MGE携带DS####
##composite_transposons##
composite_transposons_DS <- 
  DS %>%
  inner_join(composite_transposons_n, 
             by = c("contig" = "contig")) %>%
  filter(
    DSbegin >= start,
    DSend <= end
  )

# 手动选择contig，
# GCA_030426145_1-148

##plasmid##
plasmid_DS <- DS %>%
  inner_join(plasmid_n, 
             by = c("contig" = "contig")) %>%
  filter(
    DSbegin >= start,
    DSend <= end
  )
##手动选择contig，
##GCA_903838005_1-122 （Cas cluster）

##prophage##
prophage_DS <- DS %>%
  inner_join(Prophage_n, 
             by = c("contig" = "contig")) %>%
  filter(
    DSbegin >= start,
    DSend <= end
  )
##手动选择contig，
## GCA_964220685_1-285 (5个DS)
prophage_plot_use <- genomad_prophage[genomad_prophage$contig=="GCA_964220685_1-285",]



##MGE携带分泌系统####
##plasmid##
plasmid_txss <- txss_systems %>%
  inner_join(plasmid_n, 
             by = c("contig" = "contig")) %>%
  filter(
    SSbegin >= start,
    SSend <= end
  )
## 手动选择contig，
## GCA_023423535_1-96 （2种分泌系统）
GCA_023423535_1_96_ssg <- txss[txss$contig=="GCA_023423535_1-96",]


##MGE携带抗菌肽####
##integron##
integron_amp <- AMP %>%
  inner_join(integron_n, 
             by = c("contig" = "contig")) %>%
  filter(
    ampbegin >= start,
    ampend <= end
  )
## 手动选择contig，
## GCA_035648655_1-436
integron_plot_use <- integrons[integrons$ID_replicon == "GCA_035648655_1-436",]

##plasmid##
plasmid_amp <- 
  AMP %>%
  inner_join(plasmid_n, 
             by = c("contig" = "contig")) %>%
  filter(
    ampbegin >= start,
    ampend <= end
  )
## 手动选择contig，
## GCA_029880325_1-90

##prophage##
prophage_amp <- AMP %>%
  inner_join(Prophage_n, 
             by = c("contig" = "contig")) %>%
  filter(
    ampbegin >= start,
    ampend <= end
  )
## 手动选择contig，
##GCF_006402415_1-1 （这个contig巨长）
prophage_plot_use <- 
  genomad_prophage[genomad_prophage$contig=="GCF_006402415_1-1",]



##MGE携带基因注释####
selected_mge_gff_info <- read.delim2("all_mge_gff.tsv")
selected_mge_gff_info$type <- NA_character_

# 创建从 hit_id 到 hit_gene_ref 的查找表
lookup <- setNames(txss$hit_gene_ref, txss$hit_id)
# 用 match 填充 selected_mge_gff_info 的 name 列
selected_mge_gff_info$name <- lookup[selected_mge_gff_info$hit_id]
lookup <- setNames(txss$type, txss$hit_id)
selected_mge_gff_info$type <- lookup[selected_mge_gff_info$hit_id]


# 第二次填充（来自 DF_gene_info），只填充仍为 NA 的 name
na_idx <- is.na(selected_mge_gff_info$name)
lookup_gene <- setNames(DF_gene_info$gene_name, DF_gene_info$hit_id)
selected_mge_gff_info$name[na_idx] <- 
  lookup_gene[selected_mge_gff_info$hit_id[na_idx]]
lookup_gene <- setNames(DF_gene_info$type, DF_gene_info$hit_id)
selected_mge_gff_info$type[na_idx] <- 
  lookup_gene[selected_mge_gff_info$hit_id[na_idx]]


na_idx <- is.na(selected_mge_gff_info$name)
lookup_gene <- setNames(AMP$AMP_family, AMP$Access)
selected_mge_gff_info$name[na_idx] <- 
  lookup_gene[selected_mge_gff_info$hit_id[na_idx]]
selected_mge_gff_info$type[na_idx] <-   
  lookup_gene[selected_mge_gff_info$hit_id[na_idx]]

write.table(selected_mge_gff_info,"selected_mge_gff_info.tsv",
            sep = "\t", quote = F, row.names = F)

IS_plot_use <- IS[IS$contig == "GCA_030426145_1-148",]



#MGE之间的重叠携带某类基因####







save.image("./1_myxo_functional_anno.RData")
load("./1_myxo_functional_anno.RData")
