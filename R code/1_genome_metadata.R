# read data
# setwd("D:/pku/creative_work/1mobilome_predatory_bacteria/1genome_metadata")
# Genbank_genome <- read.delim2("gtdb_r226_predators_metadata_new_2.tsv")
# CRBC_genome <- read.delim2("CRBC_merged_results.tsv")
# GEM_genome <- read.delim2("GEM_merged_results.tsv")
# GOMC_genome <- read.delim2("GOMC_merged_results.tsv")
# GROWdb_genome <- read.delim2("GROWdb_merged_results.tsv")
# SMAG_genome <- read.delim2("SMAG_merged_results.tsv")
# TG2G_genome <- read.delim2("TG2G_merged_results.tsv")
# library(tidyr)
# library(dplyr)
# Genbank_genome <- Genbank_genome %>%
#   mutate(
#     source_category = case_when(
#       # 海洋环境
#       grepl("marine|seawater|ocean|hydrothermal vent|coral|sponge|sea water", 
#             ncbi_isolation_source, ignore.case = TRUE) ~ "marine",
#       
#       # 淡水环境
#       grepl("lake|river|pond|freshwater|stream|aquifer|groundwater", 
#             ncbi_isolation_source, ignore.case = TRUE) ~ "freshwater",
#       
#       # 陆地环境
#       grepl("soil|sediment|rhizosphere|root|plant|forest|grassland|sand", 
#             ncbi_isolation_source, ignore.case = TRUE) ~ "terrestrial",
#       
#       # 宿主相关
#       grepl("human|mouse|rat|pig|cow|host|gut|intestine|fecal|feces", 
#             ncbi_isolation_source, ignore.case = TRUE) ~ "host_associated",
#       
#       # 工程环境
#       grepl("bioreactor|wastewater|sludge|landfill|leachate|activated sludge|waste water", 
#             ncbi_isolation_source, ignore.case = TRUE) ~ "engineered",
#       
#       # 极端环境
#       grepl("hot spring|geothermal|alkaline|acidic|halophilic|saline", 
#             ncbi_isolation_source, ignore.case = TRUE) ~ "extreme",
#       
#       # 默认类别
#       TRUE ~ "other"
#     )
#   )
# # 查看分类分布
# category_summary <- Genbank_genome %>%
#   count(source_category) %>%
#   arrange(desc(n))
# print(category_summary)
# 
# 
# CRBC_genome$source_category <- "Rhizosphere"
# GOMC_genome$source_category <- "marine"
# GROWdb_genome$source_category <- "freshwater"
# SMAG_genome$source_category <- "soil"
# TG2G_genome$source_category <- "Glacier"
# GEM_genome$source_category <- "Earth’s microbiomes"
# 
# 
# 
# 
# 
# #1. merge####
# 
# public_data <- rbind(CRBC_genome, GEM_genome, GOMC_genome, 
#                      GROWdb_genome, SMAG_genome, TG2G_genome)
# 
# Genbank_genome <- select(Genbank_genome, 
#                          accession, checkm2_completeness, checkm2_contamination,
#                          genome_size, n50_contigs, gc_percentage, gtdb_taxonomy,
#                          source_category,ncbi_assembly_level,)
# Public_genome <- select(public_data,
#                         user_genome, Completeness, Contamination, 
#                         Genome_Size, Contig_N50, GC_Content, classification,source_category
#                         )
# Public_genome$assembly_level <- "Contig"
# 
# 
# colnames(Genbank_genome) <- c("Genome", "Completeness", 
#                               "Contamination",
#                               "genome_size", "n50_contigs", "gc_percentage", 
#                               "gtdb_taxonomy","source_category",
#                               "ncbi_assembly_level")
# colnames(Public_genome) <- c("Genome", "Completeness", 
#                               "Contamination",
#                               "genome_size", "n50_contigs", "gc_percentage", 
#                               "gtdb_taxonomy","source_category",
#                               "ncbi_assembly_level")
# Genbank_genome[, 2:6] <- lapply(Genbank_genome[, 2:6], as.numeric)
# Genbank_genome$gc_percentage <- Genbank_genome$gc_percentage / 100
# Public_genome[, 2:6] <- lapply(Public_genome[, 2:6], as.numeric)
# 
# all_4649_genome <- rbind(Genbank_genome,Public_genome)
# 
# #2. ####
# 
# library(stringr)
# 
# # 使用str_replace_all进行替换
# all_4649_genome$Genome <- str_replace_all(all_4649_genome$Genome, "\\.", "_")
# 
# all_4649_genome <- all_4649_genome %>%
#   separate(
#     gtdb_taxonomy,
#     into = c("domain", "phylum", "class", "order", "family", "genus", "species"),
#     sep = ";",
#     remove = FALSE  # 保留原始列，如果不需要可以设为TRUE
#   ) %>%
#   # 移除分类层级前缀（d__, p__, 等）
#   mutate(across(c(domain, phylum, class, order, family, genus, species),
#                 ~ str_remove(., "^[a-z]__")))
# 
# write_tsv(all_4649_genome,"all_4649_genome_metadata.tsv")
# 
setwd("D:/02_pku/01_research/06_predatory_bacteria_virome/01_genomes")
library(tidyr)
library(dplyr)

myxo_mag_metadata <- read.delim2("myxococcota_mag_metadata.tsv")

##2.1 物种组成桑基图####
# 绘制桑基图展示物种组成###################

library(tidyverse)
library(ggsankeyfier) 
library(MetBrewer)


df <- myxo_mag_metadata
df <- df[,7:13]
df[is.na(df) | df == ""] <- "Unknown"

# 统计 k 列中每个k的计数
k_counts <- df %>%
  count(domain)

# 统计 phylum 列中每个门的计数
p_counts <- df %>%
  count(phylum)
# 找出计数小于 10 的门
# low_count_p <- p_counts %>%
#   filter(n < 10) %>%
#   pull(Phylum)
# 将计数小于 1000 的门及其后续分类等级设置为 "others"
# df <- df %>%
#   mutate(across(Phylum:Genus, ~ ifelse(. == "unknown", .,  # 保留 "unknown"
#                                        ifelse(Phylum %in% low_count_p, "Others", .) # 将低计数门及后续分类设为 "Others"
#   )))

# 统计 class 列中每个门的计数
c_counts <- df %>%
  count(class)
# 找出计数小于 30 的class
low_count_c <- c_counts %>%
  filter(n < 20) %>%
  pull(class)
# 将计数小于 1000 的c及其后续分类等级设置为 "others"
df <- df %>%
  mutate(across(class:genus, ~ ifelse(. == "unknown", .,  # 保留 "unknown"
                                      ifelse(class %in% low_count_c, "Others", .) # 将低计数门及后续分类设为 "Others"
  )))

# 统计 order 列中每个门的计数
o_counts <- df %>%
  count(order)
# 找出计数小于 75 的class
low_count_o <- o_counts %>%
  filter(n < 50) %>%
  pull(order)
# 将计数小于 75 的order及其后续分类等级设置为 "others"
df <- df %>%
  mutate(across(order:genus, ~ ifelse(. == "unknown", .,  # 保留 "unknown"
                                      ifelse(order %in% low_count_o, "Others", .) # 将低计数门及后续分类设为 "Others"
  )))


# 统计 family 列中每个门的计数
f_counts <- df %>%
  count(family)
# 找出计数小于 55 的family
low_count_f <- f_counts %>%
  filter(n < 50) %>%
  pull(family)
# 将计数小于 55 的family及其后续分类等级设置为 "others"
df <- df %>%
  mutate(across(family:genus, ~ ifelse(. == "unknown", .,  # 保留 "unknown"
                                       ifelse(family %in% low_count_f, "Others", .) # 将低计数门及后续分类设为 "Others"
  )))

# 统计 genus 列中每个门的计数
g_counts <- df %>%
  count(genus)
# 找出计数小于 55 的family
low_count_g <- g_counts %>%
  filter(n < 30) %>%
  pull(genus)
# 将计数小于 55 的family及其后续分类等级设置为 "others"
df <- df %>%
  mutate(across(genus, ~ ifelse(. == "unknown", .,  # 保留 "unknown"
                                       ifelse(genus %in% low_count_g, "Others", .) # 将低计数门及后续分类设为 "Others"
  )))


library(tidyverse)
library(ggsankeyfier) 
library(MetBrewer)
# df1 <- df %>% select(2,3) %>% group_by(Domain,Phylum) %>% count() %>%
#   pivot_stages_longer(.,stages_from = c("Domain", "Phylum"),
#                       values_from = "n")
df2 <- df %>% select(3,4) %>% group_by(phylum,class) %>% count() %>%
  pivot_stages_longer(.,stages_from = c("phylum", "class"),
                      values_from = "n")
df3 <- df %>% select(4,5) %>% group_by(class,order) %>% count() %>%
  pivot_stages_longer(.,stages_from = c("class", "order"),
                      values_from = "n")
df4 <- df %>% select(5,6) %>% group_by(order,family) %>% count() %>%
  pivot_stages_longer(.,stages_from = c("order", "family"),
                      values_from = "n")
df5 <- df %>% select(6,7) %>% group_by(family,genus) %>% count() %>%
  pivot_stages_longer(.,stages_from = c("family", "genus"),
                      values_from = "n")

eukaryotes_sankey_4level <- 
  ggplot(data=df2,aes(x = stage,y =n,group = node,
                      edge_id = edge_id,connector = connector))+
  # 绘制第 1，2 层级
  geom_sankeyedge(aes(fill = node),
                  position = position_sankey(order ="ascending",v_space="auto",
                                             width = 0.01))+
  geom_sankeynode(aes(fill=node,color=node),
                  position = position_sankey(order = "ascending",v_space ="auto",
                                             width = 0.05))+
  # geom_text(data=df1 %>% filter(connector=="from"),
  #           aes(label = node),stat = "sankeynode",
  #           position = position_sankey(v_space ="auto",order="ascending",nudge_x=0.05),
  #           hjust=0,size=3,fontface="plain",color="black")+
  # 绘制第 2，3 层级
  geom_sankeyedge(data=df2,aes(fill = node),
                  position = position_sankey(order = "ascending",v_space ="auto",
                                             width = 0.01))+
  geom_sankeynode(data=df2,aes(fill=node,color=node),
                  position = position_sankey(order = "ascending",v_space ="auto",
                                             width = 0.05))+
  geom_text(data=df2 %>% filter(connector=="from"),
            aes(label = node),stat = "sankeynode",angle=0,
            position = position_sankey(v_space ="auto",order="ascending", nudge_x=0.05),
            hjust=0,size=3,vjust=0.5,color="black",fontface="plain")+
  geom_text(data=df2 %>% filter(connector=="to"),
            aes(label = node),stat = "sankeynode",angle=0,
            position = position_sankey(v_space ="auto",order="ascending",nudge_x=0.05),
            hjust=0,size=3,color="black",fontface="plain")+
  # 绘制第 3，4 层级
  geom_sankeyedge(data=df3,aes(fill = node),
                  position = position_sankey(order = "ascending",v_space ="auto",
                                             width = 0.01))+
  geom_sankeynode(data=df3,aes(fill=node,color=node),
                  position = position_sankey(order = "ascending",v_space ="auto",
                                             width = 0.05))+
  geom_text(data=df3 %>% filter(connector=="to"),
            aes(label = node),stat = "sankeynode",angle=0,
            position = position_sankey(v_space ="auto",order="ascending", nudge_x=0.05),
            hjust=0,size=3,vjust=0.5,color="black",fontface="plain")+
  # 绘制第 4，5 层级
  geom_sankeyedge(data=df4,aes(fill = node),
                  position = position_sankey(order = "ascending",v_space ="auto",
                                             width = 0.01))+
  geom_sankeynode(data=df4,aes(fill=node,color=node),
                  position = position_sankey(order = "ascending",v_space ="auto",
                                             width = 0.05))+
  geom_text(data=df4 %>% filter(connector=="to"),
            aes(label = node),stat = "sankeynode",angle=0,
            position = position_sankey(v_space ="auto",order="ascending", nudge_x=0.05),
            hjust=0,size=3,vjust=0.5,color="black",fontface="plain")+
  # 绘制第 5，6 层级
  geom_sankeyedge(data=df5,aes(fill = node),
                  position = position_sankey(order = "ascending",v_space ="auto",
                                             width = 0.01))+
  geom_sankeynode(data=df5,aes(fill=node,color=node),
                  position = position_sankey(order = "ascending",v_space ="auto",
                                             width = 0.05))+
  geom_text(data=df5 %>% filter(connector=="to"),
            aes(label = node),stat = "sankeynode",angle=0,
            position = position_sankey(v_space ="auto",order="ascending", nudge_x=0.05),
            hjust=0,size=3,vjust=0.5,color="black",fontface="plain")+
  
  coord_cartesian(clip="off")+
  scale_x_discrete(position = "bottom")+
  # scale_fill_manual(values = met.brewer("Nizami")) +
  # scale_color_manual(values = met.brewer("Nizami")) +
  theme_void()+
  theme(plot.margin = margin(0,0,0,0,unit = "in"),
        axis.text.x=element_text(color="black",face="plain",size=10,
                                 margin = margin(b=0.1,unit = "cm")),
        legend.position="none")

eukaryotes_sankey_4level


ggsave("1_MAG_tax_sankey.pdf",
       eukaryotes_sankey_4level,width = 8.5,height = 5.2,units = 'in')


##2.2 基因组大小及gc含量####
# genome_size 的复合图
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)

library(ggsignif)#geom_signif
library(gghalves)
library(ggsci)

myxo_mag_metadata <- type.convert(myxo_mag_metadata, as.is = TRUE)


# p1 <- ggplot(all_4649_genome, aes(x = phylum, y = genome_size, fill = phylum)) +
#   geom_violin(alpha = 0.6, trim = FALSE) +
#   geom_boxplot(width = 0.6, alpha = 0.8, outlier.shape = NA) +
# p1 <-   
# ggplot(myxo_mag_metadata, 
#          aes(x = class, y = genome_size)) +
#   geom_half_boxplot(aes(fill=class), color ="black",
#                     side = "l",
#                     errorbar.draw = T,
#                     outlier.shape = NA,
#                     width=0.8
#   )+
#   geom_half_point(aes(color=class),
#                   side = "r",size = 0.5,
#                   transformation_params = list(height = 0,width = 0.001,seed = 2))+
#   scale_y_continuous(labels = scales::comma) +  # 使用逗号分隔大数字
#   labs(
#        x = "Phylum",
#        y = "Genome Size (Mp)") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1),
#         legend.position = "none") +  # 移除图例，因为x轴已有标签
#   stat_compare_means(method = "kruskal.test", label = "p.format", 
#                      label.x = 1, 
#                      label.y = max(myxo_mag_metadata$genome_size, na.rm = TRUE) * 0.9)+
#   theme_minimal() +  # 使用简洁主题
#   theme(legend.position = "none")+  # 图例放在底部
#   theme(
#     panel.background = element_blank(),  # 设置背景为透明
#     panel.grid.major = element_blank(),  # 去掉主要网格线
#     panel.grid.minor = element_blank(),  # 去掉次要网格线
#     #panel.border = element_blank(),  # 去掉边界线
#     panel.border = element_rect(color = "grey20", fill = NA, size = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
#     #axis.line = element_line(color = "grey20", size = 0.5),  # 设置坐标轴线条为深灰色，更加细致
#     axis.title = element_text(face = "plain", size = 12, color = "grey20"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
#     axis.text = element_text(size = 10, color = "grey20"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
#     axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1), #横轴旋转
#     axis.ticks = element_line(color = "grey20", size = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
#     legend.position = "none",  # 设置图例位置在右侧
#     legend.title = element_text(face = "plain", size = 12, color = "grey20"),  # 设置图例标题字体为粗体，大小为12，颜色为深灰色
#     legend.text = element_text(size = 10, color = "grey20"),  # 设置图例文本字体大小为10，颜色为深灰色
#     legend.background = element_blank(),  # 设置图例背景为透明
#     legend.box.background = element_rect(color = "grey80", size = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
#     legend.key = element_blank()  # 去掉图例键的背景
#   )
# 
# print(p1)

library(ggplot2)
library(ggpubr)
library(scales)

p1 <- ggplot(myxo_mag_metadata, aes(x = class, y = genome_size)) +
  # 左侧箱线图
  geom_boxplot(aes(fill = class), color = "black",
               width = 0.6, outlier.shape = NA) +
  # 右侧散点，通过 position_nudge 向右偏移
  geom_jitter(aes(color = class), 
              width = 0.1, height = 0, size = 0.5,
              alpha = 0.7) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "Class", y = "Genome Size (Mb)") +
  stat_compare_means(method = "kruskal.test", 
                     label = "p.format",
                     label.x = 1,
                     label.y = max(myxo_mag_metadata$genome_size, na.rm = TRUE) * 0.9) +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "grey20", fill = NA, size = 0.5),
    axis.title = element_text(size = 12, color = "grey20"),
    axis.text = element_text(size = 10, color = "grey20"),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    axis.ticks = element_line(color = "grey20", size = 0.5)
  )

p1

p2 <- ggplot(myxo_mag_metadata, aes(x = class, y = gc_percentage)) +
  # 左侧箱线图
  geom_boxplot(aes(fill = class), color = "black",
               width = 0.6, outlier.shape = NA) +
  # 右侧散点，通过 position_nudge 向右偏移
  geom_jitter(aes(color = class), 
              width = 0.1, height = 0, size = 0.5,
              alpha = 0.7) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "Class", y = "GC Percentage") +
  stat_compare_means(method = "kruskal.test", 
                     label = "p.format",
                     label.x = 1,
                     label.y = max(myxo_mag_metadata$gc_percentage, na.rm = TRUE) * 0.9) +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "grey20", fill = NA, size = 0.5),
    axis.title = element_text(size = 12, color = "grey20"),
    axis.text = element_text(size = 10, color = "grey20"),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    axis.ticks = element_line(color = "grey20", size = 0.5)
  )

p2

# # gc_percentage 的复合图
# p2 <- ggplot(myxo_mag_metadata, 
#              aes(x = class, y = gc_percentage)) +
#   geom_half_boxplot(aes(fill=class), color ="black",
#                     side = "l",
#                     errorbar.draw = T,
#                     outlier.shape = NA,
#                     width=0.8
#   )+
#   geom_half_point(aes(color=class),
#                   side = "r",size = 0.5,
#                   transformation_params = list(height = 0,width = 0.001,seed = 2))+
#   labs(
#        x = "Phylum",
#        y = "GC Percentage") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1),
#         legend.position = "none") +
#   stat_compare_means(method = "kruskal.test", label = "p.format",
#                      label.x = 1, 
#                      label.y = max(myxo_mag_metadata$gc_percentage, na.rm = TRUE) * 0.9)+
#   theme_minimal() +  # 使用简洁主题
#   theme(legend.position = "none")+  # 图例放在底部
#   theme(
#     panel.background = element_blank(),  # 设置背景为透明
#     panel.grid.major = element_blank(),  # 去掉主要网格线
#     panel.grid.minor = element_blank(),  # 去掉次要网格线
#     #panel.border = element_blank(),  # 去掉边界线
#     panel.border = element_rect(color = "grey20", fill = NA, size = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
#     #axis.line = element_line(color = "grey20", size = 0.5),  # 设置坐标轴线条为深灰色，更加细致
#     axis.title = element_text(face = "plain", size = 12, color = "grey20"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
#     axis.text = element_text(size = 10, color = "grey20"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
#     axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1), #横轴旋转
#     axis.ticks = element_line(color = "grey20", size = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
#     legend.position = "none",  # 设置图例位置在右侧
#     legend.title = element_text(face = "plain", size = 12, color = "grey20"),  # 设置图例标题字体为粗体，大小为12，颜色为深灰色
#     legend.text = element_text(size = 10, color = "grey20"),  # 设置图例文本字体大小为10，颜色为深灰色
#     legend.background = element_blank(),  # 设置图例背景为透明
#     legend.box.background = element_rect(color = "grey80", size = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
#     legend.key = element_blank()  # 去掉图例键的背景
#   )
# 
# 
# print(p2)



mytheme <-   theme(
  panel.background = element_blank(),  # 设置背景为透明
  panel.grid.major = element_blank(),  # 去掉主要网格线
  panel.grid.minor = element_blank(),  # 去掉次要网格线
  panel.border = element_rect(color = "grey20", fill = NA, size = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
  axis.title = element_text( size = 12, color = "black"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
  axis.text = element_text(size = 10, color = "black"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
  axis.text.x = element_text(angle = 0, vjust = 0.8, hjust = 0.5),  # 横轴标签旋转
  axis.ticks = element_line(color = "grey20", size = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
  legend.position = 'none',
  #legend.position = c(0.2,0.75),  # 设置图例位置在右侧
  legend.title = element_text( size = 12, color = "grey20"),  # 设置图例标题字体为粗体，大小为12，颜色为深灰色
  legend.text = element_text(size = 10, color = "grey20"),  # 设置图例文本字体大小为10，颜色为深灰色
  legend.background = element_blank(),  # 设置图例背景为透明
  legend.box.background = element_rect(color = "grey80", linewidth = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
  legend.key = element_blank()  # 去掉图例键的背景
)

MAG_count <- table(myxo_mag_metadata$source_category) %>% as.data.frame()
# 在R中预览颜色
nature_colors <- c("#4E79A7", "#F28E2B", "#59A14F", "#E15759", "#79706E",
                   "#B07AA1", "#FFBE7D", "#8CD17D", "#86BCB6", "#499894", "#D37295")

# 绘制横条形图
p_compare_MAGnumber <- ggplot(MAG_count, 
                              aes(x = reorder(Var1,Freq) ,y = Freq)) +
  geom_bar(stat = "identity", aes(fill = Var1),width = 0.5) +  # 使用stat = "identity"因为我们的值已经是计数或比例
  coord_flip() +  # 翻转坐标轴，使条形图水平
  # scale_fill_gradient(low = "#FFF1EB",
  #                     high = "#ACE0F9")+
  # scale_fill_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF",
  #                              "LCJ"="#FFE2CE","s"= "#C4D8E9",
  #                              "CJ"= "#BEBCDF","HH"="#FB8C62","Ⅶ"="#FFF3CA"
  # )) +  # 手动设置颜色
  scale_fill_manual(values = nature_colors)+
  # scale_fill_viridis_c(option = "inferno") + ggtitle("'inferno'")+
#  geom_vline(xintercept = 4.5, linetype = "dashed", color = "black")+
#  geom_vline(xintercept = 5.5, linetype = "dashed", color = "black")+
  geom_text(aes(label = round(Freq,2)),vjust = 0.5,hjust = 0,size = 4)+
  labs(x = "", y = "Number of MAGs")+  # 设置标签和标题
  mytheme
p_compare_MAGnumber


library(Cairo)
library(cowplot)
#cairo_pdf(filename = "1_genome_metadata.pdf" , width = 8,height = 4,unit = 'in')
g1 <- plot_grid(p1,p2,
               labels = letters[2:3], ncol = 1
               #rel_heights = c(1, 1, 1.5),
               # label_x = 0.09,            # 标签x坐标（1为最右侧）
               #  label_y = 0.9             # 标签y坐标（1为最上方）
)
g1

g2 <- plot_grid(p_compare_MAGnumber,g1,
                labels = letters[3:4], ncol = 2,
                rel_widths = c(1, 1.2)
                #rel_heights = c(1, 1, 1.5),
                # label_x = 0.09,            # 标签x坐标（1为最右侧）
                #  label_y = 0.9             # 标签y坐标（1为最上方）
)
g2

ggsave("1_genome_metadata.pdf",g2,width = 8.5,height = 5,units = 'in')


save.image("./1_genome_metadata.RData")



