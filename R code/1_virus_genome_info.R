
setwd("D:/02_pku/01_research/06_predatory_bacteria_virome/01_phage_genome_info")
library(tidyr)
library(dplyr)
library(stringr)

#筛选出属于黏细菌的MAG####
all_mag <- read.delim2("all_4649_genome_metadata.tsv")
myxococcota_mag <- all_mag[all_mag$phylum %in% c("Myxococcota","Myxococcota_A"),]
write.table(myxococcota_mag,"myxococcota_mag_metadata.tsv",quote = F, sep = "\t",row.names = F)





#筛选一开始的病毒序列####
IMGVR <- read.delim2("spacers_IMGVR.blastn.1mis100cov.out2",header = F)
colnames(IMGVR) <- c("spacer_name","virus_name")
IMGVR$genome_id <- sub("^(.*)-.*", "\\1", IMGVR$spacer_name)
IMGVR_myxococcota <- IMGVR[IMGVR$genome_id %in% myxococcota_mag$Genome,]
write.table(IMGVR_myxococcota,"IMGVR_myxococcota.tsv",quote = F, sep = "\t",row.names = F)

VIRE <- read.delim2("spacers_VIRE.blastn.1mis100cov.out2",header = F)
colnames(VIRE) <- c("spacer_name","virus_name")
VIRE$genome_id <- sub("^(.*)-.*", "\\1", VIRE$spacer_name)
VIRE_myxococcota <- VIRE[VIRE$genome_id %in% myxococcota_mag$Genome,]
write.table(VIRE_myxococcota,"VIRE_myxococcota.tsv",quote = F, sep = "\t",row.names = F)

genome_phage <- read.delim2("genome_derived_virus.txt",header = F)
colnames(genome_phage) <- "phage_id"
genome_phage$genome_id <- sub("\\|provirus.*", "", genome_phage$phage_id)
genome_phage$genome_id <- sub("^(.*)-.*", "\\1", genome_phage$genome_id)
genome_phage_myxococcota <- genome_phage[genome_phage$genome_id %in% myxococcota_mag$Genome,]
write.table(genome_phage_myxococcota,"genome_phage_myxococcota.tsv",quote = F, sep = "\t",row.names = F)


genome_provirus <- read.delim2("filtered_validated_proviruses.tsv")
genome_provirus$genome_id <-  sub("^(.*)-.*", "\\1", genome_provirus$source_seq)
genome_provirus_myxococcota <- genome_provirus[genome_provirus$genome_id %in% myxococcota_mag$Genome,]
write.table(genome_provirus_myxococcota,"genome_provirus_myxococcota.tsv",quote = F, sep = "\t",row.names = F)


#根据iphop结果筛选病毒结果####
myxoco_all_phage_info <- read.delim2("phage_detailed_info_trim.tsv")


iphop_res_detail_info <- read.delim2("iphop_phage2hostgenome.tsv")
length(unique(iphop_res_detail_info$Virus))
#180,宿主全是myxoco，790个vOTU全部保留


#病毒质量、物种、等的统计####
myxoco_all_phage_info <- read.delim2("phage_detailed_info_trim.tsv")
myxoco_all_phage_info$taxonomy[myxoco_all_phage_info$taxonomy == "Unclassified"] <- 
  "Unclassified;Unclassified;Unclassified;Unclassified;Unclassified;Unclassified;Unclassified"

##1. 物种桑基图 ###########################

# 拆分 classification 列
myxoco_all_phage_info_tax <- myxoco_all_phage_info %>%
  separate(
    taxonomy,
    into = c("Virus", "Realm", "Kingdom", "Phylum", "Class", "Order", "Family"),
    sep = ";",
    remove = FALSE  # 如果保留原始列，设为 FALSE；否则设为 TRUE
  )

library(tidyverse)
library(ggsankeyfier) 
library(MetBrewer)

df <- myxoco_all_phage_info_tax
df <- df[,c(27:33)]
cols <- c("Virus","Realm", "Kingdom", "Phylum", "Class", "Order", "Family")
df[df$Phylum != "Uroviricota", cols] <- "Others"

df <- df[df$Virus != "Unclassified",]


# 统计 k 列中每个k的计数
realm_counts <- df %>%
  count(Realm)

# 统计 phylum 列中每个门的计数
p_counts <- df %>%
  count(Phylum)
# 统计 class 列中每个门的计数
c_counts <- df %>%
  count(Class)

# 统计 order 列中每个门的计数
o_counts <- df %>%
  count(Order)

# 统计 family 列中每个门的计数
f_counts <- df %>%
  count(Family)

#df <- df[,c(12:18)]

library(tidyverse)
library(ggsankeyfier) 
library(MetBrewer)
df1 <- df %>% select(1,2) %>% group_by(Virus,Realm) %>% count() %>%
  pivot_stages_longer(.,stages_from = c("Virus", "Realm"),
                      values_from = "n")
df2 <- df %>% select(2,4) %>% group_by(Realm,Phylum) %>% count() %>%
  pivot_stages_longer(.,stages_from = c("Realm", "Phylum"),
                      values_from = "n")
df3 <- df %>% select(4,5) %>% group_by(Phylum,Class) %>% count() %>%
  pivot_stages_longer(.,stages_from = c("Phylum", "Class"),
                      values_from = "n")

df4 <- df %>% select(5,6) %>% group_by(Class,Order) %>% count() %>%
  pivot_stages_longer(.,stages_from = c("Class", "Order"),
                      values_from = "n")
df5 <- df %>% select(6,7) %>% group_by(Order,Family) %>% count() %>%
  pivot_stages_longer(.,stages_from = c("Order", "Family"),
                      values_from = "n")

p_phage_sankey_4level <- 
  ggplot(data=df1,aes(x = stage,y =n,group = node,
                      edge_id = edge_id,connector = connector))+
  # 绘制第 1，2 层级
  geom_sankeyedge(aes(fill = node),
                  position = position_sankey(order ="ascending",v_space="auto",
                                             width = 0.01))+
  geom_sankeynode(aes(fill=node,color=node),
                  position = position_sankey(order = "ascending",v_space ="auto",
                                             width = 0.05))+
  geom_text(data=df1 %>% filter(connector=="from"),
            aes(label = node),stat = "sankeynode",
            position = position_sankey(v_space ="auto",order="ascending",nudge_x=0.05),
            hjust=0,size=3.5,fontface="plain",color="black")+
  # 绘制第 2，3 层级
  geom_sankeyedge(data=df2,aes(fill = node),
                  position = position_sankey(order = "ascending",v_space ="auto",
                                             width = 0.01))+
  geom_sankeynode(data=df2,aes(fill=node,color=node),
                  position = position_sankey(order = "ascending",v_space ="auto",
                                             width = 0.05))+
  geom_text(data=df1 %>% filter(connector=="to"),
            aes(label = node),stat = "sankeynode",angle=0,
            position = position_sankey(v_space ="auto",order="ascending", nudge_x=0.05),
            hjust=0,size=3.5,vjust=0.5,color="black",fontface="plain")+
  geom_text(data=df2 %>% filter(connector=="to"),
            aes(label = node),stat = "sankeynode",angle=0,
            position = position_sankey(v_space ="auto",order="ascending",nudge_x=0.05),
            hjust=0,size=3.5,color="black",fontface="plain")+
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
            hjust=0,size=3.5,vjust=0.5,color="black",fontface="plain")+
  # 绘制第 4，5 层级
  # geom_sankeyedge(data=df4,aes(fill = node),
  #                 position = position_sankey(order = "ascending",v_space ="auto",
  #                                            width = 0.01))+
  # geom_sankeynode(data=df4,aes(fill=node,color=node),
  #                 position = position_sankey(order = "ascending",v_space ="auto",
  #                                            width = 0.05))+
  # geom_text(data=df4 %>% filter(connector=="to"),
  #           aes(label = node),stat = "sankeynode",angle=0,
  #           position = position_sankey(v_space ="auto",order="ascending", nudge_x=0.05),
  #           hjust=0,size=3,vjust=0.5,color="black",fontface="plain")+
  # 
  coord_cartesian(clip="off")+
  scale_x_discrete(position = "bottom")+
  # scale_fill_manual(values = met.brewer("Nizami")) +
  # scale_color_manual(values = met.brewer("Nizami")) +
  theme_void()+
  theme(plot.margin = margin(0,0,0,0,unit = "in"),
        axis.text.x=element_text(color="black",face="plain",size=10,
                                 margin = margin(b=0.1,unit = "cm")),
        legend.position="none")

p_phage_sankey_4level

##2. 每个病毒对应的宿主################################################
head(IMGVR_myxococcota)
head(VIRE_myxococcota)
head(genome_phage_myxococcota)
head(genome_provirus_myxococcota)

# 使用 match()，简单高效
# 第一步：先用 IMGVR 表赋值，未匹配的默认 NA
myxoco_all_phage_info$host_genome_id <- IMGVR_myxococcota$genome_id[
  match(myxoco_all_phage_info$ID, IMGVR_myxococcota$virus_name)
]

# 第二步：对于还是 NA 的行，去 VIRE 表中查找并填充
na_idx <- is.na(myxoco_all_phage_info$host_genome_id)
myxoco_all_phage_info$host_genome_id[na_idx] <- VIRE_myxococcota$genome_id[
  match(myxoco_all_phage_info$ID[na_idx], VIRE_myxococcota$virus_name)
]

na_idx <- is.na(myxoco_all_phage_info$host_genome_id)
myxoco_all_phage_info$host_genome_id[na_idx] <- 
  sub("^(.*)-.*", "\\1", myxoco_all_phage_info$ID[na_idx])

head(myxococcota_mag)
myxoco_all_phage_info$host_class <- myxococcota_mag$class[
  match(myxoco_all_phage_info$host_genome_id, myxococcota_mag$Genome)
]

write.table(myxoco_all_phage_info, "myxoco_all_phage_info.tsv",row.names = F,
            quote = F, sep = "\t")

plot_data <- myxoco_all_phage_info %>%
  group_by(host_class) %>%
  summarise(count = n(), .groups = 'drop') %>%
  ungroup()

plot_data[is.na(plot_data)] <- "Myxococcia"


# 绘制堆叠柱状图
p_major_phage <- 
  ggplot(plot_data, aes(x = host_class, y = count, fill = host_class)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  # 添加数值标签
  geom_text(aes(label = count), 
            position = position_stack(vjust = 1),
            size = 3, color = "black", fontface = "bold") +
  # 配色方案
  #scale_fill_brewer(palette = "Set3") +
  scale_fill_manual(values = c("Myxococcia"="#5cb9cc",
                               "B64-G9"="#8e69b7",
                               "Polyangia"="#f3be82",
                               "UBA9042"="#c2b0d1",
                               "Bradymonadia"="#d47ebe",
                               "UBA9160"="#edb9d1",
                               "UBA796"="#ef9e97")
  ) +
  # 美化主题
  theme(
    panel.background = element_blank(),  # 设置背景为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    panel.border = element_rect(color = "grey20", fill = NA, linewidth = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    axis.title = element_text( size = 12, color = "black"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 10, color = "black"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    axis.text.x = element_text(angle = 90, vjust = 0.8, hjust = 1),  # 横轴标签旋转
    axis.ticks = element_line(color = "grey20", linewidth = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    legend.position = 'none',
    #legend.position = c(0.2,0.75),  # 设置图例位置在右侧
    legend.title = element_text( size = 12, color = "grey20"),  # 设置图例标题字体为粗体，大小为12，颜色为深灰色
    legend.text = element_text(size = 10, color = "grey20"),  # 设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", linewidth = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )
p_major_phage

##3. 病毒checkv序列的质量####

library(ggplot2)
library(ggbeeswarm)  # 用于更好的散点排列
# 计算每个类别的数量
quality_counts <- myxoco_all_phage_info %>%
  count(checkv_quality, name = "count") %>%
  mutate(label = paste0(checkv_quality, " (n=", count, ")"))

quality_counts <- myxoco_all_phage_info %>%
  count(checkv_quality, name = "count") %>%
  mutate(label = paste0(checkv_quality))


# 合并计数回原始数据
myxoco_all_phage_info <- myxoco_all_phage_info %>%
  left_join(quality_counts, by = "checkv_quality")

# 确保有序因子（包含计数标签）
myxoco_all_phage_info$label <- factor(
  myxoco_all_phage_info$label,
  levels = c("Complete","High-quality",
             "Medium-quality","Low-quality"),
  ordered = TRUE
)

# # 确保checkv_quality是有序因子
# checkv_phage_quality$checkv_quality <- factor(
#   checkv_phage_quality$checkv_quality,
#   levels = c("Not-determined", "Low-quality", "Medium-quality", "High-quality", "Complete"),
#   ordered = TRUE
# )
# 创建颜色映射

# 绘制横向箱线图
p_phage_quality <- 
  ggplot(myxoco_all_phage_info, 
         aes(x = contig_length, 
             y = label, 
             fill = checkv_quality,
             color = checkv_quality)) +
  scale_y_discrete(limits = rev)+
  # 添加箱线图
  geom_boxplot(alpha = 0.7, 
               outlier.shape = 1,          # 实心圆点
               outlier.size = 1,            # 离群点大小
               outlier.alpha = 0.7,         # 离群点透明度
               width = 0.6,
               size = 0.4                   # 边框线粗细
  ) +
  # 设置颜色和填充（使用相同颜色方案）
  # scale_fill_manual(values = quality_colors) +
  # scale_color_manual(values = quality_colors) +
  # 添加散点表示所有数据点
  # geom_quasirandom(
  #   aes(color = checkv_quality),  # 使用颜色区分
  #   alpha = 0.5,                 # 半透明
  #   size = 1.5,                  # 点大小
  #   groupOnX = FALSE             # 在y轴方向排列
  # ) +
  
  # 使用viridis颜色方案
  scale_fill_viridis_d(option = "D", begin = 0.2, end = 0.8) +
  scale_color_viridis_d(option = "D", begin = 0.2, end = 0.8) +
  
  # 对数变换x轴
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  # 
  # 添加标签和主题
  labs(
    title = "Contig length distribution by CheckV quality",
    x = "Contig length (bp, log10 scale)",
    y = "CheckV quality level",
    #caption = "Each point represents one contig"
  ) +
  
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "plain"),  # 标题居中，加粗，字体大小为16
    plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),  # 设置图形边距
    plot.background = element_blank(),  # 去除图形背景
    panel.background = element_blank(),  # 去除面板背景
    panel.border = element_rect(linewidth = 0.6),  # 设置面板边框
    panel.grid = element_blank(),  # 去除网格线
    axis.title.x = element_text(size = 10, face = "plain",color = 'black'),  # 设置 x 轴标题
    axis.title.y = element_text(size = 10, face = "plain",color = 'black'),  # 设置 y 轴标题
    axis.ticks.x = element_line(linewidth = 0.5),  # 设置 x 轴刻度线
    axis.ticks.y = element_line(linewidth = 0.5),  # 设置 y 轴刻度线
    axis.ticks.length.x = unit(0.1, "cm"),  # 设置 x 轴刻度长度
    axis.ticks.length.y = unit(0.1, "cm"),  # 设置 y 轴刻度长度
    axis.text.x = element_text(angle = 0, size = 10, hjust = 0.5, vjust = 0.5,color = 'black'),  # 设置 x 轴文本
    axis.text.y = element_text(size = 10, color = 'black'),  # 设置 y 轴文本
    legend.position = "none"  # 去除图例
  )
p_phage_quality

# ggsave("1_p_phage_quality2.pdf",
#        p_phage_quality ,width = 3,height = 3,units = 'in')
# 
# ggsave("1_p_phage_sankey_4level.pdf",
#        p_phage_sankey_4level ,width = 3,height = 3,units = 'in')
# 
# ggsave("1_p_major_phage.pdf",
#        p_major_phage ,width = 3,height = 3,units = 'in')
# 
# ggsave("1_p_other_phage.pdf",
#        p_other_phage ,width = 1.5,height = 1.2,units = 'in')

library(Cairo)
library(cowplot)
#cairo_pdf(filename = "1_genome_metadata.pdf" , width = 8,height = 4,unit = 'in')
g1 <- plot_grid(p_phage_quality, p_phage_sankey_4level,p_major_phage,
                labels = LETTERS[1:3], ncol = 3,
                rel_widths = c(1, 1, 1)
                # label_x = 0.09,            # 标签x坐标（1为最右侧）
                #  label_y = 0.9             # 标签y坐标（1为最上方）
)
g1

ggsave("1_virus_overview_info.pdf",g1,width = 9.9,height = 3,units = 'in')


#ggsave("1_p_other_phage.pdf",p_other_phage,width = 2.24,height = 2.74,units = 'in')

##长度分布密度图####

library(ggplot2)
library(dplyr)

# 计算统计信息
stats <- myxoco_all_phage_info %>%
  summarise(
    mean_len = mean(length, na.rm = TRUE),
    median_len = median(length, na.rm = TRUE),
    min_len = min(length, na.rm = TRUE),
    max_len = max(length, na.rm = TRUE),
    n = n()
  )

p_phage_length_density <-
  ggplot(myxoco_all_phage_info, aes(x = length)) +
  geom_density(fill = "#4DBBD5", alpha = 0.6, color = "black", linewidth = 0.6) +
  # 添加均值和分位数线
  geom_vline(aes(xintercept = mean_len), data = stats, 
             color = "red", linewidth = 0.6, linetype = "dashed") +
  geom_vline(aes(xintercept = median_len), data = stats,
             color = "blue", linewidth = 0.6, linetype = "dotted") +
  # 添加统计信息文本
  annotate("text", x = Inf, y = Inf, 
           label = paste0("N = ", stats$n, "\n",
                          "Mean = ", round(stats$mean_len, 0), "\n",
                          "Median = ", stats$median_len),
           hjust = 1.1, vjust = 1.1, size = 4) +
  theme(
    panel.background = element_blank(),  # 设置背景为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    panel.border = element_rect(color = "grey20", fill = NA, linewidth = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    axis.title = element_text( size = 12, color = "black"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 10, color = "black"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    axis.text.x = element_text(angle = 0, vjust = 0.8, hjust = 0.5),  # 横轴标签旋转
    axis.ticks = element_line(color = "grey20", linewidth = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    legend.position = 'none',
    #legend.position = c(0.2,0.75),  # 设置图例位置在右侧
    legend.title = element_text( size = 12, color = "grey20"),  # 设置图例标题字体为粗体，大小为12，颜色为深灰色
    legend.text = element_text(size = 10, color = "grey20"),  # 设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", linewidth = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )
p_phage_length_density

ggsave("1_p_phage_length_density.pdf",p_phage_length_density,
       width = 4.85,height = 3.8,units = 'in')


#complete phage############################

complete_phage <- ESKAPE_all_phage_info_tax[ESKAPE_all_phage_info_tax$topology %in% c("DTR","ITR"),]

table(complete_phage$Host_tax)


#病毒来源的环境及其地理分布####

##IMGVR####

myxo_IMGVR_phage <- read.delim2("myxo_IMGVR_phage.tsv")
myxo_IMGVR_phage$Taxon_oid <- sub(
  "^[^_]+_[^_]+_([^_]+)_.*", 
  "\\1", 
  myxo_IMGVR_phage$ID
)

IMG_genomeCart_metadata <- readxl::read_excel("IMG.xlsx")

myxo_IMGVR_phage <- myxo_IMGVR_phage %>%
  mutate(Taxon_oid = as.numeric(Taxon_oid)) %>%    # 转为数值
  left_join(IMG_genomeCart_metadata %>% 
              select(`IMG Genome ID`, Ecosystem, `Ecosystem Category`, Latitude, Longitude),
            by = c("Taxon_oid" = "IMG Genome ID"))


table(myxo_IMGVR_phage$Ecosystem)
table(myxo_IMGVR_phage$`Ecosystem Category`)

myxo_IMGVR_phage$`Ecosystem Category` <- ifelse(
  is.na(myxo_IMGVR_phage$`Ecosystem Category`), 
  "Unknown", 
  myxo_IMGVR_phage$`Ecosystem Category`
)


category_summary <- myxo_IMGVR_phage %>%
  count(`Ecosystem Category`) %>%
  arrange(desc(n))
print(category_summary)

library(ggplot2)
library(dplyr)

nature_colors <- c("#4E79A7", "#F28E2B", "#59A14F", "#E15759", "#79706E",
                   "#B07AA1", "#FFBE7D", "#8CD17D", "#86BCB6", "#499894", "#D37295")

# 水平条形图版本
p_IMGVR_phage_ecotype <- 
  ggplot(category_summary, aes(x = n, y = reorder(`Ecosystem Category`, n), 
                               fill = `Ecosystem Category`)) +
  geom_bar(stat = "identity", width = 0.7, color = "black",linewidth = 0.1) +
  geom_text(aes(label = n), hjust = -0.2, size = 3.5, color = "black") +
  labs(
    y = "",
    x = "Number of phages identified from IMG/VR"
  ) +
  scale_fill_manual(values = c("Solid waste"="#9b5de5","Terrestrial"="#FFBE7D",
                               "Plants"="#59A14F","Aquatic"="#00b4d8","Unknown"="#79706E",
                               "Bioreactor"="#E15759","Wastewater"="#D37295","WWTP"="#D37295"
  )) +
  theme(
    panel.background = element_blank(),  # 设置背景为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    panel.border = element_rect(color = "grey20", fill = NA, linewidth = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    axis.title = element_text( size = 10, color = "black"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 10, color = "black"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    axis.text.x = element_text(angle = 0, vjust = 0.8, hjust = 0.5),  # 横轴标签旋转
    axis.ticks = element_line(color = "grey20", size = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    legend.position = 'none',
    #legend.position = c(0.2,0.75),  # 设置图例位置在右侧
    legend.title = element_text( size = 10, color = "grey20"),  # 设置图例标题字体为粗体，大小为12，颜色为深灰色
    legend.text = element_text(size = 10, color = "grey20"),  # 设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", linewidth = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )
p_IMGVR_phage_ecotype
##地理信息##########


self_selected_phage_4 <- 
  myxo_IMGVR_phage[!myxo_IMGVR_phage$Latitude %in% c("none","NA","n/a"),]
self_selected_phage_4 <- 
  self_selected_phage_4[!self_selected_phage_4$Longitude %in% c("none","NA","n/a"),]
self_selected_phage_4 <-
  self_selected_phage_4
# 方法一：直接筛选两列都不是 NA 的行
self_selected_phage_4 <- self_selected_phage_4[
  !is.na(self_selected_phage_4$Latitude) & !is.na(self_selected_phage_4$Longitude),
]

self_selected_phage_4 <- self_selected_phage_4 %>%
  mutate(
    Latitude  = as.numeric(Latitude),
    Longitude = as.numeric(Longitude)
  )


redundancy_count <- self_selected_phage_4 %>%
  select(`Ecosystem Category`, Latitude, Longitude) %>%   # 只保留最后三列
  group_by(`Ecosystem Category`, Latitude, Longitude) %>% # 按三列分组
  summarise(virus_number = n(), .groups = "drop")             # 计数

#参考https://cloud.tencent.com/developer/article/2206435
library(ggplot2)
world<-map_data("world")
world.map <- 
  ggplot() +
  geom_polygon(data=world,aes(x=long,y=lat,group=group),
               fill="#e9ecef")+
  geom_point(data=redundancy_count, 
             aes(x=Longitude,y=Latitude, 
                 fill = `Ecosystem Category`, size = virus_number),
             shape = 21,      # 带边框的填充圆
             colour = "black",# 边框颜色
             stroke = 0.05     # 边框宽度（mm）
  )+
  scale_size(range = c(2,6))+
  scale_fill_manual(values = c("Solid waste"="#9b5de5","Terrestrial"="#FFBE7D",
                               "Plants"="#59A14F","Aquatic"="#00b4d8","Unknown"="#79706E",
                               "Bioreactor"="#E15759","Wastewater"="#D37295","WWTP"="#D37295"
  )) +
  labs( x = "Longitude (°)", y = "Latitude (°)") +
  theme_bw()+
  scale_y_continuous(expand = expansion(mult=c(0,0)))+
  scale_x_continuous(expand = expansion(add=c(0,0)))+
  theme(legend.position = 'none')
world.map


world.map2 <- 
  ggplot() +
  geom_polygon(data=world,aes(x=long,y=lat,group=group),
               fill="#e9ecef")+
  geom_point(data=redundancy_count, 
             aes(x=Longitude,y=Latitude, 
                 color = `Ecosystem Category`, size = virus_number),
  )+
  scale_color_manual(values = c("Built environment"="#9b5de5","Terrestrial"="#FFBE7D",
                                "Plants"="#59A14F","Aquatic"="#00b4d8","Unknown"="#79706E",
                                "Host_associated"="#E15759","Mammals: Human"="#D37295"
  )) +
  theme_bw()+
  scale_y_continuous(expand = expansion(mult=c(0,0)))+
  scale_x_continuous(expand = expansion(add=c(0,0)))+
  theme(legend.position = 'none')
world.map2




##vire病毒####

myxo_vire_phage <- read.delim2("myxo_vire_phage.tsv")

vire_database_info <- read.delim2("vire_matched_rows.tsv")

# 加载必要的包
library(dplyr)
library(tidyr)
library(stringr)
# 创建主要类别分类
classify_ecosystem_main <- function(microntology) {
  case_when(
    # 人类相关
    str_detect(microntology, "human host") & 
      str_detect(microntology, "digestive tract|intestine") ~ "Human_Gut",
    
    str_detect(microntology, "human host") & 
      str_detect(microntology, "skin") ~ "Human_Skin",
    
    str_detect(microntology, "human host") & 
      str_detect(microntology, "airways|respiratory") ~ "Human_Respiratory",
    
    str_detect(microntology, "human host") & 
      str_detect(microntology, "urogenital") ~ "Human_Urogenital",
    
    str_detect(microntology, "human host") & 
      str_detect(microntology, "mouth|oral") ~ "Human_Oral",
    
    str_detect(microntology, "human host") & 
      !str_detect(microntology, "digestive tract|intestine|skin|airways|urogenital|mouth") ~ "Human_Other",
    
    # 其他动物相关
    str_detect(microntology, "mammalian host|animal host") & 
      !str_detect(microntology, "human") ~ "Other_Animals",
    
    # 植物相关
    str_detect(microntology, "plant host") ~ "Plants",
    
    # 水生环境
    str_detect(microntology, "aquatic:marine") ~ "Marine",
    str_detect(microntology, "aquatic:fresh water") ~ "Freshwater",
    str_detect(microntology, "aquatic:brackish") ~ "Brackish",
    str_detect(microntology, "aquatic") ~ "Aquatic_Other",
    
    # 陆地环境
    str_detect(microntology, "terrestrial:soil") ~ "Soil",
    str_detect(microntology, "terrestrial:forest") ~ "Forest",
    str_detect(microntology, "terrestrial:grassland") ~ "Grassland",
    str_detect(microntology, "terrestrial") ~ "Terrestrial_Other",
    
    # 人工环境
    str_detect(microntology, "built environment") ~ "Built_Environment",
    str_detect(microntology, "wastewater") ~ "Wastewater",
    str_detect(microntology, "agriculture") ~ "Agricultural",
    str_detect(microntology, "food") ~ "Food_Related",
    str_detect(microntology, "anthropogenic") ~ "Anthropogenic_Other",
    
    # 其他
    str_detect(microntology, "hydrothermal") ~ "Hydrothermal",
    str_detect(microntology, "contaminated") ~ "Contaminated",
    str_detect(microntology, "ancient") ~ "Ancient",
    
    TRUE ~ "Unclassified"
  )
}

# 应用主要类别分类
vire_database_info <- vire_database_info %>%
  mutate(ecosystem_main = classify_ecosystem_main(microntology))

# 查看主要类别统计
cat("\n=== 主要生态系统类别统计 ===\n")
main_summary <- table(vire_database_info$ecosystem_main, useNA = "ifany") %>%
  sort(decreasing = TRUE)
print(main_summary)

category_summary_vire <- vire_database_info %>%
  count(ecosystem_main) %>%
  arrange(desc(n))
print(category_summary_vire)

nature_colors <- c("#4E79A7", "#F28E2B", "#59A14F", "#E15759", "#79706E",
                   "#B07AA1", "#FFBE7D", "#8CD17D", "#86BCB6", "#499894", "#D37295")

# 水平条形图版本
p_vire_phage_ecotype <- 
  ggplot(category_summary_vire, aes(x = n, y = reorder(ecosystem_main, n), 
                                    fill = ecosystem_main)) +
  geom_bar(stat = "identity", width = 0.7, color = "black",linewidth = 0.1) +
  geom_text(aes(label = n), hjust = -0.2, size = 3.5, color = "black") +
  labs(
    y = "",
    x = "Number of phages identified from VIRE"
  ) +
  scale_fill_manual(values = c("Anthropogenic_Other"="#9b5de5","Terrestrial_Other"="#FFBE7D",
                               "Plants"="#59A14F","Freshwater"="#00b4d8",
                               "Aquatic_Other"="#00b4d8","Marine"="#00b4d8",
                               "Unknown"="#79706E",
                               "Other_Animals"="#E15759","Soil"="#D37295"
  )) +
  theme(
    panel.background = element_blank(),  # 设置背景为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    panel.border = element_rect(color = "grey20", fill = NA, linewidth = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    axis.title = element_text( size = 10, color = "black"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 10, color = "black"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    axis.text.x = element_text(angle = 0, vjust = 0.8, hjust = 0.5),  # 横轴标签旋转
    axis.ticks = element_line(color = "grey20", size = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    legend.position = 'none',
    #legend.position = c(0.2,0.75),  # 设置图例位置在右侧
    legend.title = element_text( size = 10, color = "grey20"),  # 设置图例标题字体为粗体，大小为12，颜色为深灰色
    legend.text = element_text(size = 10, color = "grey20"),  # 设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", linewidth = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )
p_vire_phage_ecotype

library(Cairo)
library(cowplot)
#cairo_pdf(filename = "1_genome_metadata.pdf" , width = 8,height = 4,unit = 'in')
g2 <- plot_grid(p_IMGVR_phage_ecotype,world.map,p_vire_phage_ecotype,
                labels = letters[2], ncol = 3,
                rel_widths = c(1,2, 1)
                # label_x = 0.09,            # 标签x坐标（1为最右侧）
                #  label_y = 0.9             # 标签y坐标（1为最上方）
)
g2

ggsave("1_phage_geographical_distribution.pdf",
       g2,width = 10,height = 3,units = 'in')




##vire病毒的地理信息####
# 在另一个R脚本单独处理，这里得到的结果太少，不用了








save.image("./1_virus_genome_info.RData")
load("./1_virus_genome_info.RData")
