
setwd("D:/02_pku/01_research/06_predatory_bacteria_virome/03_phage_gene_anno")

all_phage_gff <- read.delim2("selected_phage_treated.gff",
                                header = F, sep = '\t')

library(dplyr)
library(stringr)
all_phage_gff <- all_phage_gff %>%
  mutate(
    # 提取ID=后面的内容，然后取最后一个_后的数字
    gene_num = str_extract(V9, "ID=[^;]+") %>% 
      str_remove("ID=") %>% 
      str_extract("[^_]+$"),
    hit_id = paste0(V1, "_", gene_num)
  ) %>%
  select(-gene_num)  # 移除临时列
all_phage_gff <- all_phage_gff[,c(1,3,4,5,7,10)]

colnames(all_phage_gff) <- c("contig","type","start","stop","strand",
                           "name")

#type和name可以改为需要的内容（type决定颜色，name决定标签）

write.table(all_phage_gff,"all_phage_gff.tsv",sep = "\t",quote = F,row.names = F)


primary_anno <- read.delim2("phage_annotation_summary.tsv")

all_phage_gff$seq <- sprintf("seq%d", 1:nrow(all_phage_gff))


all_phage_gff <- all_phage_gff %>%
  left_join(primary_anno, by = "seq")

write.table(all_phage_gff,"all_phage_gff_anno.tsv",sep = "\t",quote = F,row.names = F)
save.image("./1_phage_gene_annotation.RData")
################################################################################

all_phage_gff_manually <- read.delim2("all_phage_gff_anno.tsv")

library(dplyr)
library(stringr)

classified_gff <- all_phage_gff_manually %>%
  mutate(
    # 保留原始type列信息
    original_type = type,
    
    # 开始根据name列的内容创建新的type列
    type = case_when(
      # ========== 1. Lysis-related proteins ==========
      # 裂解相关：溶菌酶、孔蛋白、内溶素、跨膜蛋白等
      str_detect(name, regex("lysozyme|endolysin|holin|lysis|spanin|transglycosylase|amidase|L-alanyl-D-glutamate peptidase|peptidase|murein", ignore_case = TRUE)) ~ "Lysis",
      
      # ========== 2. Virion structure and morphogenesis ==========
      # 病毒颗粒结构和形态发生：衣壳、尾部、纤维、门户、组装蛋白等
      str_detect(name, regex("capsid|tail|fiber|portal|terminase|virion|morphogenesis|head|sheath|baseplate|ATPase|tape|neck|assembly|chaperone|adaptor|completion|decoration|structural protein|spike|wedge|hub|tube|major tail|minor tail|major head|minor head|tail terminator|head-tail adaptor|head maturation protease|head closure|baseplate|tail length tape|tail assembly|head decoration|scaffolding", ignore_case = TRUE)) ~ "Virion structure and morphogenesis",
      
      # ========== 3. DNA replication, recombination, repair and packaging ==========
      # DNA复制、重组、修复和包装
      str_detect(name, regex("replication|replicase|helicase|primase|polymerase|ligase|recombination|recombinase|Rec|integrase|excisionase|exonuclease|nuclease|endonuclease|terminase|packaging|Dna|parition|Par|partition|segregation|condensation|Scp|SMC|RusA|Ruv|Holliday|resolution|resolvase|ssb|single.*strand.*DNA|DNA.*binding|DNA.*methyl|methyltransferase|DNA.*repair|Rad|Sbc|nucleotidyltransferase|initiation|recombination|antirecombination", ignore_case = TRUE)) ~ "DNA replication, recombination, repair and packaging",
      
      # ========== 4. Host interaction and anti-defense systems ==========
      # 宿主互作和抗防御系统
      str_detect(name, regex("anti.*restriction|antirestriction|anti.*defense|anti.*CRISPR|DarB|Acr|immunity|superinfection|toxin|antitoxin|addiction|antibiotic|resistance|virulence|vir|Vfr|regulator|transcriptional|repressor|activator|ArsR|LysR|TetR|LacI|MerR|AraC|sigma|anti-termination|Q|CII|CIII|regulatory|FleQ|ComR|chromatin|remodeler|dps|oxidative|stress|HigB|toxin|antidote", ignore_case = TRUE)) ~ "Host interaction and anti-defense systems",
      
      # ========== 5. Metabolism and biosynthesis ==========
      # 代谢和生物合成
      str_detect(name, regex("phosphofructokinase|pyruvate.*dehydrogenase|kinase|dehydrogenase|transferase|reductase|synthase|biosynthesis|metabolism|metabolic|glycosyltransferase|methyltransferase|acetyltransferase|phosphatase|phosphorylase|peptidase|dipeptidase|dipeptidyl.*peptidase|protease|phosphoadenosine", ignore_case = TRUE)) ~ "Metabolism and biosynthesis",
      
      # ========== 6. Membrane and transport proteins ==========
      # 膜蛋白和转运蛋白
      str_detect(name, regex("membrane.*protein|periplasmic|lipoprotein|transport|transporter|secretion|Tat|YifL|channel|pore|porin|efflux", ignore_case = TRUE)) ~ "Membrane and transport proteins",
      
      # ========== 7. DNA modification and restriction-modification ==========
      # DNA修饰和限制-修饰系统
      str_detect(name, regex("restriction.*endonuclease|modification.*methylase|HhaI|RsrI|DNA.*methylase|methylase|restriction|modification", ignore_case = TRUE)) ~ "DNA modification and restriction-modification",
      
      # ========== 8. Unannotated/hypothetical proteins ==========
      # 未注释/假设蛋白
      is.na(name) | name == "" | str_detect(name, regex("^[0-9A-Z_]+$|^gp[0-9]+$|^orf[0-9]+$|hypothetical|unknown|unannotated|predicted|putative|conserved.*protein|domain.*protein|family.*protein|^protein.*$|^[A-Za-z]+[0-9]+$", ignore_case = TRUE)) ~ "Unannotated/hypothetical proteins",
      
      # ========== 9. Other and miscellaneous ==========
      # 其他和杂项
      TRUE ~ "Other function"
    ),
    
    # 添加更详细的子分类（可选）
    subtype = case_when(
      # Lysis 子分类
      str_detect(name, regex("lysozyme|endolysin|transglycosylase|amidase|peptidase", ignore_case = TRUE)) & type == "Lysis" ~ "Endolysin",
      str_detect(name, regex("holin|spanin", ignore_case = TRUE)) & type == "Lysis" ~ "Membrane disruption",
      
      # Virion 子分类
      str_detect(name, regex("capsid|head|major capsid", ignore_case = TRUE)) & type == "Virion structure and morphogenesis" ~ "Capsid/Head",
      str_detect(name, regex("tail|sheath|tube|fiber", ignore_case = TRUE)) & type == "Virion structure and morphogenesis" ~ "Tail",
      str_detect(name, regex("portal|terminase", ignore_case = TRUE)) & type == "Virion structure and morphogenesis" ~ "DNA packaging",
      str_detect(name, regex("assembly|chaperone|adaptor|completion", ignore_case = TRUE)) & type == "Virion structure and morphogenesis" ~ "Assembly/Chaperone",
      
      # DNA 子分类
      str_detect(name, regex("replication|primase|helicase|polymerase", ignore_case = TRUE)) & type == "DNA replication, recombination, repair and packaging" ~ "Replication",
      str_detect(name, regex("recombination|recombinase|integrase|excisionase", ignore_case = TRUE)) & type == "DNA replication, recombination, repair and packaging" ~ "Recombination/integration",
      str_detect(name, regex("repair|nuclease|exonuclease|endonuclease|ligase", ignore_case = TRUE)) & type == "DNA replication, recombination, repair and packaging" ~ "Repair",
      str_detect(name, regex("partition|Par|segregation", ignore_case = TRUE)) & type == "DNA replication, recombination, repair and packaging" ~ "Partition/segregation",
      
      # 其他情况的默认值
      TRUE ~ NA_character_
    )
  )

# 查看分类结果统计
category_summary <- classified_gff %>%
  count(type, sort = TRUE)

print("基因功能分类统计:")
print(category_summary)

# 查看每个类别中的示例基因
examples_by_category <- classified_gff %>%
  group_by(type) %>%
  slice_head(n = 5) %>%
  select(name, type, subtype) %>%
  arrange(type)

print("各分类示例基因:")
print(examples_by_category, n = Inf)

# 验证一些关键基因的分类
test_genes <- c(
  "Lysozyme; peptidoglycan dydrolase, cationic peptides, monomer",
  "Major Capsid",
  "Terminase large subunit",
  "DNA methyltransferase",
  "holin; lysis",
  "Integrase",
  "anti-restriction protein",
  "transcriptional regulator",
  "phosphofructokinase",
  "membrane protein"
)

test_results <- classified_gff %>%
  filter(name %in% test_genes) %>%
  select(name, type, subtype) %>%
  distinct()

print("测试基因分类结果:")
print(test_results)


write.table(classified_gff,"selected_phage_anno_withtype.tsv",sep = "\t",
            quote = F,row.names = F)








save.image("./1_phage_gene_annotation.RData")

#protospacer位置######

# setwd("D:/pku/creative_work/1mobilome_predatory_bacteria/7phage_gene_annotation")
# 
# my_17_phage <- read.delim2("Candidate_17_new_phages.tsv")
# 
# 
# phage_genomad <- read.delim2("nr_phage_rep_seq_virus_summary.tsv")
# 
# library(dplyr)
# my_17_phage <- my_17_phage %>% 
#   left_join(phage_genomad,by=c("seq_name"="seq_name"))
# 
# 
# phage_protospacer <- read.delim2("filtered_candidate17.out",header = F)
# 
# phage_protospacer <- phage_protospacer[,c(1,2,15,16)]
# colnames(phage_protospacer) <- c("host","contig","start","stop")
# phage_protospacer$type <- "Protospacer"
# 
# library(stringr)
# phage_protospacer <- phage_protospacer %>% 
#   mutate(
#     host      = str_remove(host, "-[^-]*$"),      # 删掉最后一个 “-” 及之后
#     contig = str_remove(contig, "\\|.*")   # 删掉第一个 “|” 及之后
#   )
# 
# phage_protospacer <- unique(phage_protospacer)
# phage_protospacer$strand <- "+"
# 
# #host info
# 
# all_4649_genome_metadata <- read.delim2("all_4649_genome_metadata.tsv")
# 
# phage_protospacer <- phage_protospacer %>%
#   left_join(all_4649_genome_metadata, by = c("host"="Genome"))
# 
# phage_protospacer <- phage_protospacer[,c(2:21,1)]
# 
# write.table(phage_protospacer,"phage_protospacer_anno.tsv",quote = F,
#             sep = '\t')

#基因绘制点图####
# 古菌病毒蛋白质特征点图
# 完整的R脚本，包含包安装和错误处理
library(ggplot2)

# # 设置工作目录
# current_dir <- getwd()
# cat("当前工作目录:", current_dir, "\n")
# 
# # 1. 创建病毒名称
# viruses <- c(
#   "Helarchaeota virus Nidhogg Meg22_1012",
#   "Helarchaeota virus Nidhogg Meg22_1214", 
#   "Helarchaeota virus Ratatoskr Meg22_1012",
#   "Lokiarchaeota virus Fenrir Meg22_1012",
#   "Lokiarchaeota virus Fenrir Meg22_1214",
#   "Lokiarchaeota virus Sköll Meg22_1214"
# )
# 
# # 2. 创建蛋白质特征名称
# features <- c(
#   # 病毒结构和形态发生 (12个)
#   "Major capsid protein", "Tail fibre", "Baseplate J", "Terminase large subunit",
#   "Baseplate wedge", "Minor head protein", "Tail completion protein S", 
#   "Portal protein", "Tail sheath protein", "Aminopeptidase", 
#   "HNH endonuclease", "I7/C57 endopeptidase-like",
#   
#   # DNA复制、重组和修复 (9个)
#   "AAA ATPase repC-like", "DNA ligase with TFIIS-like Zn finger", 
#   "3' flap repair endonuclease", "ERCC4 XPF DNA repair endonuclease", 
#   "PCNA", "NucS mismatch repair endonuclease", "Archaeal-eukaryote primase", 
#   "Primase-polymerase", "DNA polymerase family B",
#   
#   # RNA和转录 (7个)
#   "RNA ligase", "Homeobox-like", "RNaseH1-like", "Ro60/TROVE2", 
#   "Elf1-like glycosyltransferase", "Elp3/MiaB/NifB", 
#   "Ribonucleoside diphosphate reductase",
#   
#   # 核苷酸代谢 (4个)
#   "Deoxynucleo(tide/side) monophosphate kinase", "Cytidine deaminase", 
#   "5'(3')-deoxyribonucleotidase", "Thioredoxin",
#   
#   # 宿主-病毒相互作用 (3个)
#   "Von willebrand factor A", "Ubiquitin-activating enzyme", 
#   "Homing endonuclease-intein"
# )
# 
# # 3. 创建数据矩阵
# cat("创建数据矩阵...\n")
# create_data <- function() {
#   data_list <- list()
#   
#   for (i in 1:length(viruses)) {
#     for (j in 1:length(features)) {
#       association <- "None"
#       
#       # 根据图像描述设置关联类型
#       if (i <= 5) {  # 前5个病毒
#         if (j <= 12) {  # 病毒结构和形态发生特征
#           if (j %in% c(1, 3, 4, 7, 8)) {
#             association <- "Caudovirales"
#           } else if (j %in% c(11, 12)) {
#             association <- "NCVOG"
#           }
#         } else if (j == 18) {  # NucS mismatch repair endonuclease
#           association <- "Other"
#         } else if (j == 21) {  # DNA polymerase family B
#           association <- "NCVOG"
#         } else if (i == 2 && j == 34) {  # 第二个病毒的Ubiquitin-activating enzyme
#           association <- "Other"
#         }
#       }
#       
#       data_list[[length(data_list) + 1]] <- data.frame(
#         virus = viruses[i],
#         feature = features[j],
#         association = association,
#         feature_num = j,
#         virus_num = i,
#         stringsAsFactors = FALSE
#       )
#     }
#   }
#   
#   return(do.call(rbind, data_list))
# }
# 
# # 创建数据
# plot_data <- create_data()
# cat("数据矩阵创建完成，共", nrow(plot_data), "行数据\n")
# 
# # 4. 设置因子顺序
# plot_data$virus <- factor(plot_data$virus, levels = rev(viruses))
# plot_data$feature <- factor(plot_data$feature, levels = features)
# plot_data$association <- factor(plot_data$association, 
#                                 levels = c("Caudovirales", "NCVOG", "Other", "None"))
# 
# # 5. 创建分类背景数据
# category_ranges <- list(
#   c(1, 12),    # 病毒结构和形态发生
#   c(13, 21),   # DNA复制、重组和修复
#   c(22, 28),   # RNA和转录
#   c(29, 32),   # 核苷酸代谢
#   c(33, 35)    # 宿主-病毒相互作用
# )
# 
# category_colors <- c("#E6E6FA", "#E0F6FF", "#E0F6FF", "#E0F6FF", "#FFE6E6")
# category_names <- c("Virion structure and morphogenesis", 
#                     "DNA replication, recombination and repair",
#                     "RNA and transcription", 
#                     "Nucleotide metabolism",
#                     "Host-virus interaction")
# 
# # 创建背景矩形数据
# background_data <- data.frame()
# for (i in 1:length(category_ranges)) {
#   range <- category_ranges[[i]]
#   background_data <- rbind(background_data, data.frame(
#     xmin = range[1] - 0.5,
#     xmax = range[2] + 0.5,
#     ymin = 0.5,
#     ymax = length(viruses) + 0.5,
#     category = category_names[i],
#     fill_color = category_colors[i],
#     stringsAsFactors = FALSE
#   ))
# }
# 
# # 6. 创建颜色映射
# association_colors <- c(
#   "Caudovirales" = "#8B4B8C",  # 紫色
#   "NCVOG" = "#20B2AA",         # 青绿色
#   "Other" = "#696969",         # 深灰色
#   "None" = "white"             # 白色
# )
# 
# # 7. 绘制图形
# cat("开始绘制图形...\n")
# p_try <- ggplot() +
#   
#   # 添加分类背景
#   geom_rect(data = background_data,
#             aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
#             fill = background_data$fill_color,
#             alpha = 0.3) +
#   
#   # 添加点
#   geom_point(data = plot_data,
#              aes(x = feature_num, y = virus_num, 
#                  color = association, fill = association),
#              size = 4, shape = 21, stroke = 0.5) +
#   
#   # 设置颜色
#   scale_color_manual(values = association_colors, name = "Association Type") +
#   scale_fill_manual(values = association_colors, name = "Association Type") +
#   
#   # 设置坐标轴
#   scale_x_continuous(breaks = 1:length(features), labels = features) +
#   scale_y_continuous(breaks = 1:length(viruses), labels = rev(viruses)) +
#   
#   # 设置主题
#   theme_minimal() +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
#     axis.text.y = element_text(size = 9),
#     axis.title = element_blank(),
#     panel.grid = element_blank(),
#     legend.position = "bottom",
#     plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
#     plot.subtitle = element_text(hjust = 0.5, size = 10)
#   ) +
#   
#   # 添加标题
#   labs(
#     title = "Archaeal Virus Protein Feature Association",
#     subtitle = "Presence and association type of proteins across six archaeal viruses"
#   )
# 
# # 8. 添加分类标签
# for (i in 1:length(category_ranges)) {
#   range <- category_ranges[[i]]
#   x_pos <- mean(range)
#   p_try <- p_try + 
#     annotate("text", x = x_pos, y = length(viruses) + 0.8, 
#              label = category_names[i], angle = 45, hjust = 0, 
#              size = 3, fontface = "bold")
# }
# p_try
#MY-phage###############

my_phage <- read.delim2("selected_phage_anno_withtype.tsv")
phage_6 <- readxl::read_xlsx("novel_6_complete_phage_info.xlsx")

my_phage <- my_phage[my_phage$contig %in% phage_6$ID,]
my_phage <- my_phage %>%
  left_join(phage_6 %>% select(ID, proposed_name), 
            by = c("contig" = "ID")) %>%
  rename(virus_name = proposed_name)

my_phage <- my_phage[,c(2,8,13)]

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
MIN_VIRUS_COUNT <- 1  # 蛋白质至少需要在多少个病毒中出现
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
plot_data$virus <- factor(plot_data$virus, levels = c("MV-1","MV-2",
                                                      "MV-3","MV-4",
                                                      "MV-5","MV-6"))
#plot_data$protein <- factor(plot_data$protein, levels = proteins_filtered)

# # 根据type分类和首字母对蛋白质进行排序
# protein_order <- plot_data %>%
#   select(protein, type) %>%
#   distinct() %>%
#   arrange(type, protein) %>%
#   pull(protein)

# 方案1：添加 unique()
protein_order <- plot_data %>%
  select(protein, type) %>%
  distinct() %>%
  arrange(type, protein) %>%
  pull(protein) %>%
  unique()
#protein_order <- protein_order[1:50]
# # 方案2：直接对 protein 去重
# protein_order <- plot_data %>%
#   distinct(protein, type) %>%
#   arrange(type, protein) %>%
#   pull(protein) %>%
#   unique()
plot_data$protein <- factor(plot_data$protein, levels = protein_order)

# 5. 创建颜色映射
# n_colors <- length(categories)
# if (n_colors <= 8) {
#   colors <- RColorBrewer::brewer.pal(max(3, n_colors), "Set3")[1:n_colors]
# } else {
#   colors <- rainbow(n_colors)
# }

colors <- 
  c("#D9D9D9","#FB492F","#01CABF","#3C5BFE","#E5049C","#8020F3","#D9D9D9","#0FF67A","#A98467")
names(colors) <- categories
# 将NA值转换为字符串"NA"
plot_data$type[is.na(plot_data$type)] <- "NA"
# 为NA值添加白色
colors <- c(colors, "white")
names(colors)[length(colors)] <- "NA"

# 6. 绘制图形
cat("开始绘制图形...\n")

#plot_data <- plot_data[!plot_data$type %in% c("Unannotated/hypothetical proteins","NA"),]

p_gene_dot <- ggplot(plot_data, aes(x = virus, y = protein)) +
  
  # 添加点
  geom_point(aes( fill = type), color = "black",
             size = 4, shape = 21, stroke = 0.2) +
  
  # 设置颜色
  #  scale_color_manual(values = colors, name = "Protein Type") +
  scale_fill_manual(values = colors, name = "Protein Type",na.value = "white") +
  
  # 设置主题
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10),
    axis.text.y = element_text(size = 8),
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
print(p_gene_dot)

ggsave("p_gene_dot.pdf",plot = p_gene_dot, units = 'in', width = 5.5, height = 12)

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








save.image("./1_phage_gene_annotation.RData")
load("./1_phage_gene_annotation.RData")
