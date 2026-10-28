# =============================================================================
# Predatory-bacterial-virome -- genome catalog / 565-genome phylogeny annotation
# =============================================================================
# iTOL annotation of the 565-genome maximum-likelihood tree: class-level
# taxonomy, method used for virus-host assignment and the CRISPR-Cas
# type/subtype rings.
#
# Manuscript      : Fig. 2A
# Main outputs    : 1_MAGs_itol_anno.RData, iTOL annotation files
# Original script : 01_genomes/1_MAGs_itol_anno.R
# -----------------------------------------------------------------------------
# Input data and output files are resolved from GTR_DATA (current analysis root),
# GTR_WORK (earlier working copy) and GTR_DEF (shared defensome data). All three
# default to the current working directory; set them instead of editing files:
#
#     Sys.setenv(GTR_DATA = "/path/to/06_predatory_bacteria_virome")
#
# =============================================================================
GTR_DATA <- normalizePath(Sys.getenv("GTR_DATA", unset = "."), mustWork = FALSE)


setwd(file.path(GTR_DATA, "01_genomes"))

library(dplyr)
MAG_565 <- read.delim2("Selected_phylo_MAGs_filtered_Myxoco.txt",header = F) %>%
  as.data.frame() %>%
  setNames("genome")

all_myxoco_mag_metadata <- read.delim2("myxococcota_mag_metadata.tsv")

MAG_565 <- MAG_565 %>% 
  left_join(all_myxoco_mag_metadata,
            by = c("genome"="Genome"))
MAG_565$CRISPRtype <- "Not determined"

all_CRISPRCas <- read.delim2("CRISPR_Cas.tab")
library(stringr)
all_CRISPRCas$genome <- str_replace(all_CRISPRCas$Contig, "-[^-]*$", "")


MAG_565 <- MAG_565 %>%
  left_join(all_CRISPRCas %>% select(genome, Prediction),
            by = c("genome" = "genome")) %>%
  mutate(CRISPRtype = Prediction) %>%
  select(-Prediction)  # 删除临时列

MAG_565$CRISPRtype[is.na(MAG_565$CRISPRtype)] <- "Not determined"

MGE_host_linkage_method_1 <- read.delim2("spacer_target_virus_MAGs.txt",header = F) %>% 
  as.data.frame() %>% setNames("genome")
MGE_host_linkage_method_1$linkage_method1 <- "CRISPR based"

MGE_host_linkage_method_2 <- read.delim2("Contain_provirus_MAGs.txt",header = F) %>% 
  as.data.frame() %>% setNames("genome")
MGE_host_linkage_method_2$linkage_method2 <- "Provirus based"

intersect(MGE_host_linkage_method_1$genome, MGE_host_linkage_method_2$genome)
#15个MAGs, 也有可能provirus的CRISPR靶向了别的phage或plasmid！

MAG_565 <- MAG_565 %>%
  left_join(MGE_host_linkage_method_1 %>% select(genome, linkage_method1),
            by = c("genome" = "genome"))

MAG_565 <- MAG_565 %>%
  left_join(MGE_host_linkage_method_2 %>% select(genome, linkage_method2),
            by = c("genome" = "genome"))

MAG_565$linkage_method1[is.na(MAG_565$linkage_method1)] <- "No"
MAG_565$linkage_method2[is.na(MAG_565$linkage_method2)] <- "No"

MAG_565_2 <- MAG_565 %>%
  mutate(
    typeI = ifelse(str_detect(CRISPRtype, "(^|[, (])I-"), "yes", "no"),
    typeII = ifelse(str_detect(CRISPRtype, "(^|[, (])II-"), "yes", "no"),
    typeIII = ifelse(str_detect(CRISPRtype, "(^|[, (])III-"), "yes", "no"),
    typeV = ifelse(str_detect(CRISPRtype, "(^|[, (])V-"), "yes", "no")
  )

write.table(MAG_565_2,"MAG_565_2.tsv",sep = "\t",row.names = F,
            quote = F)


only_myxoco <- MAG_565_2[MAG_565_2$phylum %in% c("Myxococcota","Myxococcota_A"),]



#itol注释####
library(itol.toolkit)
library(data.table) # 用于元数据文件读取
library(ape) # 用于树文件读写
library(stringr) # 用于字符串处理
library(dplyr)

# 假设您的树存储在一个名为"tree"的变量中
tree <- "bin565.unrooted.tree"
phylo <- ape::read.tree(tree)
hub <- create_hub(tree)

#按门设置颜色
library(dplyr)
df_phycolor <- MAG_565_2 %>% select(genome, class)
df_phycolor <- df_phycolor %>% distinct(genome,.keep_all = TRUE)
set.seed(666)

#手动为门添加颜色
colorss <- c(
 # "#1f77b4",  # 蓝色
 # "#ff7f0e",  # 橙色
  # "#2ca02c",  # 绿色
  # "#d62728",  # 红色
  "#9467bd",  # 紫色
 # "#8c564b",  # 棕色
  "#e377c2",  # 粉色
  # "#7f7f7f",  # 灰色
 # "#bcbd22",  # 橄榄色
  "#17becf",  # 青色
  "#aec7e8",  # 浅蓝色
  "#ffbb78",  # 浅橙色
 # "#98df8a",  # 浅绿色
  "#ff9896",  # 浅红色
  "#c5b0d5",  # 浅紫色
 # "#c49c94",  # 浅棕色
  "#f7b6d2"  # 浅粉色
 # "#dbdb8d",  # 浅橄榄色
 # "#9edae5",  # 浅青色
 # "#393b79",  # 深蓝色
 # "#637939"   # 深绿色
)
element_counts <- df_phycolor %>%
  count(class)

element_count_sorted <- element_counts[order(-element_counts$n), ]
element_count_sorted <- element_count_sorted %>%
  mutate(class = ifelse(row_number() >= 8, "Others", class)) %>%
  group_by(class) %>%
  summarise(n = sum(n)) %>%
  ungroup()

element_count_sorted$color <- c(colorss )

#
df_phycolor <- merge(df_phycolor,
                     element_count_sorted,
                     by.x = "class",by.y = "class",all.x = T)
df_phycolor <- df_phycolor[, c("genome", "class","color")]
df_phycolor$color[is.na(df_phycolor$color)] <- "#aec7e8"

#subtype = "range"才能大面积改变颜色
unit_2 <- create_unit(data = df_phycolor, 
                      key = "phylo_color_range", 
                      type = "TREE_COLORS", 
                      subtype = "range", 
                      #color = "wesanderson",  #需要的颜色超过了集合中有的颜色
                      # line_type = c(rep("normal",4),"dashed"),
                      # size_factor = 5, 
                      tree = phylo)
write_unit(unit_2)

unit_3 <- create_unit(data = df_phycolor, 
                      key = "phylo_color_branch", 
                      type = "TREE_COLORS", 
                      subtype = "branch", 
                      #color = "wesanderson",  #需要的颜色超过了集合中有的颜色
                      # line_type = c(rep("normal",4),"dashed"),
                      # size_factor = 5, 
                      tree = phylo)
write_unit(unit_3)




#typeI
Cas_typeI <- MAG_565_2 %>% select(genome, typeI)
Cas_typeI <- Cas_typeI %>%
  group_by(genome) %>%
  summarise(
    typeI = ifelse("yes" %in% typeI, "yes", "no")
  ) %>%
  ungroup()
Cas_typeI <- Cas_typeI %>%
  mutate(typeI = case_when(
    typeI == "yes" ~ 1,
    typeI == "no" ~ 0,
    TRUE ~ as.numeric(NA)  # 处理其他情况（如果有的话）
  ))


unit_4 <- create_unit(data = Cas_typeI,
                      key = 'cas_typeI',
                      type = 'DATASET_EXTERNALSHAPE',
                      tree = tree
)
write_unit(unit_4)

#typeII
Cas_typeII <- MAG_565_2 %>% select(genome, typeII)
Cas_typeII <- Cas_typeII %>%
  group_by(genome) %>%
  summarise(
    typeII = ifelse("yes" %in% typeII, "yes", "no")
  ) %>%
  ungroup()
Cas_typeII <- Cas_typeII %>%
  mutate(typeII = case_when(
    typeII == "yes" ~ 1,
    typeII == "no" ~ 0,
    TRUE ~ as.numeric(NA)  # 处理其他情况（如果有的话）
  ))


unit_4 <- create_unit(data = Cas_typeII,
                      key = 'cas_typeII',
                      type = 'DATASET_EXTERNALSHAPE',
                      tree = tree
)
write_unit(unit_4)


#typeIII
Cas_typeIII <- MAG_565_2 %>% select(genome, typeIII)
Cas_typeIII <- Cas_typeIII %>%
  group_by(genome) %>%
  summarise(
    typeIII = ifelse("yes" %in% typeIII, "yes", "no")
  ) %>%
  ungroup()
Cas_typeIII <- Cas_typeIII %>%
  mutate(typeIII = case_when(
    typeIII == "yes" ~ 1,
    typeIII == "no" ~ 0,
    TRUE ~ as.numeric(NA)  # 处理其他情况（如果有的话）
  ))


unit_4 <- create_unit(data = Cas_typeIII,
                      key = 'cas_typeIII',
                      type = 'DATASET_EXTERNALSHAPE',
                      tree = tree
)
write_unit(unit_4)

#typeV
Cas_typeV <- MAG_565_2 %>% select(genome, typeV)
Cas_typeV <- Cas_typeV %>%
  group_by(genome) %>%
  summarise(
    typeV = ifelse("yes" %in% typeV, "yes", "no")
  ) %>%
  ungroup()
Cas_typeV <- Cas_typeV %>%
  mutate(typeV = case_when(
    typeV == "yes" ~ 1,
    typeV == "no" ~ 0,
    TRUE ~ as.numeric(NA)  # 处理其他情况（如果有的话）
  ))


unit_4 <- create_unit(data = Cas_typeV,
                      key = 'cas_typeV',
                      type = 'DATASET_EXTERNALSHAPE',
                      tree = tree
)
write_unit(unit_4)


#
CRISPR_based <- MAG_565_2 %>% select(genome, linkage_method1)
CRISPR_based <- CRISPR_based %>%
  group_by(genome) %>%
  summarise(
    linkage_method1 = ifelse("CRISPR based" %in% linkage_method1, "yes", "no")
  ) %>%
  ungroup()
CRISPR_based <- CRISPR_based %>%
  mutate(linkage_method1 = case_when(
    linkage_method1 == "yes" ~ 1,
    linkage_method1 == "no" ~ 0,
    TRUE ~ as.numeric(NA)  # 处理其他情况（如果有的话）
  ))

unit_4 <- create_unit(data = CRISPR_based,
                      key = 'CRISPR_based',
                      type = 'DATASET_EXTERNALSHAPE',
                      tree = tree
)
write_unit(unit_4)

#
provirus_based <- MAG_565_2 %>% select(genome, linkage_method2)
provirus_based <- provirus_based %>%
  group_by(genome) %>%
  summarise(
    linkage_method2 = ifelse("Provirus based" %in% linkage_method2, "yes", "no")
  ) %>%
  ungroup()
provirus_based <- provirus_based %>%
  mutate(linkage_method2 = case_when(
    linkage_method2 == "yes" ~ 1,
    linkage_method2 == "no" ~ 0,
    TRUE ~ as.numeric(NA)  # 处理其他情况（如果有的话）
  ))

unit_4 <- create_unit(data = provirus_based,
                      key = 'provirus_based',
                      type = 'DATASET_EXTERNALSHAPE',
                      tree = tree
)
write_unit(unit_4)




#统计
class1 <- Cas_typeI
class1$n <- Cas_typeI$typeI + Cas_typeIII$typeIII
sum(class1$n >0)
class1 <- class1 %>% 
  left_join(df_phycolor,by = c("genome" = "genome")) 
class1 <- class1[class1$n>0,]
table(class1$class)

class2 <- Cas_typeII
class2$n <- Cas_typeII$typeII + Cas_typeV$typeV
sum(class2$n >0)
class2 <- class2 %>% 
  left_join(df_phycolor,by = c("genome" = "genome")) 
class2 <- class2[class2$n>0,]
table(class2$class)



save.image("./1_MAGs_itol_anno.RData")
load("./1_MAGs_itol_anno.RData")

