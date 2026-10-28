# =============================================================================
# Predatory-bacterial-virome -- marker phylogeny / Type V effector (Cas12m) phylogeny
# =============================================================================
# iTOL annotation of the collected Cas12m proteins and TnpB nucleases
# (taxonomy and protein features).
#
# Manuscript      : Fig. 6D
# Main outputs    : iTOL annotation files
# Original script : 1mobilome_predatory_bacteria/9special_protein/cas12m/cas12_phylo/1_Cas12m_itol_anno.R
# -----------------------------------------------------------------------------
# Input data and output files are resolved from GTR_DATA (current analysis root),
# GTR_WORK (earlier working copy) and GTR_DEF (shared defensome data). All three
# default to the current working directory; set them instead of editing files:
#
#     Sys.setenv(GTR_DATA = "/path/to/06_predatory_bacteria_virome")
#
# =============================================================================
GTR_DATA <- normalizePath(Sys.getenv("GTR_DATA", unset = "."), mustWork = FALSE)
GTR_WORK <- normalizePath(Sys.getenv("GTR_WORK", unset = "."), mustWork = FALSE)


setwd(file.path(GTR_WORK, "1mobilome_predatory_bacteria/9special_protein/cas12m/cas12_phylo"))

library(dplyr)

Cas12m_data <- readxl::read_excel("collected_cas12.xlsx")

#itol注释####
library(itol.toolkit)
library(data.table) # 用于元数据文件读取
library(ape) # 用于树文件读写
library(stringr) # 用于字符串处理
library(dplyr)

# 假设您的树存储在一个名为"tree"的变量中
tree <- "all_typeV_effector2_aln.trim.afa.treefile"
phylo <- ape::read.tree(tree)
hub <- create_hub(tree)

#按门设置颜色
library(dplyr)
df_phycolor <- Cas12m_data %>% select(SeqID, Type)
df_phycolor <- df_phycolor %>% distinct(SeqID,.keep_all = TRUE)
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
  "#f7b6d2",  # 浅粉色
  "#dbdb8d"  # 浅橄榄色
 # "#9edae5",  # 浅青色
 # "#393b79",  # 深蓝色
 # "#637939"   # 深绿色
)
element_counts <- df_phycolor %>%
  count(Type)

element_count_sorted <- element_counts[order(-element_counts$n), ]
element_count_sorted$color <- c(colorss )

#
df_phycolor <- merge(df_phycolor,element_count_sorted,by.x = "Type",by.y = "Type")
df_phycolor <- df_phycolor[, c("SeqID", "Type","color")]

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


unit_21 <- create_unit(data = df_phycolor,
                       key = "phylo_color_strip",
                       type = "DATASET_COLORSTRIP",
                       tree = phylo)
write_unit(unit_21)


df_data <- Cas12m_data[,c(1,2)]
unit_30 <- create_unit(data = df_data,
                       key = "phylo_symbol_1",
                       type = "DATASET_SYMBOL",
                       position = 1,
                       tree = tree)
write_unit(unit_30)


