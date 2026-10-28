# =============================================================================
# Predatory-bacterial-virome -- marker phylogeny / Marker-gene phylogeny annotation
# =============================================================================
# Annotate the 39-phage ViPTree maximum-likelihood tree: gene maps aligned
# to the tips plus the marker-gene presence strips.
#
# Manuscript      : Fig. S4
# Main outputs    : 1_phage_gene_annotation.RData, tree annotation figures
# Original script : 09_special_phage_marker_phylogeny/1_phage_gene_annotation.R
# -----------------------------------------------------------------------------
# Input data and output files are resolved from GTR_DATA (current analysis root),
# GTR_WORK (earlier working copy) and GTR_DEF (shared defensome data). All three
# default to the current working directory; set them instead of editing files:
#
#     Sys.setenv(GTR_DATA = "/path/to/06_predatory_bacteria_virome")
#
# =============================================================================
GTR_DATA <- normalizePath(Sys.getenv("GTR_DATA", unset = "."), mustWork = FALSE)


setwd(file.path(GTR_DATA, "09_special_phage_marker_phylogeny"))

all_39phage_gff <- read.delim2("viptree_selected_39_phage_clean_treated.gff",
                                header = F, sep = '\t')

library(dplyr)
library(stringr)
all_39phage_gff <- all_39phage_gff %>%
  mutate(
    # 提取ID=后面的内容，然后取最后一个_后的数字
    gene_num = str_extract(V9, "ID=[^;]+") %>% 
      str_remove("ID=") %>% 
      str_extract("[^_]+$"),
    hit_id = paste0(V1, "_", gene_num)
  ) %>%
  select(-gene_num)  # 移除临时列
all_39phage_gff <- all_39phage_gff[,c(1,3,4,5,7,10)]

colnames(all_39phage_gff) <- c("contig","type","start","stop","strand",
                           "name")

#type和name可以改为需要的内容（type决定颜色，name决定标签）

#write.table(all_phage_gff,"all_phage_gff.tsv",sep = "\t",quote = F)


primary_anno <- read.delim2("39phage_annotation_summary.tsv")

all_39phage_gff$seq <- sprintf("seq%d", 1:nrow(all_39phage_gff))


all_39phage_gff <- all_39phage_gff %>%
  left_join(primary_anno, by = "seq")

write.table(all_39phage_gff,"all_39phage_gff_anno.tsv",sep = "\t",quote = F)
save.image("./1_phage_gene_annotation.RData")
################################################################################

#major capsid protein

# 忽略大小写进行匹配（推荐，因为注释中可能出现"Major capsid"等变体）
filtered_39phage_MCP <- 
  subset(all_39phage_gff, 
         grepl("major capsid", detail_anno, ignore.case = TRUE) &
         !grepl("head decoration", detail_anno, ignore.case = TRUE) &
         !grepl("head maturation protease", detail_anno, ignore.case = TRUE) &
         !grepl("5WK1_M", detail_anno, ignore.case = TRUE)
           )

filtered_39phage_TerL <- 
  subset(all_39phage_gff, 
         grepl("terminase large|large terminase|large subunit", 
               detail_anno, ignore.case = TRUE))

filtered_39phage_portal <- 
  subset(all_39phage_gff, 
         grepl("portal protein|portal", detail_anno, ignore.case = TRUE) &
           !grepl("Tail terminator|Tail tubular|Gene 7|Major tail protein|Minor capsid protein", detail_anno, ignore.case = TRUE)
         )

save.image("./1_phage_gene_annotation.RData")


write.table(filtered_39phage_MCP, 
            "filtered_39phage_MCP.tsv", 
            sep = "\t",           # TSV使用制表符分隔
            row.names = FALSE,    # 不写入行名
            col.names = TRUE,     # 写入列名
            quote = FALSE        # 不添加引号
            )
write.table(filtered_39phage_portal, 
            "filtered_39phage_portal.tsv", 
            sep = "\t",           # TSV使用制表符分隔
            row.names = FALSE,    # 不写入行名
            col.names = TRUE,     # 写入列名
            quote = FALSE        # 不添加引号
            )
write.table(filtered_39phage_TerL, 
            "filtered_39phage_TerL.tsv", 
            sep = "\t",           # TSV使用制表符分隔
            row.names = FALSE,    # 不写入行名
            col.names = TRUE,     # 写入列名
            quote = FALSE        # 不添加引号
            )

selected_7_phage <- read.delim2("selected_phage_anno_withtype.tsv")
selected_7_phage_MCP <- 
  subset(selected_7_phage, 
         grepl("major capsid", detail_anno, ignore.case = TRUE) &
           !grepl("head decoration", detail_anno, ignore.case = TRUE) &
           !grepl("head maturation protease", detail_anno, ignore.case = TRUE) &
           !grepl("5WK1_M", detail_anno, ignore.case = TRUE)
        )
write.table(selected_7_phage_MCP, 
            "selected_7_phage_MCP.tsv", 
            sep = "\t",           # TSV使用制表符分隔
            row.names = FALSE,    # 不写入行名
            col.names = TRUE,     # 写入列名
            quote = FALSE        # 不添加引号
        )
selected_7_phage_TerL <- 
  subset(selected_7_phage, 
         grepl("terminase large|large terminase|large subunit", 
               detail_anno, ignore.case = TRUE))
write.table(selected_7_phage_TerL, 
            "selected_7_phage_TerL.tsv", 
            sep = "\t",           # TSV使用制表符分隔
            row.names = FALSE,    # 不写入行名
            col.names = TRUE,     # 写入列名
            quote = FALSE        # 不添加引号
        )
selected_7_phage_portal <- 
  subset(selected_7_phage, 
         grepl("portal protein|portal", detail_anno, ignore.case = TRUE) &
           !grepl("Tail terminator|Tail tubular|Gene 7|Major tail protein|Minor capsid protein", detail_anno, ignore.case = TRUE)
  )

write.table(selected_7_phage_portal, 
            "selected_7_phage_portal.tsv", 
            sep = "\t",           # TSV使用制表符分隔
            row.names = FALSE,    # 不写入行名
            col.names = TRUE,     # 写入列名
            quote = FALSE        # 不添加引号
        )

save.image("./1_phage_gene_annotation.RData")


################################################################################
all_phage_gff <- read.delim2("ESKAPE_complete_phage_treated.gff",
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

#write.table(all_phage_gff,"all_phage_gff.tsv",sep = "\t",quote = F)


primary_anno <- read.delim2("phage_annotation_summary.tsv")

all_phage_gff$seq <- sprintf("seq%d", 1:nrow(all_phage_gff))


all_phage_gff <- all_phage_gff %>%
  left_join(primary_anno, by = "seq")

write.table(all_phage_gff,"all_phage_gff_anno.tsv",sep = "\t",quote = F)
save.image("./1_phage_gene_annotation.RData")
#######################

#major capsid protein

# 忽略大小写进行匹配（推荐，因为注释中可能出现"Major capsid"等变体）
filtered_phage_MCP <- 
  subset(all_phage_gff, 
         grepl("major capsid", detail_anno, ignore.case = TRUE) &
           !grepl("head decoration", detail_anno, ignore.case = TRUE) &
           !grepl("head maturation protease", detail_anno, ignore.case = TRUE) &
           !grepl("5WK1_M", detail_anno, ignore.case = TRUE)
  )

filtered_phage_TerL <- 
  subset(all_phage_gff, 
         grepl("terminase large|large terminase|large subunit", 
               detail_anno, ignore.case = TRUE))

filtered_phage_portal <- 
  subset(all_phage_gff, 
         grepl("portal protein|portal", detail_anno, ignore.case = TRUE) &
           !grepl("Tail terminator|Tail tubular|Gene 7|Major tail protein|Minor capsid protein", detail_anno, ignore.case = TRUE)
  )





save.image("./1_phage_gene_annotation.RData")

ESKAPE_phage_metadata <- read.delim2("ESKAPE_all_phage_info_filter.tsv")

filtered_phage_MCP <- filtered_phage_MCP %>%
  left_join(ESKAPE_phage_metadata, by = c("contig"="ID"))

filtered_phage_portal <- filtered_phage_portal %>%
  left_join(ESKAPE_phage_metadata, by = c("contig"="ID"))

filtered_phage_TerL <- filtered_phage_TerL %>%
  left_join(ESKAPE_phage_metadata, by = c("contig"="ID"))

write.table(filtered_phage_MCP, 
            "filtered_phage_MCP.tsv", 
            sep = "\t",           # TSV使用制表符分隔
            row.names = FALSE,    # 不写入行名
            col.names = TRUE,     # 写入列名
            quote = FALSE        # 不添加引号
)
write.table(filtered_phage_portal, 
            "filtered_phage_portal.tsv", 
            sep = "\t",           # TSV使用制表符分隔
            row.names = FALSE,    # 不写入行名
            col.names = TRUE,     # 写入列名
            quote = FALSE        # 不添加引号
)
write.table(filtered_phage_TerL, 
            "filtered_phage_TerL.tsv", 
            sep = "\t",           # TSV使用制表符分隔
            row.names = FALSE,    # 不写入行名
            col.names = TRUE,     # 写入列名
            quote = FALSE        # 不添加引号
)

